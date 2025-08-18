from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import obspy
import pandas as pd

# =============== 공통 헬퍼 ===============


def _ensure_dir(p: str | Path) -> Path:
    d = Path(p).expanduser().resolve()
    if not d.is_dir():
        raise FileNotFoundError(f"Directory not found: {d}")
    return d


def _list_z_sac_files(dirpath: Path) -> list[Path]:
    files = [
        p
        for p in dirpath.iterdir()
        if p.is_file() and re.search(r"Z\.sac$", p.name, re.IGNORECASE)
    ]
    files.sort(key=lambda x: (len(x.name), x.name))
    return files


def _num_ok(v) -> bool:
    try:
        v = float(v)
        return np.isfinite(v) and v > -1.0e4  # SAC 미설정(-12345) 배제
    except Exception:
        return False


def _str_ok(val) -> bool:
    if val is None:
        return False
    s = str(val).strip()
    return s != "" and s != "-12345"


def _normalize_code(s: str) -> str:
    return str(s).strip().upper()


# =============== saclst 실행/파싱 ===============


def run_saclst_shell(directory: str | Path) -> str:
    directory = _ensure_dir(directory)
    cmd = "saclst kstnm knetwk stla stlo stel f *Z.sac"
    r = subprocess.run(
        ["bash", "-lc", cmd],
        cwd=str(directory),
        capture_output=True,
        text=True,
        check=True,
    )
    return r.stdout


def parse_saclst_stdout(stdout: str) -> pd.DataFrame:
    """
    saclst 출력(한 줄 = 파일 하나)을 파싱.
    형식 예) FNAME kstnm STA knetwk NET stla LAT stlo LON stel ELEV
    """
    rows: list[dict] = []
    for line in stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        kv: dict[str, str] = {}
        i = 0
        while i + 1 < len(parts):
            key = parts[i].lower()
            val = parts[i + 1]
            kv[key] = val
            i += 2
        # saclst는 보통 '... f <filename>' 형태 → f 키에서 파일명 취득
        fname = kv.get("f")
        if not fname:
            # 혹시 포맷이 달라도 마지막 토큰이 *.sac이면 사용
            last = parts[-1]
            fname = last if last.lower().endswith(".sac") else None
        if not fname:
            continue
    
        sta = kv.get("kstnm")
        net = kv.get("knetwk")
        try:
            stla = float(kv.get("stla", "nan"))
            stlo = float(kv.get("stlo", "nan"))
            stel = float(kv.get("stel", "nan"))
        except ValueError:
            continue

        if not (sta and net):
            continue
        rows.append(
            dict(
                NETWORK=str(net).strip(),
                STATION=str(sta).strip(),
                LATITUDE=stla if _num_ok(stla) else np.nan,
                LONGITUDE=stlo if _num_ok(stlo) else np.nan,
                ELEVATION=stel if _num_ok(stel) else np.nan,
                FILE=fname,
            )
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        # 숫자 칼럼 정리
        for c in ["LATITUDE", "LONGITUDE", "ELEVATION"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# =============== station 메타(헤더 폴백) ===============


def _build_station_from_headers(waveform_dir: Path) -> pd.DataFrame:
    rows: list[dict] = []
    for p in _list_z_sac_files(waveform_dir):
        try:
            st = obspy.read(str(p), headonly=True)
            tr = st[0]
            sac = tr.stats.sac
            sta = getattr(sac, "kstnm", None)
            net = getattr(sac, "knetwk", None) or tr.stats.network
            stla = getattr(sac, "stla", None)
            stlo = getattr(sac, "stlo", None)
            stel = getattr(sac, "stel", None)
            if sta and net and _num_ok(stla) and _num_ok(stlo):
                rows.append(
                    dict(
                        NETWORK=str(net).strip(),
                        STATION=str(sta).strip(),
                        LATITUDE=float(stla),
                        LONGITUDE=float(stlo),
                        ELEVATION=float(stel) if _num_ok(stel) else 0.0,  # <- 일관 폴백
                    )
                )
        except Exception as e:
            print(f"[SKIP] header read failed: {p.name} ({e})")
            continue
    df = pd.DataFrame(rows)
    if not df.empty:
        df = (
            df.sort_values(["NETWORK", "STATION"])
            .drop_duplicates(["NETWORK", "STATION"], keep="first")
            .reset_index(drop=True)
        )
    return df


# =============== 대표 파일 선택 ===============


def _choose_unit_files(waveform_dir: Path) -> list[str]:
    """
    우선순위: (has_S, has_P, chan_prio) 내림차순
              chan_prio: HHZ(3) > HGZ(2) > ELZ(1)
    """
    files = _list_z_sac_files(waveform_dir)
    prio = {"HHZ": 3, "HGZ": 2, "ELZ": 1}
    best: dict[str, dict] = {}

    for p in files:
        try:
            st = obspy.read(str(p), headonly=True)
        except Exception as e:
            print(f"[SKIP] read failed: {p.name} ({e})")
            continue

        tr = st[0]
        sac = getattr(tr.stats, "sac", None)
        sta = str(getattr(sac, "kstnm", None) or p.name.split(".")[0]).strip()
        chan = str(
            getattr(tr.stats, "channel", None)
            or (p.name.split(".")[2] if "." in p.name else "")
        ).upper()

        has_p = _pick_ok(getattr(sac, "a", None))
        has_s = _pick_ok(getattr(sac, "t0", None))
        score = (1 if has_s else 0, 1 if has_p else 0, prio.get(chan, 0))

        prev = best.get(sta)
        if (prev is None) or (score > prev["score"]):
            best[sta] = dict(file=p.name, score=score)

    return [best[k]["file"] for k in sorted(best.keys(), key=lambda s: (len(s), s))]


def _pick_ok(val) -> bool:
    return _num_ok(val)


# =============== init ===============


def init(waveform_directory: str | Path) -> tuple[pd.DataFrame, list[str]]:
    """
    station 메타와 대표파일(unit_array) 생성:
    1) saclst → 파싱
    2) station.txt (dir/CWD)
    3) SAC 헤더 폴백
    """
    waveform_dir = _ensure_dir(waveform_directory)
    station: pd.DataFrame | None = None

    if shutil.which("saclst"):
        try:
            df = parse_saclst_stdout(run_saclst_shell(waveform_dir))
            if not df.empty:
                station = df[
                    ["NETWORK", "STATION", "LATITUDE", "LONGITUDE", "ELEVATION"]
                ].copy()
        except Exception as e:
            print(f"[WARN] saclst failed: {e}")

    if station is None or station.empty:
        for path in (waveform_dir / "station.txt", Path("./station.txt")):
            if path.exists():
                try:
                    station = pd.read_csv(
                        path,
                        sep=r"\s+",
                        engine="python",
                        header=None,
                        names=[
                            "NETWORK",
                            "STATION",
                            "LATITUDE",
                            "LONGITUDE",
                            "ELEVATION",
                        ],
                    )
                    break
                except Exception as e:
                    print(f"[WARN] read station.txt failed ({path}): {e}")

    if station is None or station.empty:
        station = _build_station_from_headers(waveform_dir)

    if station is None:
        station = pd.DataFrame(
            columns=["NETWORK", "STATION", "LATITUDE", "LONGITUDE", "ELEVATION"]
        )

    print("=" * 30)
    print("파일을 정렬중입니다... ")
    print("=" * 30)
    unit_array = _choose_unit_files(waveform_dir)
    print("picked : ", unit_array)

    if not station.empty:
        station["STATION"] = station["STATION"].astype(str).str.strip()

    return station, unit_array


# =============== hypo.phase 생성 ===============


def make_hypo(filename: str | Path, wavedata_direction: str | Path) -> None:
    """
    hypo.phase 작성: 관측소별 대표 Z성분 파일에서
    P(필수)/S(있으면) 도달시각을 포맷팅하여 라인 출력.
    """
    waveform_dir = _ensure_dir(wavedata_direction)
    out_path = Path(filename).expanduser().resolve()
    _, unit_array = init(waveform_dir)

    def _ptime_from_epoch(ts: float) -> tuple[str, float]:
        t = obspy.UTCDateTime(ts)
        head = t.strftime("%y%m%d%H%M")  # 분까지
        sec = t.second + t.microsecond / 1e6
        return f"{head}{sec:05.2f}", sec

    with out_path.open("w") as f:
        for name in unit_array:
            sac_path = waveform_dir / name
            try:
                st = obspy.read(str(sac_path), headonly=True)
            except Exception as e:
                print(f"[SKIP] read failed: {sac_path.name} ({e})")
                continue

            sac = st[0].stats.sac
            tnm = str(getattr(sac, "kstnm", "")).strip()
            ka = getattr(sac, "ka", None)

            # P: a 유효 && ka 유효
            if (not _pick_ok(getattr(sac, "a", None))) or (not _str_ok(ka)):
                print(
                    f"[SKIP] invalid P pick/label: {sac_path.name} a={getattr(sac,'a',None)} ka={ka!r}"
                )
                continue

            # 기준 epoch: 헤더 기준시각 - b
            req = ("nzyear", "nzjday", "nzhour", "nzmin", "nzsec", "nzmsec", "b")
            if any(not hasattr(sac, k) for k in req) or (
                not _num_ok(getattr(sac, "b", None))
            ):
                print(f"[SKIP] missing required headers: {sac_path.name}")
                continue

            base_ts = obspy.UTCDateTime(
                year=int(sac.nzyear),
                julday=int(sac.nzjday),
                hour=int(sac.nzhour),
                minute=int(sac.nzmin),
                second=int(sac.nzsec),
                microsecond=int(sac.nzmsec) * 1000,
            ).timestamp - float(sac.b)

            p_ts = base_ts + float(sac.a)
            ptime, p_sec = _ptime_from_epoch(p_ts)
            ka = str(ka).strip()

            # S가 있으면 stime/kt0 출력
            if _pick_ok(getattr(sac, "t0", None)):
                s_ts = base_ts + float(sac.t0)
                pstim = s_ts - p_ts
                stime = f"{(p_sec + pstim):05.2f}"
                kt0 = (
                    str(getattr(sac, "kt0", "")).strip()
                    if _str_ok(getattr(sac, "kt0", None))
                    else "S"
                )
                line = f"{tnm:<4}{ka:<5}{ptime:<22}{stime}{kt0}"
            else:
                line = f"{tnm:<4}{ka:<5}{ptime}"

            f.write(line + "\n")


# =============== station.sta 생성 ===============


def _deg_min_str(lat_deg: float, lon_deg: float) -> tuple[int, str, str, int, str, str]:
    ns = "N" if lat_deg >= 0 else "S"
    ew = "E" if lon_deg >= 0 else "W"
    alat, alon = abs(float(lat_deg)), abs(float(lon_deg))
    ilat, ilon = int(alat), int(alon)
    mlat = round((alat - ilat) * 60.0, 2)
    mlon = round((alon - ilon) * 60.0, 2)
    if mlat >= 60.0:
        ilat += 1
        mlat = 0.0
    if mlon >= 60.0:
        ilon += 1
        mlon = 0.0
    return ilat, ns, f"{mlat:05.2f}", ilon, ew, f"{mlon:05.2f}"


def make_station(filename: str | Path, wavedata_direction: str | Path) -> None:
    """
    station.sta 작성:
    NETWORK STATION LAT/LON(도·분) ELEV 형식으로 두 줄(두 번째 줄은 '*') 출력.
    - station DF에 관측소가 없으면 해당 SAC 헤더(STLA/STLO/STEL)로 폴백.
    """
    waveform_dir = _ensure_dir(wavedata_direction)
    out_path = Path(filename).expanduser().resolve()
    station, unit_array = init(waveform_dir)

    # DF 사전화
    meta: dict[str, tuple[float, float, float]] = {}
    if station is not None and not station.empty:
        for _, row in station.iterrows():
            try:
                sta = _normalize_code(row["STATION"])
                la, lo, el = (
                    float(row["LATITUDE"]),
                    float(row["LONGITUDE"]),
                    float(row["ELEVATION"]),
                )
                meta[sta] = (la, lo, el)
            except Exception:
                continue

    with out_path.open("w") as f:
        for name in unit_array:
            sac_path = waveform_dir / name
            try:
                st = obspy.read(str(sac_path), headonly=True)
            except Exception as e:
                print(f"[SKIP] read failed: {sac_path.name} ({e})")
                continue

            sac = st[0].stats.sac
            tnm = _normalize_code(getattr(sac, "kstnm", ""))

            la = lo = el = None
            if tnm in meta:
                la, lo, el = meta[tnm]

            if not (_num_ok(la) and _num_ok(lo)):
                la_hdr, lo_hdr, el_hdr = (
                    getattr(sac, "stla", None),
                    getattr(sac, "stlo", None),
                    getattr(sac, "stel", None),
                )
                if _num_ok(la_hdr) and _num_ok(lo_hdr):
                    la = float(la_hdr)
                    lo = float(lo_hdr)
                    el = float(el_hdr) if _num_ok(el_hdr) else 0.0
                else:
                    print(f"[SKIP] station meta not found (DF/HEADER): {tnm}")
                    continue

            ilat, ns, lat_s, ilon, ew, lon_s = _deg_min_str(la, lo)
            # el이 NaN/None이면 0으로 처리
            try:
                elev = int(round(float(el))) if _num_ok(el) else 0
            except Exception:
                elev = 0            
            line1 = f"{tnm:>4}{ilat:>2}{ns}{lat_s}{ilon:>4}{ew}{lon_s}{elev:>5}"
            line2 = f"{tnm:>4}*"
            f.write(line1 + "\n" + line2 + "\n")


# =============== 예시 실행 ===============
if __name__ == "__main__":
    wf_dir = Path("~/mingyu_JS/cc_final/200427110710").expanduser()
    try:
        make_hypo("hypo.phase", wf_dir)
        make_station("station.sta", wf_dir)
        print("Done.")
    except Exception as e:
        print(f"Error: {e}")
