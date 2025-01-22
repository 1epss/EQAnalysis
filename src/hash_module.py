import numpy as np
import pandas as pd
from obspy import *
from obspy.imaging.beachball import beachball
from folium import plugins
import os

# degree to dms
def degree_to_dms(deg):
    degree = int(deg)
    minute = int((deg - degree) * 60)
    second = (deg - degree - minute/60) * 3600
    second = int(round(second, 0))
    return degree, minute, second


# Write station.hash file for hash
def hash_write_station(instrument, filename):
    csv = pd.read_csv(instrument)
    csv['ELEV'] = csv['ELEV'].astype(int)
    with open(filename, 'a') as f:
        for _, row in csv.iterrows():
            if len(row['STA']) != 5:
                line = "{:<4s}{:>4s}{:>42.5f}{:>11.5f}{:>6d}{:>25}\n".format(row['STA'], row['CHANNEL'], row['LAT'], row['LON'], row['ELEV'], row['NET'])
                f.write(line)
    with open(filename, 'r') as f:
        lines = f.readlines()
    unique_lines = list(dict.fromkeys(lines))
    sort_lines = sorted(unique_lines)
    with open(filename, 'w') as f:
        f.writelines(sort_lines)
    return f


# Write saclst.hash file for phase.hash
def hash_write_saclst(event_directory, filename):
    event_list = sorted(os.listdir(event_directory))
    with open(filename, 'a') as f:
        for event in event_list:
            df = pd.DataFrame(columns=[event, 'DIST', 'P', 'S'])            
            station_list = sorted(os.listdir(os.path.join(event_directory, event)))
            for station in station_list:
                st = read(os.path.join(event_directory, event, station))
                try:
                    if 'Z' in st[0].stats.channel:
                        new_row = {event: station, 'DIST': st[0].stats.sac.dist, 'P': st[0].stats.sac.a, 'S': st[0].stats.sac.t0}
                        df = pd.concat([df, pd.DataFrame(new_row, index=[0])], ignore_index=True)
                except:
                    continue
            df['STATION'] = df[event].apply(lambda x : x.split('.')[0])
            df['CHANNEL'] = df[event].apply(lambda x : x.split('.')[2])
            priority_order = {'HHZ': 0, 'ELZ': 1, 'HGZ': 2}
            df['priority'] = df['CHANNEL'].map(priority_order)
            df = df.sort_values(by=['STATION', 'priority']).drop_duplicates(subset=['STATION'], keep='first').drop(columns=['priority', 'CHANNEL', 'STATION']).reset_index(drop=True)
            header_line = f"{event} DIST P S\n"
            f.write(header_line)
            for _, row in df.iterrows():
                phase_line = "{:<15s}{:>13.4f}{:>12.4f}{:>12.4f}\n".format(row[event], row['DIST'], row['P'], row['S'])
                f.write(phase_line)


# Write phase.hash file for hash
def hash_write_phase(event_csv, saclst_file, event_directory, filename):
    csv = pd.read_csv(event_csv)
    saclst = pd.read_csv(saclst_file, sep='\s+', header = None)
    event_list = sorted(os.listdir(event_directory))

    for idx, event in csv.iterrows():
        with open(filename, 'a') as f:
            event_num = 100001 + idx
            org = UTCDateTime.strptime(f'{int(event.YR)}.{int(event.MO)}.{int(event.DY)} {int(event.HR)}:{int(event.MI)}:{round(event.SC, 4)}', '%Y.%m.%d %H:%M:%S.%f')
            org_format = UTCDateTime.strftime(org, '%Y%m%d%H%M%S.%f')[:-4]
            lat1, lat2, lat3 = degree_to_dms(event.LAT)
            lon1, lon2, lon3 = degree_to_dms(event.LON)
            depth, hor_uncer, ver_uncer, mag = event.DEP, event.EH1, event.EZ, event.MAG
            event_line = "{:<17s}{:<2d}N{:<2d}.{:<2d}{:<3d}E{:<2d}.{:<2d} {:>.2f}{:>54.2f}{:>6.2f}{:>44.2f}{:>22d}\n".format(org_format, lat1, lat2, lat3, lon1, lon2, lon3, depth, hor_uncer, ver_uncer, mag, event_num)
            f.write(event_line)

            sections = []
            current_section = None
            for _, row in saclst.iterrows():
                if row[1] == 'DIST':
                    if current_section is not None:
                        sections.append(current_section)
                    current_section = pd.DataFrame(columns = saclst.columns)
                if current_section is not None:
                    current_section = pd.concat([current_section, pd.DataFrame([row], columns=saclst.columns)], ignore_index=True)
            if current_section is not None:
                sections.append(current_section)

            for idx3, section in sections[idx].iterrows():
                if idx3 == 0:
                    continue
                st = read(f'{event_directory}/{event_list[idx]}/{section[0]}')
                sta, net, chan = st[0].stats.station, st[0].stats.network, st[0].stats.channel
                try:
                    ka = st[0].stats.sac.ka
                    kt0 = st[0].stats.sac.kt0 # S파 피킹 안된 파형들 제거하는 목적
                    onset = ka[0]
                    polarity = ka[2]
                    polarity_line = "{:<4s}{:>3s}{:>5s}{:>2s}{:>2s}\n".format(sta, net, chan, onset, polarity)
                    f.write(polarity_line)
                except:
                    continue
            conclude_line = "{:>72d}\n".format(event_num)
            f.write (conclude_line)
            

# Write amp.hash file for hash
def hash_write_amp(event_csv, saclst_file, event_directory, filename, method = 'default'):
    csv = pd.read_csv(event_csv)
    event_list = sorted(os.listdir(event_directory))
    saclst = pd.read_csv(saclst_file, sep='\s+', header = None)
    sections = []
    current_section = None
    for _, row in saclst.iterrows():
        if row[1] == 'DIST':
            if current_section is not None:
                sections.append(current_section)
            current_section = pd.DataFrame(columns = saclst.columns)
        if current_section is not None:
            current_section = pd.concat([current_section, pd.DataFrame([row], columns=saclst.columns)], ignore_index=True)
    if current_section is not None:
        sections.append(current_section)

    with open(filename, 'a') as f:
        event_num = 100001
        for idx, row in csv.iterrows():
            org = UTCDateTime(int(row.YR), int(row.MO), int(row.DY), int(row.HR), int(row.MI), float(row.SC))
            observed_num = len(sections[idx])
            event_line = "{:<6d}{:>24d}\n".format(event_num, observed_num)
            f.write(event_line)
            for idx2, i in sections[idx].iterrows():
                if idx2 == 0:
                    continue
                st = read(f'{event_directory}/{event_list[idx]}/{i[0]}')
                sta, net, chan, starttime = st[0].stats.station, st[0].stats.network, st[0].stats.channel, st[0].stats.starttime

                if method == 'default':     
                    try:
                        a = int(np.round(org + st[0].stats.sac.a - starttime, 2) * 100)
                        t0 = int(np.round(org + st[0].stats.sac.t0 - starttime, 2) * 100)
                        a_amp = max(abs(st[0].data[a - 30 : a + 30]))
                        t0_amp = max(abs(st[0].data[t0 - 30 : t0 + 30]))
                        a_noise_amp = st[0].data[a - 20]
                        t0_noise_amp = st[0].data[t0 - 20] 
                        amplitude_line = "{:<4s}{:>4s}{:>3s}{:>27.3f}{:>11.3f}{:>11.3f}{:>11.3f}\n".format(sta, chan, net, a_noise_amp, a_noise_amp, a_amp, t0_amp)   
                        f.write(amplitude_line)
                    except:
                        continue

                elif method == 'vector_sum':
                    st2 = read(f'{event_directory}/{event_list[idx]}/{i[0][:-5]}N.sac')     
                    st3 = read(f'{event_directory}/{event_list[idx]}/{i[0][:-5]}E.sac')                
                    try:
                        a = int(np.round(org + st[0].stats.sac.a - starttime, 2) * 100)
                        t0 = int(np.round(org + st[0].stats.sac.t0 - starttime, 2) * 100)
                        a_amp1 = max(abs(st[0].data[a - 30 : a + 30]))
                        a_amp2 = max(abs(st2[0].data[a - 30 : a + 30]))
                        a_amp3 = max(abs(st3[0].data[a - 30 : a + 30]))
                        a_amp = np.sqrt(a_amp1 ** 2 + a_amp2 ** 2 + a_amp3 ** 2)
                        t0_amp1 = max(abs(st[0].data[t0 - 30 : t0 + 30]))
                        t0_amp2 = max(abs(st2[0].data[t0 - 30 : t0 + 30]))
                        t0_amp3 = max(abs(st3[0].data[t0 - 30 : t0 + 30]))
                        t0_amp = np.sqrt(t0_amp1 ** 2 + t0_amp2 ** 2 + t0_amp3 ** 2)
                        a_noise_amp1 = np.median(st[0].data[a - 20 : a - 1])
                        a_noise_amp2 = np.median(st2[0].data[a - 20: a - 1])
                        a_noise_amp3 = np.median(st3[0].data[a - 20: a - 1])
                        a_noise_amp = np.sqrt(a_noise_amp1 ** 2 + a_noise_amp2 ** 2 + a_noise_amp3 ** 2)
                        t0_noise_amp1 = np.median(st[0].data[t0 - 20 : t0 - 1]) 
                        t0_noise_amp2 = np.median(st2[0].data[t0 - 20 : t0 - 1]) 
                        t0_noise_amp3 = np.median(st3[0].data[t0 - 20 : t0 - 1]) 
                        t0_noise_amp = np.sqrt(t0_noise_amp1 ** 2 + t0_noise_amp2 ** 2 + t0_noise_amp3 ** 2)
                        amplitude_line = "{:<4s}{:>4s}{:>3s}{:>27.3f}{:>11.3f}{:>11.3f}{:>11.3f}\n".format(sta, chan, net, a_noise_amp, a_noise_amp, a_amp, t0_amp)   
                        f.write(amplitude_line)
                    except:
                        continue
            event_num += 1
    return saclst