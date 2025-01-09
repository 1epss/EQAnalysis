import numpy as np
import pandas as pd
from obspy import *

# Write staion.dat for ph2dt
def write_station_file(instrument, filename):
    csv = pd.read_csv(instrument)
    with open(filename, 'a') as f:
        for idx, i in csv.iterrows():
            line = "{:s} {:.4f} {:.4f}\n".format(i['STA'], i['LAT'], i['LON']) 
            f.write(line)
    return f

# phase.dat 생성
def write_phase_file(event, filename):
    csv = pd.read_csv(event)
    
    with open(filename, 'a') as f:
        for idx, i in csv.iterrows():
            # 연 월 일 시 분 초 위 경 깊 규 평오 직오 RMS ID
            travel_line = "{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.2f} {:.4f} {:.4f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:d}\n".format(
                i['YR'], i['MO'], i['DY'], i['HR'], i['MI'], i['SC'], i['LAT'], i['LON'], i['DEP'], i['MAG'], i['EH1'], i['EZ'], i['RMS'], idx + 1) 
            
            f.write(            travel_line = "{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.2f} {:.4f} {:.4f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:d}\n".format(
)
            for
                nobs_line = 

    return f


# ph2dt.inp 생성




# 실행 코드
#write_station_file(instrument = './korea_instrument_2024_Ahn.csv', filename = 'station.dat')
write_phase_file(event = './event.csv', filename = 'phase.dat')