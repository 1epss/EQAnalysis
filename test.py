import numpy as np
import pandas as pd
from obspy import *
import os

# Write staion.dat for ph2dt
def write_station_file(instrument, filename):
    csv = pd.read_csv(instrument)
    with open(filename, 'a') as f:
        for idx, i in csv.iterrows():
            line = "{:s} {:.4f} {:.4f}\n".format(i['STA'], i['LAT'], i['LON']) 
            f.write(line)
    return f

# Write phase.dat for ph2dt
def write_phase_file(event, waveform_directory, traveltime, filename):
    csv = pd.read_csv(event)
    eventlist = os.listdir(waveform_directory)
    
    with open(filename, 'a') as f:
        for idx, i in csv.iterrows():
            travel_line = "{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.2f} {:.4f} {:.4f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:d}\n".format(
                i['YR'], i['MO'], i['DY'], i['HR'], i['MI'], i['SC'], i['LAT'], i['LON'], i['DEP'], i['MAG'], i['EH1'], i['EZ'], i['RMS'], idx + 1) 
            f.write(travel_line)
            
        for event in eventlist:
            waveformlist = os.listdir(f'{waveform_directory}/{event}')
            for waveform in waveformlist:
                if not waveform.endswith('.sac'):
                    continue
                st = read(f'{waveform_directory}/{event}/{waveform}')
                sta = st[0].stats.station
                #tt = 
                wght = 1.0
                nobs_line = "{:>7s} {:.0f} {:.0f} {:.0f}\n".format(sta, i['MO'], i['DY'], i['HR']) 
            
    return f

#write_station_file(instrument = './korea_instrument_2024_Ahn.csv', filename = 'station.dat')
#write_phase_file(event = './event.csv', waveform_directory = '../cc_final', traveltime = 'kimetalhypoel.out', filename = 'phase.dat')