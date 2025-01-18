from obspy import *
import numpy as np
import pandas as pd
import os


# Write station.dat for ph2dt
def write_station_file(instrument, filename):
    csv = pd.read_csv(instrument).drop_duplicates(subset=['NET', 'STA'], keep='last').sort_values(by=['STA'])
    with open(filename, 'a') as f:
        for idx, i in csv.iterrows():
            line = "{:s} {:.4f} {:.4f}\n".format(i['STA'], i['LAT'], i['LON']) 
            f.write(line)
    return f


# Write traveltime.dat from hypoel.out
def write_traveltime_file(hypoel_out, filename):
    with open(hypoel_out, 'r') as f:
        a = f.readlines()
    
    start_condition = np.array([ 'travel times and delays' in line for line in a ])
    end_condition = np.array([ 'magnitude data' in line for line in a ])
    start = np.where(start_condition == True) 
    end = np.where(end_condition == True)

    for i in range(len(start[0])):
        d = pd.DataFrame(a[start[0][i]:end[0][i]][:-1])
        for idx, j in d.iterrows():
            component = j[0].split()
            if idx == 0:
                continue
            elif idx == 1 :
                df = pd.DataFrame(columns = component)
            else :
                comp_extended = component + [np.nan] * (df.shape[1] - len(component))
                dummy_df = pd.DataFrame([comp_extended], columns=df.columns)
                df = pd.concat([df, dummy_df], ignore_index=True)
                df = df.loc[:, ~df.columns.duplicated()]

        tt_list = []
        with open(filename, 'a') as g:
            for idx2, k in df.iterrows():
                if 'P' in str(k['pha']):
                    tt_list.append((k['stn'].upper(), float(k['ain']), 1.0, 'P'))
                if 'S' in str(k['remk']):  
                    tt_list.append((k['stn'].upper(), float(k['dist']), 1.0, 'S'))    

                nobs_line = "{:<7s} {:>7.3f} {:>.0f} {:s}\n".format(tt_list[-1][0], tt_list[-1][1], tt_list[-1][2], tt_list[-1][3]) 
                g.write(nobs_line)
            g.write("\n")
    return g


# Write phase.dat for ph2dt
def write_phase_file(event, traveltime, filename):
    csv = pd.read_csv(event)

    with open(filename, 'a') as f:
        with open(traveltime, "r") as g:
            line = g.readlines()
            
        blank_line = 0
        for idx, i in csv.iterrows():
            travel_line = "# {:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.2f} {:.4f} {:.4f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:d}\n".format(
                i['YR'], i['MO'], i['DY'], i['HR'], i['MI'], i['SC'], i['LAT'], i['LON'], i['DEP'], i['MAG'], i['EH1'], i['EZ'], i['RMS'], idx + 1) 
            f.write(travel_line)
            
            for idx2, nobs_line in enumerate(line[blank_line:]):
                if nobs_line.strip() == "":
                    blank_line += idx2 + 1
                    break
                f.write(nobs_line)         
    return f


# Compute Cross-correlation and write cc.dat for hypoDD
def write_cc_file(event_directory, station, filename, phase, btime, channel_list, resampling_rate, freqmin, freqmax, minsec, maxsec, threshold):
    
    # Input station.dat from write_station_file module
    station_list = pd.read_csv(station, sep='\s+', header = None)
    station_list.columns = ['STATION', 'LATITUDE', "LONGITUDE"]
    # Input event_directory
    event_list = sorted([event for event in os.listdir(event_directory)])
    acc = channel_list[0][-1]

    # Compute cross-correlation between two different events
    for idx, event in enumerate(event_list):
        for idx2, event2 in enumerate(event_list):
            if idx2 == idx:
                continue
            
            with open(filename, 'a') as f:
                otc = 0.0
                header_line = "# {:d} {:d} {:.1f}\n".format(idx + 1, idx2 + 1, otc) 
                f.write(header_line)
                
                for station in station_list['STATION']:
                    for channel in channel_list[0]:

                        try:
                            slave = read(f'{event_directory}/{event}/{station}.KS.{channel}.sac')
                            master = read(f'{event_directory}/{event2}/{station}.KS.{channel}.sac')

                            if channel == acc and os.path.exists(os.path.join(f'{event_directory}/{event}', f'{station}.KS.HHZ.sac')):
                                continue

                            if channel == acc and os.path.exists(os.path.join(f'{event_directory}/{event}', f'{station}.KS.ELZ.sac')):
                                continue

                            slave = slave.resample(resampling_rate)
                            master = master.resample(resampling_rate)

                            slave = slave.filter('bandpass', freqmin = freqmin, freqmax = freqmax)
                            master = master.filter('bandpass', freqmin = freqmin, freqmax = freqmax)

                            if phase == 'P':
                                time_slave = UTCDateTime(slave[0].stats.starttime + btime + slave[0].stats.sac.a) 
                                time_master = UTCDateTime(master[0].stats.starttime + btime + master[0].stats.sac.a)
                            elif phase == 'S':
                                time_slave = UTCDateTime(slave[0].stats.starttime + 5 + slave[0].stats.sac.t0)
                                time_master = UTCDateTime(master[0].stats.starttime + 5 + master[0].stats.sac.t0)

                            wave_slave = slave[0].slice(time_slave - minsec, time_slave + maxsec).data
                            wave_master = master[0].slice(time_master - minsec, time_master + maxsec).data

                            #numerator / denominator = np.corrcoef
                            numerator = np.sum(wave_slave * wave_master)
                            denominator = np.sqrt(np.sum(wave_slave ** 2) * np.sum(wave_master ** 2))

                            corr = np.correlate(wave_slave, wave_master, 'full') / denominator
                            maxcorr = round(max(corr) ** 2, 3) 
                            lag = np.where(corr == max(corr))
                            # if 100Hz waveform : btime * 100
                            lag_time = (lag[0][0] - btime * 100) / resampling_rate
                            
                            if phase == 'P':
                                obs_line = "{:s} {:.3f} {:.3f} P\n".format(station, lag_time, maxcorr) 
                            elif phase == 'S':
                                obs_line = "{:s} {:.3f} {:.3f} S\n".format(station, lag_time, maxcorr)
                            if maxcorr >= threshold:
                                f.write(obs_line)

                        except:
                            continue