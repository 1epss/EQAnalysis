import pandas as pd
import numpy as np
from obspy import *
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings(action='ignore')

# Plot seismic phases per stations    
def plot_phase(instrument, event_directory, channel, method = 'both'):
    
    # Get all available stations from instrument file
    instrument = pd.read_csv(instrument)
    station_list = instrument['STA'].unique() 

    # Make sure there's no addional directory in event_directory
    event_list = sorted([event for event in os.listdir(event_directory)])

    # Compare the station list in the station_file with SAC files in the event_directory and extract distances only for matching files
    for station in station_list:
        matched_list = [os.path.join(event_directory, event, f'{station}.KS.{channel}.sac') for event in event_list if os.path.isfile(os.path.join(event_directory, event, f'{station}.KS.{channel}.sac'))]
        dist_matched_list = [matched for matched in matched_list if read(matched)[0].stats.sac.dist <= 100]
        dist_list = [read(matched)[0].stats.sac.dist for matched in dist_matched_list]
        df_filename = pd.DataFrame(data = dist_matched_list, columns = ['Filename'])
        df_dist = pd.DataFrame(data = dist_list, columns = ['Dist'])
        df = pd.concat([df_filename, df_dist], axis = 1).reset_index()
        
        plt.figure(figsize = (10,10))
        idxx = 0
        for idx, row in df.iterrows():
            st = read(row['Filename'])
            if method == 'both': 
                try : 
                    tr = st[0].trim(UTCDateTime(st[0].stats.starttime + st[0].stats.sac.a), UTCDateTime(st[0].stats.starttime + st[0].stats.sac.a + 20))
                    p_trv = tr.stats.sac.a
                    plt.scatter(500, idxx, facecolor = 'None' , color = 'blue', s = 1000, zorder = 10, lw = 1)
                    plt.plot(tr.data / np.max(np.abs(tr.data))*3 + idxx, 'black')
                    s_trv = tr.stats.sac.get('t0', None)
                    if s_trv is not None:
                        plt.scatter((s_trv - p_trv) * 100 + 500, idxx, facecolor='None', color='red', s=1000, zorder=10, lw=1)        
                except :
                    pass
                plt.title(st[0].stats.network + '.' + st[0].stats.station + '.' + st[0].stats.channel, fontsize=20)                
                plt.xlim(0,2000)
                plt.ylim(-65,5)
                plt.xticks([0,250,500,750,1000,1250,1500,1750,2000],[-5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15], fontsize=12)
                plt.vlines(x=500, ymin = -65, ymax = 5, colors='b', linestyles='--', linewidth = 0.5 )
                plt.xlabel('Time Since P Arrival (s)', fontsize = 15)
                plt.yticks([])
                plt.text(-405, idxx, "          Event  " + str(idx + 1) + "\n" + st[0].stats.starttime.strftime('%Y.%m.%d %H:%M:%S'), fontsize=10)
                plt.text(2020, idxx, "M$_L$ " + str(st[0].stats.sac.mag), fontsize=10)
                plt.savefig(f'{station}.png')

            elif method == 'P':  
                try : 
                    tr = st[0].trim(UTCDateTime(st[0].stats.starttime + st[0].stats.sac.a + 4.5), UTCDateTime(st[0].stats.starttime + st[0].stats.sac.a + 6))
                    tr.filter('highpass', freq = 2)
                    p_trv = tr.stats.sac.a
                    plt.plot(tr.data / np.max(np.abs(tr.data))*3 + idxx, 'black')
                    # Plot additional picks if available
                    try : 
                        t3_trv = tr.stats.sac.t3 - tr.stats.sac.a
                        plt.vlines(x = 50 + t3_trv*100, ymin = idxx-3, ymax = idxx +3, colors = 'magenta', linestyles='-', linewidth = 0.5)
                    except :
                        pass
                except :
                    pass
                plt.title(st[0].stats.network + '.' + st[0].stats.station + '.' + st[0].stats.channel, fontsize=20)
                plt.xlim(0,150)
                plt.ylim(-65,5)
                plt.xticks([0,25,50,75,100,125,150],[-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], fontsize=12)
                plt.vlines(x=50, ymin = -65, ymax = 5, colors='b', linestyles='--', linewidth = 0.5)
                plt.xlabel('Time Since P Arrival (s)', fontsize = 15)
                plt.yticks([])
                plt.text(-30, idxx, "          Event  " + str(idx + 1) + "\n" + st[0].stats.starttime.strftime('%Y.%m.%d %H:%M:%S'), fontsize=10)
                plt.text(152, idxx, "M$_L$ " + str(st[0].stats.sac.mag), fontsize=10)
                plt.savefig(f'{station}.png')

            elif method == 'S' :
                try : 
                    tr = st[0].trim(UTCDateTime(st[0].stats.starttime + st[0].stats.sac.t0 +4.5), UTCDateTime(st[0].stats.starttime + st[0].stats.sac.t0 + 6))
                    tr.filter('highpass', freq = 2)
                    s_trv = tr.stats.sac.t0
                    plt.plot(tr.data/np.max(np.abs(tr.data))*3 + idxx, 'black')
                    # Plot additional picks if available
                    try : 
                        t4_trv = tr.stats.sac.t4 - tr.stats.sac.t0
                        plt.vlines(x = 50 + t4_trv*100, ymin = idxx-3, ymax = idxx +3, colors = 'magenta', linestyles='-', linewidth = 0.5)
                    except :
                        pass
                except :
                    pass     
                plt.title(st[0].stats.network + '.' + st[0].stats.station + '.' + st[0].stats.channel, fontsize = 20)
                plt.xlim(0,150)
                plt.ylim(-65,5)
                plt.xticks([0,25,50,75,100,125,150],[-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], fontsize=12)
                plt.vlines(x=50, ymin = -65, ymax = 5, colors='r', linestyles='--', linewidth = 0.5)
                plt.xlabel('Time Since S Arrival (s)', fontsize = 15)
                plt.yticks([])
                plt.text(-30, idxx, "          Event  " + str(idx + 1) + "\n" + st[0].stats.starttime.strftime('%Y.%m.%d %H:%M:%S'), fontsize=10)
                plt.text(152, idxx, "M$_L$ " + str(st[0].stats.sac.mag), fontsize=10) 
                plt.savefig(f'{station}.png')

            idxx -= 5

# Write header.file for plotting traveltime curve 
def write_header_file(event_directory, hypoel_arc, filename):
    
    #========================================================================
    #=========================== waveform header  ===========================
    #========================================================================
    
    data = pd.DataFrame(columns=['EVENTNUM', 'EVENTID', 'NETWORK', 'STATION', 'CHANNEL', 'EVLA', 'EVLO', 'EVDP', 'MAG', 'DIST', 'AZ', 'BAZ', 'P', 'S', 'T3', 'T4', 'T5', 'T6', 'KA', 'KT0'])
    event_list = os.listdir(event_directory)
    eventnum = 1
    for event in event_list:
        waveform_list = os.listdir(os.path.join(event_directory, event))
        for waveform in waveform_list:
            st = read(os.path.join(event_directory, event, waveform))
            sac = st[0].stats.sac
            eventid = UTCDateTime(year=sac.nzyear, julday=sac.nzjday, hour=sac.nzhour, minute=sac.nzmin, second=sac.nzsec).strftime('%Y%m%d%H%M%S')
            newdata = pd.DataFrame({
                'EVENTNUM': [eventnum], 'EVENTID': [eventid], 'NETWORK': [st[0].stats.network], 'STATION': [st[0].stats.station], 'CHANNEL': [st[0].stats.channel],
                'EVLA': [sac.evla], 'EVLO': [sac.evlo], 'EVDP': [sac.evdp], 'MAG': [sac.mag], 'DIST': [sac.dist], 'AZ': [sac.az], 'BAZ': [sac.baz],
                'P': [getattr(sac, 'a', -12345)], 'S': [getattr(sac, 't0', -12345)], 'T3': [getattr(sac, 't3', -12345)], 'T4': [getattr(sac, 't4', -12345)],
                'T5': [getattr(sac, 't5', -12345)], 'T6': [getattr(sac, 't6', -12345)], 'KA': [getattr(sac, 'ka', -12345)], 'KT0': [getattr(sac, 'kt0', -12345)]})
            data = pd.concat([data, newdata], ignore_index=True)
        eventnum += 1
    data = data[(data.P != -12345) | (data.S != -12345) | (data.T3 != -12345) | (data.T4 != -12345)]
    data.reset_index(drop=True, inplace=True)
    
    #========================================================================
    #============================== hypoel.arc ==============================
    #========================================================================
    
    hypoel = pd.read_csv(hypoel_arc, header = None)
    hypoel['ORIGIN TIME'] = np.NaN
    hypoel['DEPTH'] = np.NaN

    # Extract origin time based on the format of hypoel.arc
    filtered_rows = hypoel[hypoel[0].astype(str).str.startswith('2')]
    filtered_rows = filtered_rows.apply(lambda row: ' '.join(row.astype(str)), axis=1)
    filtered_rows = filtered_rows.str[:16].str.replace(' ', '0')
    filtered_rows = filtered_rows.str[:8] + 'T' + filtered_rows.str[8:]
    filtered_rows = filtered_rows.str[:15] + '.' + filtered_rows.str[15:]
    
    filtered_rows2 = hypoel[hypoel[0].astype(str).str.startswith('2')]
    filtered_rows2 = filtered_rows2.apply(lambda row: ' '.join(row.astype(str)), axis=1)
    filtered_rows2 = filtered_rows2.str[31:36].astype(float) / 100

    for i in range(len(filtered_rows.index)):
        if i == len(filtered_rows.index) - 1:
            hypoel.loc[filtered_rows.index[i]:, 'ORIGIN TIME'] = filtered_rows.iloc[i]
            hypoel.loc[filtered_rows2.index[i]:, 'DEPTH'] = filtered_rows2.iloc[i]

            break
        hypoel.loc[filtered_rows.index[i]:filtered_rows.index[i+1], 'ORIGIN TIME'] = filtered_rows.iloc[i]
        hypoel.loc[filtered_rows2.index[i]:filtered_rows2.index[i+1], 'DEPTH'] = filtered_rows2.iloc[i]
        
    hypoel.drop(filtered_rows.index, inplace=True)
    hypoel.reset_index(inplace = True, drop = True)
    hypoel['STATION'] = hypoel[0].str[0:4]
    hypoel['S wave'] = hypoel[0].str[31:36]
    hypoel['P wave'] = hypoel[0].apply(lambda x: '20' + x[9:15] + 'T' + x[15:24])
    hypoel['Distance'] = hypoel[0].str[24:28]
    hypoel['Azimuth'] = hypoel[0].str[28:31]
    unique_values = hypoel['ORIGIN TIME'].unique()
    value_to_num = {value: num for num, value in enumerate(unique_values, start=1)}
    hypoel['EVENTNUM'] = hypoel['ORIGIN TIME'].map(value_to_num)
    
    for i in range(len(hypoel)):               
        if hypoel['S wave'][i].strip() == '': 
            hypoel['P wave'][i] = UTCDateTime(hypoel['P wave'][i]) - UTCDateTime(hypoel['ORIGIN TIME'][i])
        else:
            hypoel['S wave'][i] = float(hypoel['S wave'][i]) - float(hypoel['P wave'][i][-5:])
            hypoel['P wave'][i] = UTCDateTime(hypoel['P wave'][i]) - UTCDateTime(hypoel['ORIGIN TIME'][i])
            hypoel['S wave'][i] = hypoel['P wave'][i] + hypoel['S wave'][i]        

    hypoel['P wave'] = hypoel['P wave'].astype(float)
    hypoel.loc[hypoel['S wave'].str.strip() == '', 'S wave'] = np.NaN
    hypoel['S wave'] = hypoel['S wave'].astype(float)
    
    data = data.merge(hypoel, on = ['EVENTNUM', 'STATION'])
    data['HYPODEPTH'] = np.sqrt(data['DIST']**2 + data['DEPTH']**2)
    data.to_csv(filename)
    
    return data, hypoel

#plot_phase(instrument = 'korea_instrument_2024_Ahn.csv', event_directory = '../cc_final/', channel = 'HHZ', method = 'both')
#d, f = write_header_file(event_directory = './cc_final', hypoel_arc = './hypoel.arc', filename = 'ok.csv')
