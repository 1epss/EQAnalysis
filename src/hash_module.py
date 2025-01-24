import numpy as np
import pandas as pd
from obspy import *
from obspy.imaging.beachball import beachball
import folium
from folium.plugins import BeautifyIcon
import pygmt
import os
import re

# Function for hash_write_phase
def degree_to_dms(deg):
    degree = int(deg)
    minute = int((deg - degree) * 60)
    second = (deg - degree - minute/60) * 3600
    second = int(round(second, 0))
    return degree, minute, second

# Function for hash_draw_beachball
def normalize(series, new_min=0.25, new_max=0.5):
    min_val = series.min()
    max_val = series.max()
    normalized = (series - min_val) / (max_val - min_val) * (new_max - new_min) + new_min
    return normalized


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


# Draw and save beachball diagrams for each seismic event with obspy.imaging.beachball.beachball
def hash_draw_beachball_obspy(hash_output, output_dir):
    columns = [
        "Event_number", "Year", "Month", "Day", "Event_type", "Magnitude", "Magnitude_type", "Latitude", "Longitude", "Depth", "Location_quality", "RMS", "Horizontal_error", "Depth_error", "Origin_time_error", "Number_of_picks", "Number_of_P_picks", "Number_of_S_picks", "Strike", "Dip", "Rake", "Fault_plane_uncertainty", "Auxiliary_plane_uncertainty", "Number_of_P_first_motion_polarities", "Weighted_percent_misfit_of_first_motions", "Quality", "Probability_mechanism", "Station_distribution_ratio", "Number_of_S/P_ratios", "Average_log10(S/P)_ratio", "Multiple_flag"
    ]
    csv = pd.read_csv(hash_output, sep='\s+', names=columns)
    os.makedirs(output_dir, exist_ok=True)

    for _, row in csv.iterrows():
        event_num, strike, dip, rake = int(row['Event_number']), int(row['Strike']), int(row['Dip']), int(row['Rake'])
        existing_files = [f for f in os.listdir(output_dir) if f.startswith(f"{event_num}") and f.endswith(".png")]
        file_name = f"{event_num}-{len(existing_files)}.png" if existing_files else f"{event_num}.png"
        file_path = os.path.join(output_dir, file_name)
        beachball([strike, dip, rake], size=200, linewidth=2, facecolor='b', outfile=file_path)


# Draw and save beachball diagrams for each seismic event with pygmt.Figure.meca
def hash_draw_beachball_pygmt(hypoel_arc, amp, hash_output, output_dir, event_num = 100000):
    # *============= hypoel.arc ================*
    columns = ['Station', 'P_polarity', 'P_arrivaltime', 'Distance', 'Azimuth', 'S_arrivaltime', 'S_polarity', 'Takeoff', 'Event_number']
    data = []
    with open(hypoel_arc, 'r') as f:
        for arc_line in f:
            if re.match(r'^[A-Za-z]', arc_line):
                data.append({
                    'Station': arc_line[0:4].strip(),
                    'P_polarity': arc_line[4:8].strip(),
                    'P_arrivaltime': arc_line[9:24].strip(),
                    'Distance': arc_line[24:28].strip(),
                    'Azimuth': arc_line[28:31].strip(),
                    'S_arrivaltime': arc_line[31:36].strip(),
                    'S_polarity': arc_line[36:40].strip(),
                    'Takeoff': arc_line[40:43].strip(),
                    'Event_number' : str(event_num) 
                })
            elif re.match(r'^[0-9]', arc_line) :
                data.append({
                    'Station' : np.nan,
                    'P_polarity': np.nan,
                    'P_arrivaltime': np.nan,
                    'Distance': np.nan,
                    'Azimuth': np.nan,
                    'S_arrivaltime': np.nan,
                    'S_polarity': np.nan,
                    'Takeoff': np.nan,
                    'Event_number': np.nan
                })
                event_num += 1 
            else :
                continue
    arc_df = pd.DataFrame(data, columns=columns)
    arc_sections = []
    arc_current_section = None
    for _, row in arc_df.iterrows():
        if pd.isna(row['Station']):
            if arc_sections is not None:
                arc_sections.append(arc_current_section)
            arc_current_section = pd.DataFrame(columns = columns)
        if arc_current_section is not None:
            arc_current_section = pd.concat([arc_current_section, pd.DataFrame([row], columns=columns)], ignore_index=True)
        arc_current_section = arc_current_section.dropna(subset=['Station'])
    if arc_current_section is not None:
        arc_sections.append(arc_current_section)
    arc_sections = arc_sections[1:]
    # *============= amp.hash ================*
    columns = ['Station', 'Channel', 'Network', 'P_noise_level', 'S_noise_level', 'P_level', 'S_level', 'Event_number']
    data = []
    with open(amp, 'r') as f:
        for amp_line in f:
            if re.match(r'^[A-Za-z]', amp_line):
                data.append({
                    'Station': amp_line[0:4].strip(),
                    'Channel': amp_line[6:8].strip(),
                    'Network': amp_line[10:11].strip(),
                    'P_noise_level': amp_line[29:38].strip(),
                    'S_noise_level': amp_line[40:49].strip(),
                    'P_level': amp_line[51:60].strip(),
                    'S_level': amp_line[62:71].strip(),
                    'Event_number': str(event_num)
                    })
            elif re.match(r'^[0-9]', amp_line) :
                data.append({
                    'Station': np.nan,
                    'Channel': np.nan,
                    'Network': np.nan,
                    'P_noise_level': np.nan,
                    'S_noise_level': np.nan,
                    'P_level': np.nan,
                    'S_level': np.nan,
                    'Event_number': np.nan
                    })
                event_num = amp_line.split()[0]
            else :
                continue
    amp_df = pd.DataFrame(data, columns=columns)
    amp_sections = []
    amp_current_section = None
    for _, row in amp_df.iterrows():
        if pd.isna(row['Station']):
            if amp_current_section is not None:
                amp_sections.append(amp_current_section)
            amp_current_section = pd.DataFrame(columns = columns)
        if amp_current_section is not None:
            amp_current_section = pd.concat([amp_current_section, pd.DataFrame([row], columns=columns)], ignore_index=True)
        amp_current_section = amp_current_section.dropna(subset=['Station'])
    if amp_current_section is not None:
        amp_sections.append(amp_current_section)
    # normalize
    for i in range(len(amp_sections)):
        amp_sections[i]['P_level'] = amp_sections[i]['P_level'].astype('float')
        amp_sections[i]['Normalized'] = normalize(amp_sections[i]['P_level'])
    # *============= plotting ================*
    columns = [
        "Event_number", "Year", "Month", "Day", "Hour", "Minute", "Second", "Event_type", "Magnitude", "Magnitude_type", "Latitude", "Longitude", "Depth", "Location_quality", "RMS", "Horizontal_error", "Depth_error", "Origin_time_error", "Number_of_picks", "Number_of_P_picks", "Number_of_S_picks", "Strike", "Dip", "Rake", "Fault_plane_uncertainty", "Auxiliary_plane_uncertainty", "Number_of_P_first_motion_polarities", "Weighted_percent_misfit_of_first_motions", "Quality", "Probability_mechanism", "Station_distribution_ratio", "Number_of_S/P_ratios", "Average_log10(S/P)_ratio", "Multiple_flag"
    ]
    csv = pd.read_csv(hash_output, sep = '\s+', names = columns)
    os.makedirs(output_dir, exist_ok=True)

    for idx, row in csv.iterrows():
        event_num, strike, dip, rake = int(row['Event_number']), int(row['Strike']), int(row['Dip']), int(row['Rake'])
        fig = pygmt.Figure()
        pygmt.config(FONT_TITLE="14p,Helvetica,black", FORMAT_GEO_MAP="+D")
        fig.basemap(
            region=[0, 360, 0, 1],
            projection="P5c+a",
            frame=["xa45f", "+gwhite"])
        # 주향을 270도 회전 (360도를 넘지 않도록 조정)
        if strike >= 90:
            rotated_strike = (strike + 270) % 360
        else :
            rotated_strike = (strike + 270)
        fig.meca(spec=np.array([0, 0, float(row['Depth']), rotated_strike, dip, rake, float(row['Magnitude'])]), scale="5c+m", convention="aki")
        if idx == 0:
            current_event_num = np.nan        
        if current_event_num == event_num : 
            fig.text(position="TC", text=f"Event {event_num}(2)", offset="0/1.5c", no_clip=True)
            multiple = True
        else :
            fig.text(position="TC", text=f"Event {event_num}", offset="0/1.5c", no_clip=True)
            multiple = False
        current_event_num = event_num
        
        # Station take-off angle and azimuth for each event
        arc_combined = pd.concat(arc_sections, ignore_index=True)
        amp_combined = pd.concat(amp_sections, ignore_index=True)
        arc_filtered = arc_combined.loc[arc_combined['Event_number'] == str(event_num)]
        amp_filtered = amp_combined.loc[amp_combined['Event_number'] == str(event_num)]

        for _, arc in arc_filtered.iterrows():
            x = int(arc['Azimuth'])
            y = float(arc['Takeoff']) / 100
            sta = arc['Station']
            if y >= 0.90:
                x += 180 
                y = abs(0.90 - (y - 0.90))
            try:
                if 'U' in arc['P_polarity']:
                    fig.plot(x = x, y = y, size = [amp_filtered.loc[amp_filtered['Station'] == sta]['Normalized'].values[0]], style='cc', fill='black', pen='white')  # x : azimuth, y : takeoff
                    fig.text(x = x, y = y, text= sta, font = "2.5p,white", transparency = 0)
                elif 'D' in arc['P_polarity']:
                    fig.plot(x = x, y = y, size = [amp_filtered.loc[amp_filtered['Station'] == sta]['Normalized'].values[0]], style='cc', fill='white', pen='black')  # x : azimuth, y : takeoff
                    fig.text(x = x, y = y, text= sta, font = "2.5p", transparency = 0)
            except:
                continue
        
        if multiple == True:
            fig.savefig(f'{output_dir}/{event_num}(2).png')
        else :
            fig.savefig(f'{output_dir}/{event_num}.png')


# Project beachball diagrams onto a map using Folium
def hash_map_beachball_folium(hash_output, output_dir):
    '''
    df_KS_station = pd.read_csv('station.txt', sep = ' ', header = None)
    df_KS_station.columns = ['NETWORK', 'STATION', 'LATITUDE', 'LONGITUDE', 'ELEVATION']




    # 지도 생성 (중심 좌표 설정)
    center = [35.8034 ,127.5311]

    #m = folium.Map(width=900, height=900, location=center, zoom_start=10, control_scale = True, tiles="cartodbpositron")
    #m = folium.Map(width=900, height=900, location=center, zoom_start=10, control_scale = True, tiles="cartodbdark_matter") 

    m = folium.Map(width=900, height=900, location=center, zoom_start=16, control_scale = True, 
                   tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                   attr='Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community',
                   name='Esri World Imagery')

    event = pd.read_csv('event.csv')
    image_list = os.listdir('./image')

    for i in range(1, 14):
        image = './image/' + str(i) + '.png'

        if not os.path.isfile(image):
            print(f"Could not find {image}")

        else:
            img = folium.raster_layers.ImageOverlay(
                name="Mercator projection SW",
                image=image,
                bounds=[[event.iloc[i - 1][1]-0.0005, event.iloc[i - 1][2]-0.0005], 
                        [event.iloc[i - 1][1]+0.0005, event.iloc[i - 1][2]+0.0005]],
                opacity=0.6,
                interactive=True,
                cross_origin=False,
                zindex=1).add_to(m)

        icon_circle = BeautifyIcon(icon='circle', inner_icon_style='color:red;font-size:5px;', background_color='transparent', border_color='transparent')
        folium.map.Marker((event.iloc[i - 1][1],event.iloc[i - 1][2]), tooltip = f'{i}', icon=folium.features.DivIcon(icon_size=(0,0),icon_anchor=(0,0),html=f'<div style="font-size: 8pt; color: {"white"}">{i}</div>')).add_to(m)  


    for idx, row in df_KS_station.iterrows() :
        folium.features.RegularPolygonMarker(location=(row.LATITUDE, row.LONGITUDE), tooltip = f'station:{row.STATION}<br/>Network:{row.NETWORK}<br/>Location:{row.LATITUDE:.4f},{row.LONGITUDE:.4f}', 
                                             color='yellow', fill_color='green',number_of_sides=6, rotation=30, radius=5, fill_opacity=1).add_to(m)
        folium.map.Marker((row.LATITUDE,row.LONGITUDE), icon=folium.features.DivIcon(icon_size=(0,0),icon_anchor=(0,-20),html=f'<div style="font-size: 8pt; color: {"white"}">{row.STATION}</div>')).add_to(m)



    folium.Circle(location=center, color='black', fill_opacity=0, radius=50*1e+3, label = '50km').add_to(m)
    folium.Circle(location=center, color='black', fill_opacity=0, radius=100*1e+3, label = '100km').add_to(m)

    plugins.Fullscreen(
    position='topright',
    title='Expand me',
    title_cancel='Exit me',
    force_separate_button=True ).add_to(m)    

    m
    # 지도를 HTML 파일로 저장
    m.save('jangsu_focalmechanism.html')
    '''