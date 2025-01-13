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