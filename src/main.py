# src/main.py

# Import all modules from *_module.py
from hash_module import *
from hypoDD_module import *

def main():
    print("For running Hash, these modules are available:")
    print("========================================")
    print("1. hash_write_station")
    print("2. hash_write_saclst")
    print("3. hash_write_phase")
    print("4. hash_write_amp")
    print("5. hash_draw_beachball_obspy")
    print("6. hash_draw_beachball_pygmt")
    print("========================================")
    print("========================================")
    print("For running HypoDD, these modules are available:")
    print("========================================")
    print("1. hypoDD_write_station")
    print("2. hypoDD_write_traveltime")
    print("3. hypoDD_write_phase")
    print("4. hypoDD_write_cc")
    print("========================================")

    # Get the module name from the user
    module_name = input("Enter the module name to execute (e.g., write_station_file): ").strip()
    
    if module_name == "hash_write_station":
        instrument = input("Enter the path to the instrument CSV file: ").strip()
        filename = input("Enter the output filename: ").strip()
        hash_write_station(instrument, filename)
        print(f"Station file created successfully as '{filename}'.")

    elif module_name == "hash_write_saclst":
        event_directory = input("Enter the path to the event_directory: ").strip()
        filename = input("Enter the output filename: ").strip()
        hash_write_saclst(event_directory, filename)
        print(f"Saclst file created successfully as '{filename}'.")

    elif module_name == "hash_write_phase":
        event_csv = input("Enter the path to the event CSV file: ").strip()
        saclst_file = input("Enter the path to the saclst file: ").strip()
        event_directory = input("Enter the path to the event_directory: ").strip()
        filename = input("Enter the output filename: ").strip()
        hash_write_phase(event_csv, saclst_file, event_directory, filename)
        print(f"Phase file created successfully as '{filename}'.")

    elif module_name == "hash_write_amp":
        event_csv = input("Enter the path to the event CSV file: ").strip()
        saclst_file = input("Enter the path to the saclst file: ").strip()
        event_directory = input("Enter the path to the event_directory: ").strip()
        filename = input("Enter the output filename: ").strip()
        method = input("Enter the method: ")
        hash_write_amp(event_csv, saclst_file, event_directory, filename, method)
        print(f"Amplitude file created successfully as '{filename}'.")

    elif module_name == "hash_draw_beachball_obspy":
        hash_output = input("Enter the path to the hash.output file: ").strip()
        output_dir = input("Enter the output directory: ").strip()
        hash_draw_beachball_obspy(hash_output, output_dir)
        print(f"Beachball files created successfully.")

    elif module_name == "hash_draw_beachball_pygmt":
        hypoel_arc = input("Enter the path to the hypoel.arc file: ").strip()
        amp = input("Enter the path to the amp.hash file: ").strip()
        hash_output = input("Enter the path to the hash.output file: ").strip()
        output_dir = input("Enter the output directory: ").strip()
        event_num = int(input("Enter the number of first event in hash format: ").strip())
        hash_draw_beachball_pygmt(hypoel_arc, amp, hash_output, output_dir, event_num)
        print(f"Beachball files created successfully.")

    #==============================================================================

    elif module_name == "hypoDD_write_station":
        instrument = input("Enter the path to the instrument CSV file: ").strip()
        filename = input("Enter the output filename: ").strip()
        hypoDD_write_station(instrument, filename)
        print(f"Station file created successfully as '{filename}'.")

    elif module_name == "hypoDD_write_traveltime":
        hypoel_out = input("Enter the path to the hypoel.out file: ").strip()
        filename = input("Enter the output filename: ").strip()
        hypoDD_write_traveltime(hypoel_out, filename)
        print(f"Traveltime file created successfully as '{filename}'.")

    elif module_name == "hypoDD_write_phase":
        event_csv = input("Enter the path to the event CSV file: ").strip()
        traveltime = input("Enter the path to the traveltime file: ").strip()
        filename = input("Enter the output filename: ").strip()
        hypoDD_write_phase(event_csv, traveltime, filename)
        print(f"Phase file created successfully as '{filename}'.")
    
    elif module_name == "hypoDD_write_cc":
        event_directory = input("Enter the path to the event_directory: ").strip()
        station = input("Enter the path to the station file: ").strip()
        phase = input("Enter the specific phase to be computed: ").strip()
        btime = float(input("Enter the btime of waveform: ").strip())
        filename = input("Enter the output filename: ").strip()
        channel_list = input("Enter the channel_list: ")
        resampling_rate = input("Enter the resampling_rate: ").strip()
        freqmin = input("Enter the freqmin: ").strip()
        freqmax = input("Enter the freqmax: ").strip()
        p_minsec = input("Enter the p_minsec: ").strip()
        p_maxsec = input("Enter the p_maxsec: ").strip()
        s_minsec = input("Enter the s_minsec: ").strip()
        s_maxsec = input("Enter the s_maxsec: ").strip()
        threshold = float(input("Enter the threshold: ").strip())
        if channel_list == "":
            channel_list = ['HHZ', 'ELZ', 'HGZ'],
        if resampling_rate == "":
            resampling_rate = 1000.0
        else:
            resampling_rate = float(resampling_rate)
        if freqmin == "":
            freqmin = 4.0
        else:
            freqmin = float(freqmin)
        if freqmax == "":
            freqmax = 20.0
        else:
            freqmax = float(freqmax)
        if p_minsec == "":
            p_minsec = 0.1
        else:
            p_minsec = float(p_minsec)
        if p_maxsec == "":
            p_maxsec = 0.3
        else:
            p_maxsec = float(p_maxsec)
        if s_minsec == "":
            s_minsec = 0.1
        else:
            s_minsec = float(s_minsec)
        if s_maxsec == "":
            s_maxsec = 0.4
        else:
            s_maxsec = float(s_maxsec)
        hypoDD_write_cc(event_directory, station, filename, phase, btime, channel_list, resampling_rate, freqmin, freqmax, p_minsec, p_maxsec, s_minsec, s_maxsec, threshold)
        print(f"Cross-correlation file created successfully as '{filename}'.")

    else:
        print(f"Invalid module name: {module_name}. Please select a valid module.")

if __name__ == "__main__":
    main()