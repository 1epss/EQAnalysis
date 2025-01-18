# src/main.py

# Import all modules from *_module.py
from hypoDD_module import *

def main():
    print("========================================")
    print("For running HypoDD, these modules are available:")
    print("========================================")
    print("1. write_station_file")
    print("2. write_traveltime_file")
    print("3. write_phase_file")
    print("4. write_cc_file")
    print("========================================")

    # Get the module name from the user
    module_name = input("Enter the module name to execute (e.g., write_station_file): ").strip()
    
    if module_name == "write_station_file":
        instrument = input("Enter the path to the instrument CSV file: ").strip()
        filename = input("Enter the output filename: ").strip()
        write_station_file(instrument, filename)
        print(f"station file created successfully as '{filename}'.")

    elif module_name == "write_traveltime_file":
        hypoel_out = input("Enter the path to the hypoel.out file: ").strip()
        filename = input("Enter the output filename: ").strip()
        write_traveltime_file(hypoel_out, filename)
        print(f"traveltime file created successfully as '{filename}'.")

    elif module_name == "write_phase_file":
        event = input("Enter the path to the event CSV file: ").strip()
        traveltime = input("Enter the path to the traveltime file: ").strip()
        filename = input("Enter the output filename: ").strip()
        write_phase_file(event, traveltime, filename)
        print(f"phase file created successfully as '{filename}'.")
    
    elif module_name == "write_cc_file":
        event_directory = input("Enter the path to the event_directory: ").strip()
        station = input("Enter the path to the station file: ").strip()
        phase = input("Enter the specific phase to be computed: ").strip()
        btime = float(input("Enter the btime of waveform: ").strip())
        filename = input("Enter the output filename: ").strip()
        channel_list = input("Enter the channel_list: ")
        resampling_rate = input("Enter the resampling_rate: ").strip()
        freqmin = input("Enter the freqmin: ").strip()
        freqmax = input("Enter the freqmax: ").strip()
        minsec = input("Enter the minsec: ").strip()
        maxsec = input("Enter the maxsec: ").strip()
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
        if minsec == "":
            minsec = 0.1
        else:
            minsec = float(minsec)
        if maxsec == "":
            maxsec = 0.4
        else:
            maxsec = float(maxsec)
        write_cc_file(event_directory, station, filename, phase, btime, channel_list, resampling_rate, freqmin, freqmax, minsec, maxsec, threshold)
        print(f"cc file created successfully as '{filename}'.")

    else:
        print(f"Invalid module name: {module_name}. Please select a valid module.")

if __name__ == "__main__":
    main()