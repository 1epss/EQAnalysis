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
        print("traveltime file created successfully as '{filename}'.")

    elif module_name == "write_phase_file":
        event = input("Enter the path to the event CSV file: ").strip()
        traveltime = input("Enter the path to the traveltime file: ").strip()
        filename = input("Enter the output filename: ").strip()
        write_phase_file(event, traveltime, filename)
        print(f"phase file created successfully as '{filename}'.")

    else:
        print(f"Invalid module name: {module_name}. Please select a valid module.")

if __name__ == "__main__":
    main()
