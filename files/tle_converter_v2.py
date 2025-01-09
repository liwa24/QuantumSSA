import math
import pandas as pd
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import random
import pickle

# Define the reference start date
REFERENCE_DATE = datetime(1995, 1, 1)
PLOT_CUTTED = False

# Convert a date to seconds since the reference date
def to_seconds_since_1995(year, day_of_year, hour=0, minute=0, second=0):
    date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1, hours=hour, minutes=minute, seconds=second)
    return (date - REFERENCE_DATE).total_seconds()

# Parse the anomalies file
def parse_anomalies(file_path):
    anomalies = []
    with open(file_path, "r") as f:
        for line in f.readlines():
            parts = line.split()
            year = int(parts[1])
            start_seconds = to_seconds_since_1995(year, int(parts[2]), int(parts[3]), int(parts[4]))
            end_seconds = to_seconds_since_1995(int(parts[5]), int(parts[6]), int(parts[7]), int(parts[8]))
            anomalies.append({"Start Seconds": start_seconds, "End Seconds": end_seconds})
    return pd.DataFrame(anomalies)

# Parse the TLE file
def parse_tle(tle_lines):
    tle_data = []
    for i in range(0, len(tle_lines), 2):
        line1 = tle_lines[i].strip()
        line2 = tle_lines[i + 1].strip()

        year = int("20" + line1[18:20]) if int(line1[18:20]) < 25 else int("19" + line1[18:20]) #to correctly assign the year
        day_of_year = float(line1[20:32])
        inclination = float(line2[8:16]) * math.pi / 180
        omega = float(line2[17:25]) * math.pi / 180
        eccentricity = float("0." + line2[26:33])
        argument_of_perigee = float(line2[34:42]) * math.pi / 180
        mean_anomaly = float(line2[43:51]) * math.pi / 180
        mean_motion = float(line2[52:63])
        orbital_period = 86400 / mean_motion

        seconds_since_1995 = to_seconds_since_1995(year, day_of_year)
        tle_data.append({
            "Seconds Since 1995": seconds_since_1995,
            "Inclination (rad)": inclination,
            "Omega (rad)": omega,
            "Eccentricity": eccentricity,
            "Argument of Perigee (rad)": argument_of_perigee,
            "Mean Anomaly (rad)": mean_anomaly,
            "Mean Motion (rev/day)": mean_motion,
            "Orbital Period (s)": orbital_period,
        })
    return tle_data

# Main processing
def main(code=1, to_cut=False):
    tle_file_path = f"files/tle_dataset/tle_{code}.txt"
    anomalies_file_path = f"files/man_dataset/man_{code}.txt"

    # Parse TLE and anomaly data
    with open(tle_file_path, "r") as f:
        tle_lines = f.readlines()
    tle_data = parse_tle(tle_lines)
    tle_dataframe = pd.DataFrame(tle_data)

    anomalies_df = parse_anomalies(anomalies_file_path)

    # Compute derived columns
    tle_dataframe["Semi-major axis (km)"] = (
        (((tle_dataframe["Orbital Period (s)"] / (2 * math.pi)) ** 2) * 398600) ** (1 / 3)
    )
    sincos_columns = ["Inclination (rad)", "Omega (rad)", "Argument of Perigee (rad)", "Mean Anomaly (rad)"]
    for column in sincos_columns:
        tle_dataframe[f"{column} (sin)"] = tle_dataframe[column].apply(math.sin)
        tle_dataframe[f"{column} (cos)"] = tle_dataframe[column].apply(math.cos)

    # Mark maneuvers
    tle_dataframe["Maneuver"] = 0
    for _, anomaly in anomalies_df.iterrows():
        start_seconds = anomaly["Start Seconds"]
        end_seconds = anomaly["End Seconds"]
        tle_dataframe.loc[
            (tle_dataframe["Seconds Since 1995"] >= start_seconds) & 
            (tle_dataframe["Seconds Since 1995"] <= end_seconds), 
            "Maneuver"
        ] = 1

        for i in range(1, len(tle_dataframe)):
                if (
                    tle_dataframe["Seconds Since 1995"].iloc[i - 1] < start_seconds and
                    tle_dataframe["Seconds Since 1995"].iloc[i] > end_seconds
                ):
                    tle_dataframe.at[i, "Maneuver"] = 2

    # Save TLE data
    save_path = "files/tle_delta_dataframe"
    os.makedirs(save_path, exist_ok=True)
    tle_dataframe = tle_dataframe.sort_values("Seconds Since 1995")
    tle_dataframe.to_csv(os.path.join(save_path, f"tle_{code}_delta.csv"), index=False)

     # Plotting
    plot_columns = ["Semi-major axis (km)", "Eccentricity", "Inclination (rad) (sin)", "Omega (rad) (sin)", "Argument of Perigee (rad) (sin)"]
    os.makedirs("files/tle_plots", exist_ok=True)

    PLOT_CUTTED = to_cut  # Set to True to enable cutting
    cut_percentage = 0.05  # Fraction of the total time range to include (e.g., 0.5 for 50%)

    if PLOT_CUTTED:
        # Calculate the total time range
        min_time = tle_dataframe["Seconds Since 1995"].min()
        max_time = tle_dataframe["Seconds Since 1995"].max()
        total_range = max_time - min_time

        # Determine the cut range
        cut_range = total_range * cut_percentage
        start_cut = random.uniform(min_time, max_time - cut_range)
        end_cut = start_cut + cut_range

        # Filter the DataFrame for the selected range
        cut_tle_dataframe = tle_dataframe[
            (tle_dataframe["Seconds Since 1995"] >= start_cut) &
            (tle_dataframe["Seconds Since 1995"] <= end_cut)
        ]
    else:
        # Use the full DataFrame if PLOT_CUTTED is False
        cut_tle_dataframe = tle_dataframe 
        
     # Plotting
    fig = plt.figure(figsize=(12, 8))
    for i, column in enumerate(plot_columns):
        plt.subplot(3, 2, i + 1)
        cut_tle_dataframe["Formatted Date"] = cut_tle_dataframe["Seconds Since 1995"].apply(
            lambda x: (REFERENCE_DATE + timedelta(seconds=x))
            )
        
        plt.plot(cut_tle_dataframe["Formatted Date"], cut_tle_dataframe[column], label=column, alpha=0.7)

        # Highlight maneuver points
        maneuver_points = cut_tle_dataframe[cut_tle_dataframe["Maneuver"] == 1]
        plt.scatter(maneuver_points["Formatted Date"], maneuver_points[column], color="red", label="Maneuver", zorder=5)

        maneuver_points = cut_tle_dataframe[cut_tle_dataframe["Maneuver"] == 2]
        plt.scatter(maneuver_points["Formatted Date"], maneuver_points[column], color="orange", label="Maneuver", zorder=5)
        plt.xlabel("Date")
        plt.ylabel(column)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%y'))
        plt.xticks(rotation=45)
        plt.legend()

    plt.tight_layout()
    
    plt.savefig(f"files/tle_plots/tle_{code}_with_maneuvers{'_cut' if PLOT_CUTTED else ''}.pdf")
    with open(f"files/tle_plots/tle_{code}_with_maneuvers{'_cut' if PLOT_CUTTED else ''}.pkl", 'wb') as f:
        pickle.dump(fig, f)
    plt.show() 
    
output_file_summary = f"files/tle_plots/summary_{'_cut' if PLOT_CUTTED else ''}.txt"

# Check if the file exists
if not os.path.exists(output_file_summary):
    # If the file doesn't exist, create it and write a header (optional)
    with open(output_file_summary, 'w') as f:
        f.write("Maneuver Summary of Analysis:\n\n")  # Optional header

for code in range(5, 8):
    main(code=code, to_cut=PLOT_CUTTED)
    # count the number of maneuvers where the maneuver is marked as 1 or 2
    tle_dataframe = pd.read_csv(f"files/tle_delta_dataframe/tle_{code}_delta.csv")
    print(f"Number of maneuvers for TLE {code}: {tle_dataframe[tle_dataframe['Maneuver'] > 0].shape[0]}")
    # count the number of maneuvers in the "man_dataset" file
    anomalies_df = parse_anomalies(f"files/man_dataset/man_{code}.txt")
    print(f"Number of maneuvers in the anomaly file for TLE {code}: {anomalies_df.shape[0]}") 
    
    with open(output_file_summary, 'a') as f:
        f.write(f"Number of maneuvers for TLE {code}: {tle_dataframe[tle_dataframe['Maneuver'] > 0].shape[0]}\n")
        f.write(f"Number of maneuvers in the anomaly file for TLE {code}: {anomalies_df.shape[0]}\n")
        f.write("\n")  # Add a newline for better readability

#reading a file pkl
print("reading a file")
with open(f"files/tle_plots/tle_{1}_with_maneuvers{'_cut' if False else ''}.pkl", 'rb') as f:
    fig = pickle.load(f)
    plt.show()
print("done")