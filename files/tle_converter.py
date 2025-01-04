import math
import pandas as pd
import os

from datetime import datetime, timedelta

def compute_decimal_day(day, hour, minute, second):
    return day + hour / 24 + minute / 1440 + second / 86400

# Parse the anomaly file
def parse_anomalies(file_path):
    anomalies = []
    with open(file_path, "r") as f:
        for line in f.readlines():
            parts = line.split()
            year = int(parts[1])
            start_day = int(parts[2])
            start_hour = int(parts[3])
            start_min = int(parts[4])
            start_sec = float(parts[5])
            end_day = int(parts[6])
            end_hour = int(parts[7])
            end_min = int(parts[8])
            end_sec = float(parts[9])
            
            start_decimal_day = compute_decimal_day(start_day, start_hour, start_min, start_sec)
            end_decimal_day = compute_decimal_day(end_day, end_hour, end_min, end_sec)
            
            anomalies.append({
                "Year": year,
                "Start Decimal Day": start_decimal_day,
                "End Decimal Day": end_decimal_day,
            })
    return pd.DataFrame(anomalies)

# Parse the anomalies
#anomalies_file_path = "files/man_1.txt"
#anomalies_df = parse_anomalies(anomalies_file_path)

def compute_julian_date(year, day_of_year):
    # Ensure year is an int
    year = int(year)

    date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
    starting_date = datetime(1995, 1, 1)
    date -= starting_date
    date = date.total_seconds()
    return date / (24 * 3600)



def parse_tle(tle_lines):
    tle_data = []
    
    for i in range(0, len(tle_lines), 2):
        line1 = tle_lines[i].strip()
        line2 = tle_lines[i+1].strip()
        
        name = line1[2:7]  # Name of the satellite
        # Extracting data from Line 1
        year = int('20' + line1[18:20])  # Epoch year
        day_of_year = float(line1[20:32])  # Epoch day of the year and fractional part
        
        # Extracting data from Line 2
        inclination = float(line2[8:16]) * math.pi / 180  # Inclination in radians
        omega = float(line2[17:25]) * math.pi / 180  # Right Ascension of Ascending Node in radians
        eccentricity = float('0.' + line2[26:33])  # Eccentricity
        argument_of_perigee = float(line2[34:42]) * math.pi / 180  # Argument of perigee in radians
        mean_anomaly = float(line2[43:51]) * math.pi / 180  # Mean anomaly in radians
        mean_motion = float(line2[52:63])  # Mean motion (revolutions per day)
        orbital_period = 86400 / mean_motion  # Orbital period in seconds

        tle_data.append({
            "Catalog Number": name,
            "Year": year,
            "Day of Year": day_of_year,
            "Inclination (rad)": inclination,
            "Omega (rad)": omega,
            "Eccentricity": eccentricity,
            "Argument of Perigee (rad)": argument_of_perigee,
            "Mean Anomaly (rad)": mean_anomaly,
            "Mean Motion (rev/day)": mean_motion,
            "Orbital Period (s)": orbital_period,
        })
    
    return tle_data



# read path from . tel_dataset\tel_1.txt
tle_lines = []
with open("files/tle_dataset/tle_1.txt", "r") as f:
    tle_lines = f.readlines()
    
parsed_tle_data = parse_tle(tle_lines)

for i, tle in enumerate(parsed_tle_data):
    print(f"TLE {i+1}:")
    for key, value in tle.items():
        print(f"  {key}: {value}")
    print()
    if i > 0:
        break
tle_dataframe = pd.DataFrame(parsed_tle_data)

tle_dataframe["Semi-major axis (km)"] = (((tle_dataframe["Orbital Period (s)"]/2/math.pi) ** 2) * 398600)** (1/3) #costante sbagliata
# Add Julian Date column
tle_dataframe["Julian Date"] = tle_dataframe.apply(
    lambda row: compute_julian_date(row["Year"], row["Day of Year"]), axis=1
)
sincos_columns = ["Inclination (rad)", "Omega (rad)", "Argument of Perigee (rad)", "Mean Anomaly (rad)"]
for column in sincos_columns:
    tle_dataframe[f"{column} (sin)"] = tle_dataframe[column].apply(math.sin)
    tle_dataframe[f"{column} (cos)"] = tle_dataframe[column].apply(math.cos)

save_path = "files/tle_dataframe"
# create folder save_path

if not os.path.exists(save_path):
    os.makedirs(save_path)
tle_dataframe.to_csv(save_path + "/tle_1.csv", index=False)

print("TLE data saved to CSV file")
print(tle_dataframe.head())

# Calculate delta values for relevant columns
relevant_columns = [
    "Julian Date",
    "Inclination (rad)",
    "Omega (rad)",
    "Eccentricity",
    "Argument of Perigee (rad)",
    "Mean Anomaly (rad)",
    "Mean Motion (rev/day)",
    "Orbital Period (s)",
    "Semi-major axis (km)"
]

#create a deepcopy of the dataframe
tle_delta_dataframe = tle_dataframe.copy()

# Add delta columns
for column in relevant_columns:
    if column not in sincos_columns:
        tle_delta_dataframe[f"Delta {column}"] = tle_delta_dataframe[column].diff().fillna(0)
    else:
        # attenzione: da sistemare per le differenze di angoli
        value = tle_delta_dataframe[column].diff().fillna(0)
        tle_delta_dataframe[f"Delta {column}"] = value

tle_delta_dataframe["Maneuver"] = 0

anomalies_file_path = "files/man_dataset/man_1.txt"
anomalies_df = parse_anomalies(anomalies_file_path)
print(anomalies_df)

for _, anomaly in anomalies_df.iterrows():
    year = anomaly["Year"]
    start_decimal_day = anomaly["Start Decimal Day"]
    end_decimal_day = anomaly["End Decimal Day"]

    start_julian_date = compute_julian_date(year, start_decimal_day)
    end_julian_date = compute_julian_date(year, end_decimal_day)
    
    # Check if a maneuver occurred within the time frame
    tle_delta_dataframe.loc[
        (tle_delta_dataframe["Julian Date"] >=  start_julian_date) &
        (tle_delta_dataframe["Julian Date"] <= end_julian_date), "Maneuver"
    ] = 1

# Save the updated dataframe to a new folder
delta_save_path = "files/tle_delta_dataframe"
if not os.path.exists(delta_save_path):
    os.makedirs(delta_save_path)

tle_delta_dataframe_path = os.path.join(delta_save_path, "tle_1_delta.csv")
tle_delta_dataframe.to_csv(tle_delta_dataframe_path, index=False)

print("TLE data with delta values saved to CSV file")
print(tle_delta_dataframe.head())

tle_delta_normalized = tle_delta_dataframe.copy()
# normalize the non-delta columns
#for column in relevant_columns:
#    max_val = tle_delta_normalized[column].max()
#    min_val = tle_delta_normalized[column].min()
#    tle_delta_normalized[column] = (tle_delta_normalized[column] - min_val) / (max_val - min_val)

# to plot: a e i equal to 
plot_vals = ["Semi-major axis (km)", "Eccentricity", "Inclination (rad) (sin)", "Omega (rad) (sin)", "Argument of Perigee (rad) (sin)"] #Omega RAAN

# select the first 10 elements
# plot the values for x equal to the day of year
import matplotlib.pyplot as plt
import numpy as np

#create tlt_plots folder
if not os.path.exists("files/tle_plots"):
    os.makedirs("files/tle_plots")

plt.figure(figsize=(12, 8))
for i, column in enumerate(plot_vals):
    plt.subplot(3, 2, i+1)
    plt.plot(tle_delta_dataframe["Julian Date"], tle_delta_dataframe[column], label=column, alpha=0.7)
    # Highlight maneuver points
    maneuver_points = tle_delta_dataframe[tle_delta_dataframe["Maneuver"] == 1]
    plt.scatter(maneuver_points["Julian Date"], maneuver_points[column], color="red", label="Maneuver", zorder=5)
    plt.xlabel("Julian Date")
    plt.ylabel(column)
    plt.legend()

plt.tight_layout()
plt.savefig("files/tle_plots/tle_1_with_maneuvers.pdf")
plt.show()