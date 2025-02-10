
def convert_to_seconds(time_str):
    days, time = time_str.split('-')
    hours, minutes, seconds = time.split(':')  
    days = int(days)
    hours = int(hours)
    minutes = int(minutes)
    seconds = int(seconds)
    total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
    return total_seconds

import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('speedup.dat', delimiter=" ")
data["Time"]=data['Time'].apply(convert_to_seconds)
data['Speedup'] = data['Time'].iloc[0] / data['Time']

plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(data['CPU'], data['CPU'], label="ideal", linewidth=3, color="gray")
plt.plot(data['CPU'], data['Speedup'], 'o-', label="real", linewidth=3, markersize=30, color="blue")
plt.xlabel('Number of CPUs [-]')
plt.ylabel('Speed-up [-]')
plt.legend()
plt.grid()
plt.savefig('speedup.jpg', bbox_inches='tight')
