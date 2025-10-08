import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from Lab 2 - Part 1a.txt
def read_lab_data(filename):
    """Read lab data from text file, skipping the header information"""
    # Read the file, skipping the first 3 lines of metadata, line 4 has headers
    data = pd.read_csv(filename, sep='\t', skiprows=3, header=0)
    # Remove any empty rows
    data = data.dropna()
    return data

# Load the data
data = read_lab_data('Lab 2 - Part 1a.txt')

# Display basic information about the data
print("Data shape:", data.shape)
print("\nColumn names:")
print(data.columns.tolist())
print("\nFirst few rows:")
print(data.head())

# Extract the columns for plotting

label1 = 'P1(PSI)'
label2 = 'T1(Deg C)'

time = data['Time(s)']
X1 = data[label1]
X2 = data[label2]

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(time, X1, 'b-', label=label1, linewidth=2)
plt.plot(time, X2, 'r-', label=label2, linewidth=2)

# Add labels and title
plt.xlabel('Time (s)')
plt.ylabel(label1)
plt.title(label1 + 'and' + label2 + ' vs Time')
plt.legend()
plt.grid(True, alpha=0.3)

# Improve the layout
plt.tight_layout()

# Display the plot
plt.show()

# Optional: Save the plot
# plt.savefig('pressure_vs_time.png', dpi=300, bbox_inches='tight')

print(f"\nData summary:")
print(f"Time range: {time.min():.2f} to {time.max():.2f} seconds")
print(f"P1 range: {P1.min():.2f} to {P1.max():.2f} PSI")
print(f"P2 range: {P2.min():.2f} to {P2.max():.2f} PSI")