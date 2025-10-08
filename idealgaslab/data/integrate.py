import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import trapezoid

# Load data file
data_path = Path("../rawdata/Lab 2 - Part 2a.txt")
with data_path.open(encoding="utf-8") as f:
    lines = f.readlines()

# Find header and read data
header_idx = next(i for i, line in enumerate(lines) if line.strip().startswith("Time(s)"))
df = pd.read_csv(data_path, sep=r"\s*\t\s*", engine="python", skiprows=header_idx)

# Extract and convert data
t = df["Time(s)"].to_numpy()
m_dot = df["Mass Flowrate(g/min)"].to_numpy() / 60000  # Convert g/min to kg/s

# Calculate total mass using different methods
riemann_mass = np.sum(m_dot[:-1] * np.diff(t))  # Left Riemann sum
trap_mass = trapezoid(m_dot, t)  # Trapezoidal rule

print(f"Total mass (Riemann): {riemann_mass:.6f} kg")
print(f"Total mass (Trapezoidal): {trap_mass:.6f} kg")

# Visualization
plt.figure(figsize=(10, 6))
plt.plot(t, m_dot, 'b-', label="Mass flow rate")
plt.fill_between(t, 0, m_dot, step="pre", alpha=0.3, label="Riemann rectangles")
plt.xlabel("Time (s)")
plt.ylabel("Mass flow rate (kg/s)")
plt.title("Mass Flow Rate Integration")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
