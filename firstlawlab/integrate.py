from process_data import load_data_file
import numpy as np
from scipy import integrate


def integrate_mass_flow(filepath, t_start=None, t_end=None):
    """
    Load data file and integrate mass flow rate with respect to time.
    
    Parameters:
    -----------
    filepath : str
        Path to the data file
    t_start : float, optional
        Start time for integration (seconds). If None, uses first time point.
    t_end : float, optional
        End time for integration (seconds). If None, uses last time point.
        
    Returns:
    --------
    total_mass : float
        Integrated mass (in grams)
    time_range : tuple
        (t_start, t_end) actual time range used for integration
    """
    # Load the data
    df, metadata = load_data_file(filepath)
    
    # Extract time and mass flow rate
    time = df['Time(s)'].values
    mass_flow_rate = df['Mass Flow rate(g/min)'].values
    
    # Convert mass flow rate from g/min to g/s for integration with time in seconds
    mass_flow_rate_per_sec = mass_flow_rate / 60.0
    
    # Determine integration limits
    if t_start is None:
        t_start = time[0]
    if t_end is None:
        t_end = time[-1]
    
    # Filter data to integration range
    mask = (time >= t_start) & (time <= t_end)
    time_filtered = time[mask]
    mass_flow_filtered = mass_flow_rate_per_sec[mask]
    
    # Perform integration using trapezoidal rule
    total_mass = integrate.trapezoid(mass_flow_filtered, time_filtered)
    
    return total_mass, (t_start, t_end)


if __name__ == "__main__":
    # Example usage
    import os
    
    # Example file path - adjust as needed
    example_file = os.path.join('rawdata', 'Lab 2 - Part 1d.txt')
    
    if os.path.exists(example_file):
        # Integrate over entire time range
        total_mass, time_range = integrate_mass_flow(example_file)
        print(f"Total mass (entire range): {total_mass:.2f} g")
        print(f"Time range: {time_range[0]:.2f} to {time_range[1]:.2f} seconds")
        
        # Example with custom time limits
        total_mass_custom, time_range_custom = integrate_mass_flow(example_file, t_start=10, t_end=48)
        print(f"\nTotal mass (10-100s): {total_mass_custom:.2f} g")
        print(f"Time range: {time_range_custom[0]:.2f} to {time_range_custom[1]:.2f} seconds")
    else:
        print(f"Example file not found: {example_file}")