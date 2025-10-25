import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path


def load_data_file(filepath):
    """
    Load a data file from firstlawlab/rawdata directory.
    
    Handles two formats:
    1. Files with headers (ambient conditions + column names)
    2. Files without headers (only numeric data)
    
    Parameters:
    -----------
    filepath : str
        Path to the data file
        
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame with properly formatted data
    metadata : dict
        Dictionary containing ambient conditions (if present)
    """
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
    
    metadata = {}
    
    # Check if file has headers (starts with "Ambient")
    if first_line.startswith("Ambient"):
        # Read metadata from first two lines
        with open(filepath, 'r') as f:
            line1 = f.readline().strip()
            line2 = f.readline().strip()
            
            # Parse ambient temperature
            if "Ambient Temperature" in line1:
                metadata['ambient_temp_C'] = float(line1.split(':')[1].strip())
            
            # Parse ambient pressure
            if "Ambient Pressure" in line2:
                metadata['ambient_pressure_inHg'] = float(line2.split(':')[1].strip())
        
        # Define standard column names (header in file is split across 2 lines and hard to parse)
        column_names = [
            'Time(s)',
            'T1(Deg C)',
            'Wall Temp(Deg C)',
            'Top Plate(Deg C)',
            'Bottom Plate(Deg C)',
            'P1(PSI)',
            'P2(PSI)',
            'Mass Flow rate(g/min)',
            'Heater Energy (kJ)'
        ]
        
        # Read data starting from line 4 (skip first 4 lines: 2 metadata + 2 header lines)
        # Use whitespace separator to handle irregular spacing
        df = pd.read_csv(filepath, sep=r'\s+', skiprows=4, names=column_names, engine='python')
        
    else:
        # File has no headers, define standard column names
        column_names = [
            'Time(s)',
            'T1(Deg C)',
            'Wall Temp(Deg C)',
            'Top Plate(Deg C)',
            'Bottom Plate(Deg C)',
            'P1(PSI)',
            'P2(PSI)',
            'Mass Flow rate(g/min)',
            'Heater Energy (kJ)'
        ]
        
        df = pd.read_csv(filepath, sep=r'\s+', names=column_names, engine='python')
    
    # Clean column names (remove extra spaces)
    df.columns = df.columns.str.strip()
    
    return df, metadata


def plot_variables(df, var1, var2, title=None, save_path=None, metadata=None):
    """
    Plot two variables against time.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the data
    var1 : str
        First variable to plot (y-axis left)
    var2 : str
        Second variable to plot (y-axis right)
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save the plot
    metadata : dict, optional
        Metadata to display in the plot
    """
    # Get time column (should be first column)
    time_col = df.columns[0]
    
    # Create figure with two y-axes
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Plot first variable
    color = 'tab:blue'
    ax1.set_xlabel(time_col, fontsize=12)
    ax1.set_ylabel(var1, color=color, fontsize=12)
    ax1.plot(df[time_col], df[var1], color=color, linewidth=2, label=var1)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    
    # Create second y-axis
    ax2 = ax1.twinx()
    color = 'tab:orange'
    ax2.set_ylabel(var2, color=color, fontsize=12)
    ax2.plot(df[time_col], df[var2], color=color, linewidth=2, label=var2)
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Add title
    if title:
        plt.title(title, fontsize=14, fontweight='bold')
    else:
        plt.title(f'{var1} and {var2} vs {time_col}', fontsize=14)
    
    # Add metadata text if available
    if metadata:
        metadata_text = '\n'.join([f"{key}: {value}" for key, value in metadata.items()])
        plt.text(0.02, 0.98, metadata_text, transform=fig.transFigure,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Add legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
    
    plt.tight_layout()
    
    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()


def list_available_columns(df):
    """
    Print available columns in the DataFrame.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to inspect
    """
    print("\nAvailable columns:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")


# Example usage
if __name__ == "__main__":
    # Define paths
    data_dir = Path(__file__).parent / "rawdata"
    
    # Example 1: Load a file with headers
    file1 = data_dir / "Lab 2 - Part 2d.txt"
    if file1.exists():
        print(f"Loading {file1.name}...")
        df1, metadata1 = load_data_file(file1)
        
        print(f"\nDataFrame shape: {df1.shape}")
        print(f"Metadata: {metadata1}")
        list_available_columns(df1)
        
        print("\nFirst few rows:")
        print(df1.head())

        # Plot temperature and pressure
        plot_variables(
            df1, 
            'T1(Deg C)', 
            'Mass Flow rate(g/min)',
            title=f'Temperature and Heater Energy vs Time - {file1.name}',
            metadata=metadata1
        )