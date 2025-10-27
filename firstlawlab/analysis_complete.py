"""
Complete Analysis Script for CHE260 First Law of Thermodynamics Lab
Consolidates all data processing, integration, and visualization.
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
import matplotlib.colors as mcolors
from pathlib import Path
from uncertainties import ufloat
from uncertainties.umath import log as ulog
from scipy.stats import linregress
from scipy import integrate
from scipy.ndimage import gaussian_filter1d
import os

# ============================================================================
# TEXT STYLING HELPERS FOR READABLE LABELS
# ============================================================================

def _auto_text_color(color):
    """Auto-select black or white text based on background luminance."""
    r, g, b = mcolors.to_rgb(color)
    # Relative luminance formula (WCAG)
    L = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return 'white' if L < 0.6 else 'black'

TEXT_STROKE = [pe.withStroke(linewidth=3, foreground='white')]
BBOX = dict(facecolor='white', alpha=0.85, boxstyle='round,pad=0.4', edgecolor='gray', linewidth=1)


# ============================================================================
# DATA LOADING AND PROCESSING
# ============================================================================

def load_data_file(filepath):
    """
    Load a data file from firstlawlab/rawdata directory.
    
    Handles two formats:
    1. Files with headers (ambient conditions + column names)
    2. Files without headers (only numeric data)
    
    NOTE: This function is designed for potential future expansion or batch processing.
    It is currently used by analyze_part1() to extract metadata (ambient pressure).
    For Part 2 analysis, direct pd.read_csv() is used in analyze_part2_basic() for efficiency.
    
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
        
        # Define standard column names
        column_names = [
            'Time(s)',
            'T1(Deg C)',
            'Wall Temp(Deg C)',
            'Top Plate(Deg C)',
            'Bottom Plate(Deg C)',
            'P1(PSI)',
            'P2(PSI)',
            'Mass Flowrate(g/min)',
            'Heater Energy (kJ)'
        ]
        
        # Read data starting from line 4 (skip first 4 lines: 2 metadata + 2 header lines)
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
            'Mass Flowrate(g/min)',
            'Heater Energy (kJ)'
        ]
        
        df = pd.read_csv(filepath, sep=r'\s+', names=column_names, engine='python')
    
    # Clean column names (remove extra spaces)
    df.columns = df.columns.str.strip()
    
    return df, metadata


def list_available_columns(df):
    """Print available columns in the DataFrame."""
    print("\nAvailable columns:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")


# ============================================================================
# MASS FLOW INTEGRATION
# ============================================================================

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
    mass_flow_rate = df['Mass Flowrate(g/min)'].values
    
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


# ============================================================================
# PART 1: MASS ANALYSIS
# ============================================================================

def analyze_part1(trial):
    """
    Analyze Part 1 data: mass balance and flow rates.
    
    Parameters:
    -----------
    trial : str
        Trial identifier ('a', 'b', 'c', or 'd')
        
    Returns:
    --------
    dict : Results dictionary with calculated values
    """
    filepath1 = f"rawdata/Lab 2 - Part 1{trial}.txt"
    
    # Try to load with metadata first
    try:
        df1, metadata = load_data_file(filepath1)
        ambient_inHg = metadata.get('ambient_pressure_inHg', 29.9)
    except:
        # Fallback to direct read
        df1 = pd.read_csv(filepath1, sep="\t", skiprows=2, header=0)
        ambient_inHg = 29.9
    
    # Get average initial + final values
    P1_gauge_initial = df1.head(10)['P1(PSI)'].mean()
    T1_initial = df1.head(10)['T1(Deg C)'].mean()
    P1_gauge_final = df1.tail(10)['P1(PSI)'].mean()
    T1_final = df1.tail(10)['T1(Deg C)'].mean()
    
    # Convert to Kelvin
    T1_initial_K = T1_initial + 273.15
    T1_final_K = T1_final + 273.15
    
    # Adjust for ambient pressure
    ambient_pressure_psi = ambient_inHg * 0.491154
    P1_abs_initial = P1_gauge_initial + ambient_pressure_psi
    P1_abs_final = P1_gauge_final + ambient_pressure_psi
    
    # Mass flow rate is in g/min, time is in seconds
    time_steps = df1['Time(s)'].to_numpy()
    mass_flow_rates = df1['Mass Flowrate(g/min)'].to_numpy()
    time_minutes = time_steps / 60.0
    
    # Integrate mass flow rate over time using trapezoidal rule
    positive_flow_mask = mass_flow_rates > 0
    if np.any(positive_flow_mask):
        positive_times = time_minutes[positive_flow_mask]
        positive_flows = mass_flow_rates[positive_flow_mask]
        m_added_grams = np.trapezoid(positive_flows, positive_times)
    else:
        m_added_grams = 0
        
    m_addedtotank = m_added_grams / 1000.0
    
    # Calculate mass left in tank using ideal gas law
    ratio = (P1_abs_final * T1_initial_K) / (P1_abs_initial * T1_final_K)
    
    # Guard against ratio-1 ≈ 0
    if abs(ratio - 1) < 1e-6:
        raise ValueError(f"Pressure/temperature ratio too close to 1 (ratio={ratio:.6f}); check windows or data quality.")
    
    m_lefttank = m_addedtotank * (1 + 1/(ratio - 1))
    
    results = {
        'trial': trial,
        'P1_abs_initial': P1_abs_initial,
        'T1_initial': T1_initial,
        'T1_initial_K': T1_initial_K,
        'P1_abs_final': P1_abs_final,
        'T1_final': T1_final,
        'T1_final_K': T1_final_K,
        'm_added_grams': m_added_grams,
        'm_added_kg': m_addedtotank,
        'm_lefttank_kg': m_lefttank,
        'm_escaped_kg': m_addedtotank - m_lefttank,
        'pressure_ratio': ratio,
        'time_steps': time_steps,
        'mass_flow_rates': mass_flow_rates,
        'time_minutes': time_minutes
    }
    
    return results


def print_part1_results(results):
    """Print Part 1 analysis results."""
    trial = results['trial']
    print(f"\n{'='*70}")
    print(f"PART 1 ANALYSIS - Trial {trial}")
    print(f"{'='*70}")
    print(f"Initial conditions (avg of first 10 points):")
    print(f"  P1 (abs): {results['P1_abs_initial']:.3f} PSI")
    print(f"  T1: {results['T1_initial']:.2f}°C ({results['T1_initial_K']:.2f}K)")
    print(f"\nFinal conditions (avg of last 10 points):")
    print(f"  P1 (abs): {results['P1_abs_final']:.3f} PSI")
    print(f"  T1: {results['T1_final']:.2f}°C ({results['T1_final_K']:.2f}K)")
    print(f"\nMass Balance:")
    print(f"  Mass added to tank: {results['m_added_kg']:.4f} kg ({results['m_added_grams']:.2f} g)")
    print(f"  Mass left in tank: {results['m_lefttank_kg']:.4f} kg")
    print(f"  Mass that escaped: {results['m_escaped_kg']:.4f} kg")
    print(f"\nPressure ratio (P2*T1)/(P1*T2): {results['pressure_ratio']:.4f}")


def plot_part1_mass_flow(results):
    """Create integration plot for Part 1 mass flow with enhanced readability."""
    trial = results['trial']
    time_steps = results['time_steps']
    mass_flow_rates = results['mass_flow_rates']
    time_minutes = results['time_minutes']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True, facecolor='white')
    
    # Plot mass flow rate with enhanced styling
    ax1.plot(time_steps, mass_flow_rates, color='#2980b9', linewidth=2.5, label='Mass Flow Rate')
    ax1.fill_between(time_steps, mass_flow_rates, 0, where=(mass_flow_rates > 0),
                     alpha=0.25, color='#2980b9')
    ax1.axhline(y=0, color='#555', linestyle='--', alpha=0.6, linewidth=1)
    ax1.set_ylabel('Mass Flow Rate (g/min)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Trial {trial}: Mass Flow Analysis', fontsize=13, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(fontsize=11, loc='upper right', framealpha=0.95)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Calculate and plot cumulative mass
    cumulative_mass = np.cumsum(np.where(mass_flow_rates > 0, mass_flow_rates * np.gradient(time_minutes), 0.0))
    ax2.plot(time_steps, cumulative_mass, color='#e74c3c', linewidth=2.5, label='Cumulative Mass')
    ax2.fill_between(time_steps, cumulative_mass, alpha=0.2, color='#e74c3c')
    ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cumulative Mass Added (g)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(fontsize=11, loc='upper left', framealpha=0.95)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'graphs/part1_trial_{trial}_mass_flow.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: graphs/part1_trial_{trial}_mass_flow.png")
    plt.close()


# ============================================================================
# PART 2: HEAT LOSS AND SPECIFIC HEAT ANALYSIS
# ============================================================================

def analyze_part2_basic(trial):
    """
    Analyze Part 2 data: heat loss and steady-state analysis.
    
    Parameters:
    -----------
    trial : str
        Trial identifier ('a', 'b', 'c', or 'd')
        
    Returns:
    --------
    dict : Results dictionary with calculated values
    """
    filepath2 = f"rawdata/Lab 2 - Part 2{trial}.txt"
    df2 = pd.read_csv(filepath2, sep="\t", skiprows=2, header=0)
    
    # Convert Heater Energy from kJ to J
    df2['Heater Energy (J)'] = df2['Heater Energy (kJ)'] * 1000
    
    # Identify steady-state (last 20% of data)
    steady_state_df = df2.tail(int(len(df2) * 0.2))
    
    # Linear regression on steady-state data
    time_steady = steady_state_df['Time(s)'].to_numpy()
    energy_steady = steady_state_df['Heater Energy (J)'].to_numpy()
    slope, intercept, r_value, p_value, std_err = linregress(time_steady, energy_steady)
    
    heat_loss_power = slope
    avg_steady_temp = steady_state_df['T1(Deg C)'].mean()
    
    results = {
        'trial': trial,
        'df': df2,
        'steady_state_start': time_steady[0],
        'steady_state_end': time_steady[-1],
        'avg_steady_temp': avg_steady_temp,
        'heat_loss_power': heat_loss_power,
        'heat_loss_power_std_err': std_err
    }
    
    return results


def print_part2_basic_results(results):
    """Print Part 2 basic analysis results."""
    trial = results['trial']
    print(f"\n{'='*70}")
    print(f"PART 2 BASIC ANALYSIS - Trial {trial}")
    print(f"{'='*70}")
    print(f"Steady-state period: t={results['steady_state_start']:.1f}s to t={results['steady_state_end']:.1f}s")
    print(f"Average temperature during steady state: {results['avg_steady_temp']:.2f}°C")
    print(f"Calculated Heat Loss Rate (Input Power): {results['heat_loss_power']:.2f} ± {results['heat_loss_power_std_err']:.3f} W")


# ============================================================================
# HEAT LOSS SPLIT: ACRYLIC WALLS VS PLATES
# ============================================================================

def conduction_split(df_ss, k=0.185, L=0.28575, r2=0.1016/2, thickness=0.009525, T_amb=None):
    """
    Split heat loss between acrylic wall conduction and residual plate losses.
    
    Parameters:
    -----------
    df_ss : pd.DataFrame
        Steady-state data (last 20%)
    k : float
        Thermal conductivity of acrylic (W/m·K), default 0.185
    L : float
        Height of cylinder (m), default 0.28575
    r2 : float
        Outer radius (m), default 0.1016/2 = 0.0508
    thickness : float
        Acrylic thickness (m), default 0.009525
    T_amb : float, optional
        Ambient temperature (°C) for outer wall. If None, uses first 30 points average.
        
    Returns:
    --------
    dict : {
        'Qdot_acrylic': heat loss through acrylic walls (W),
        'Qdot_acrylic_u': uncertainty in Qdot_acrylic (W),
        'T_gas': average gas temperature during steady state (°C),
        'T_wall_out': average outer wall temperature (°C),
        'r1': inner radius (m),
        'r2': outer radius (m),
        'ln_ratio': ln(r2/r1)
    }
    """
    r1 = r2 - thickness
    T_gas = df_ss['T1(Deg C)'].mean()
    
    # Prefer measured wall temp if available
    if 'Wall Temp(Deg C)' in df_ss.columns:
        T_wall_out = df_ss['Wall Temp(Deg C)'].mean()
    else:
        # Fallback: use first 30 points or ambient
        if T_amb is not None:
            T_wall_out = T_amb
        else:
            T_wall_out = df_ss['T1(Deg C)'].head(30).mean() - 5  # rough estimate
    
    dT_wall = T_gas - T_wall_out
    ln_ratio = math.log(r2 / r1)
    
    # Heat transfer through acrylic (cylindrical conduction)
    Qdot_acrylic = 2 * math.pi * k * L * (dT_wall / ln_ratio)
    
    # Simple uncertainty estimate (propagate thermal conductivity and geometry)
    sigma_k = 0.015  # ±0.015 W/m·K
    sigma_L = 0.0001  # ±0.1 mm
    sigma_dT = math.sqrt(0.05**2 + 0.05**2)  # ±0.05°C for each temp sensor
    
    # Partial derivatives
    partial_k = 2 * math.pi * L * (dT_wall / ln_ratio)
    partial_L = 2 * math.pi * k * (dT_wall / ln_ratio)
    partial_dT = 2 * math.pi * k * L / ln_ratio
    
    sigma_Qdot = math.sqrt(
        (partial_k * sigma_k)**2 +
        (partial_L * sigma_L)**2 +
        (partial_dT * sigma_dT)**2
    )
    
    return {
        'Qdot_acrylic': Qdot_acrylic,
        'Qdot_acrylic_u': sigma_Qdot,
        'T_gas': T_gas,
        'T_wall_out': T_wall_out,
        'r1': r1,
        'r2': r2,
        'ln_ratio': ln_ratio
    }


# ============================================================================
# HEAT LOSS SPLIT WITH UNCERTAINTY PROPAGATION
# ============================================================================

def analyze_heat_loss_split(results_p2, trial):
    """
    Analyze heat loss split (walls vs plates) with uncertainty propagation.
    
    Parameters:
    -----------
    results_p2 : dict
        Part 2 analysis results (from analyze_part2_basic)
    trial : str
        Trial identifier
        
    Returns:
    --------
    dict : {
        'trial': trial,
        'Qdot_total': total heat loss (W),
        'Qdot_total_u': uncertainty (W),
        'Qdot_acrylic': wall loss (W),
        'Qdot_acrylic_u': uncertainty (W),
        'Qdot_plates': residual plate loss (W),
        'Qdot_plates_u': uncertainty (W),
        'T_gas_ss': steady-state gas temp (°C)
    }
    """
    df = results_p2['df']
    ss = df.tail(int(len(df) * 0.2))
    
    Q_total = results_p2['heat_loss_power']
    Q_total_u = results_p2['heat_loss_power_std_err']
    
    # Get wall split
    wall_split = conduction_split(ss)
    Q_acrylic = wall_split['Qdot_acrylic']
    Q_acrylic_u = wall_split['Qdot_acrylic_u']
    
    # Use ufloat for uncertainty propagation
    Q_total_uf = ufloat(Q_total, Q_total_u)
    Q_acrylic_uf = ufloat(Q_acrylic, Q_acrylic_u)
    Q_plates_uf = Q_total_uf - Q_acrylic_uf
    
    return {
        'trial': trial,
        'Qdot_total': Q_total,
        'Qdot_total_u': Q_total_u,
        'Qdot_acrylic': Q_acrylic,
        'Qdot_acrylic_u': Q_acrylic_u,
        'Qdot_plates': float(Q_plates_uf.nominal_value),
        'Qdot_plates_u': float(Q_plates_uf.std_dev),
        'T_gas_ss': wall_split['T_gas']
    }


# ============================================================================
# HEAT TRANSFER AND PLATE ANALYSIS
# ============================================================================

def calculate_heat_transfer_with_uncertainty(inner_wall_temp, net_heat_transfer, net_heat_transfer_uncertainty):
    """
    Calculate heat transfer through walls and plates with uncertainty propagation.
    
    Parameters:
    -----------
    inner_wall_temp : float
        Inner wall temperature (°C)
    net_heat_transfer : float
        Net heat transfer rate (W)
    net_heat_transfer_uncertainty : float
        Uncertainty in net heat transfer (W)
        
    Returns:
    --------
    tuple : (heat_walls, heat_walls_unc, heat_plates, heat_plates_unc)
    """
    # Constants with uncertainties
    k_acrylic = 0.185
    k_acrylic_uncertainty = 0.015
    l = 0.2858
    l_uncertainty = 0.000001
    r1 = 0.09208
    r1_uncertainty = 0.0006
    r2 = 0.1016
    r2_uncertainty = 0.001
    
    outer_wall_temp = 18.5
    outer_wall_temp_uncertainty = 0.05
    inner_wall_temp_uncertainty = 0.05
    
    # Calculate temperature difference and uncertainty
    delta_T = inner_wall_temp - outer_wall_temp
    delta_T_uncertainty = math.sqrt(inner_wall_temp_uncertainty**2 + outer_wall_temp_uncertainty**2)
    
    # Heat transfer through walls
    heat_transfer_walls = 2 * k_acrylic * math.pi * l * (delta_T / math.log(r2 / r1))
    
    # Partial derivatives for uncertainty propagation
    partial_k = 2 * math.pi * l * (delta_T / math.log(r2 / r1))
    partial_l = 2 * k_acrylic * math.pi * (delta_T / math.log(r2 / r1))
    partial_delta_T = 2 * k_acrylic * math.pi * l / math.log(r2 / r1)
    partial_r1 = -2 * k_acrylic * math.pi * l * delta_T / (r1 * (math.log(r2 / r1)**2))
    partial_r2 = 2 * k_acrylic * math.pi * l * delta_T / (r2 * (math.log(r2 / r1)**2))
    
    heat_transfer_walls_uncertainty = math.sqrt(
        (partial_k * k_acrylic_uncertainty)**2 +
        (partial_l * l_uncertainty)**2 +
        (partial_delta_T * delta_T_uncertainty)**2 +
        (partial_r1 * r1_uncertainty)**2 +
        (partial_r2 * r2_uncertainty)**2
    )
    
    # Heat transfer through plates
    heat_transfer_plates = net_heat_transfer - heat_transfer_walls
    heat_transfer_plates_uncertainty = math.sqrt(
        net_heat_transfer_uncertainty**2 + heat_transfer_walls_uncertainty**2
    )
    
    return heat_transfer_walls, heat_transfer_walls_uncertainty, heat_transfer_plates, heat_transfer_plates_uncertainty


# ============================================================================
# SPECIFIC HEAT ANALYSIS
# ============================================================================

def estimate_cv(Q_in_J, Q_loss_W, dt_s, m_kg, dT_K):
    """
    Estimate specific heat capacity at constant volume.
    
    Parameters:
    -----------
    Q_in_J : float
        Heat input from heater (J)
    Q_loss_W : float
        Heat loss rate (W) during the measurement period
    dt_s : float
        Duration of measurement (s)
    m_kg : float
        Total mass of air (kg)
    dT_K : float
        Temperature change (K)
        
    Returns:
    --------
    float : c_v estimate (J/kg·K)
    """
    if dT_K < 0.1:
        return float('nan')
    
    Q_net_J = Q_in_J - Q_loss_W * dt_s
    c_v = Q_net_J / (m_kg * dT_K)
    return c_v


def estimate_cv_with_uncertainty(Q_in_J, Q_in_u, Q_loss_W, Q_loss_u, dt_s, m_kg, m_u, dT_K, dT_u):
    """
    Estimate c_v with uncertainty propagation.
    
    Parameters:
    -----------
    Q_in_J : float
        Heat input (J)
    Q_in_u : float
        Heat input uncertainty (J)
    Q_loss_W : float
        Heat loss rate (W)
    Q_loss_u : float
        Heat loss rate uncertainty (W)
    dt_s : float
        Duration (s)
    m_kg : float
        Total mass (kg)
    m_u : float
        Mass uncertainty (kg)
    dT_K : float
        Temperature change (K)
    dT_u : float
        Temperature uncertainty (K)
        
    Returns:
    --------
    tuple : (c_v, sigma_c_v)
    """
    if dT_K < 0.1:
        return float('nan'), float('nan')
    
    Q_net_J = Q_in_J - Q_loss_W * dt_s
    c_v = Q_net_J / (m_kg * dT_K)
    
    # Uncertainty components
    dQ_net_dQ_in = 1.0
    dQ_net_dQ_loss = -dt_s
    dc_v_dQ_net = 1.0 / (m_kg * dT_K)
    dc_v_dm = -Q_net_J / (m_kg**2 * dT_K)
    dc_v_dT = -Q_net_J / (m_kg * dT_K**2)
    
    sigma_Q_net = math.sqrt((dQ_net_dQ_in * Q_in_u)**2 + (dQ_net_dQ_loss * Q_loss_u)**2)
    sigma_c_v = math.sqrt(
        (dc_v_dQ_net * sigma_Q_net)**2 +
        (dc_v_dm * m_u)**2 +
        (dc_v_dT * dT_u)**2
    )
    
    return c_v, sigma_c_v


# ============================================================================
# PROPELLER WORK ANALYSIS
# ============================================================================

def calculate_P2_and_uncertainty(p2, T2, sigma_p2, sigma_T2):
    """
    Calculate propeller power with uncertainty propagation.
    
    Parameters:
    -----------
    p2 : float
        Pressure (kPa)
    T2 : float
        Temperature (°C)
    sigma_p2 : float
        Pressure uncertainty (kPa)
    sigma_T2 : float
        Temperature uncertainty (°C)
        
    Returns:
    --------
    tuple : (P2, sigma_P2) - Power in Watts and its uncertainty
    """
    P1 = 0.07457  # Reference power in Watts
    p1 = 101.325  # Reference pressure in kPa
    T1_Celsius = 25
    T1_Kelvin = T1_Celsius + 273.15
    
    T2_Kelvin = T2 + 273.15
    
    P2 = 9.261 * P1 * (p2 * T1_Kelvin) / (p1 * T2_Kelvin)
    
    partial_p2 = 9.261 * P1 * T1_Kelvin / (p1 * T2_Kelvin)
    partial_T2 = -9.261 * P1 * p2 * T1_Kelvin / (p1 * T2_Kelvin**2)
    
    sigma_P2 = math.sqrt((partial_p2 * sigma_p2)**2 + (partial_T2 * sigma_T2)**2)
    
    return P2, sigma_P2


def fan_power_from_similarity(P1_W=0.001*745.7, n1=4200, D1=0.0635, rho1=1.204, 
                              n2=4200, D2=0.0635, rho2=1.2):
    """
    Calculate propeller power using similarity law.
    
    Power scales as: P ∝ ρ * n³ * D⁵
    
    Parameters:
    -----------
    P1_W : float
        Reference power at reference conditions (W)
    n1 : float
        Reference speed (RPM)
    D1 : float
        Reference diameter (m)
    rho1 : float
        Reference density (kg/m³)
    n2 : float
        Actual speed (RPM)
    D2 : float
        Actual diameter (m)
    rho2 : float
        Actual density (kg/m³)
        
    Returns:
    --------
    float : Estimated power (W)
    """
    P2 = P1_W * (rho2 / rho1) * (n2 / n1)**3 * (D2 / D1)**5
    return P2


def analyze_part3_fan(df_ss, m_total_kg, c_v_avg, T_amb=25):
    """
    Analyze Part 3: propeller work and temperature rise rate.
    
    Parameters:
    -----------
    df_ss : pd.DataFrame
        Steady-state data (last 20%)
    m_total_kg : float
        Total mass of air (kg)
    c_v_avg : float
        Average c_v from Part 2 analysis (J/kg·K)
    T_amb : float
        Ambient temperature (°C)
        
    Returns:
    --------
    dict : {
        'P2_W': propeller power (W),
        'P2_u': uncertainty (W),
        'dT_dt': predicted temperature rise rate (K/s),
        'dT_dt_u': uncertainty (K/s)
    }
    """
    # Get average steady-state conditions
    P2_gauge = df_ss['P2(PSI)'].mean()
    T2 = df_ss['T1(Deg C)'].mean()
    
    # Convert P2 to absolute pressure
    P2_abs_psia = P2_gauge + 14.696  # ~1 atm in psia
    P2_abs_kpa = P2_abs_psia * 6.89476  # convert psia to kPa
    
    # Use similarity law for fan power
    P2_W = fan_power_from_similarity(n2=4200, rho2=1.2)
    
    # Estimate uncertainty (±10% on power)
    P2_u = 0.10 * P2_W
    
    # Temperature rise rate = P / (m * c_v)
    dT_dt = P2_W / (m_total_kg * c_v_avg)
    
    # Propagate uncertainty
    dT_dt_u = dT_dt * math.sqrt(
        (P2_u / P2_W)**2 +
        (0.06 / m_total_kg)**2 +  # ±6% mass uncertainty
        (0.15 / c_v_avg)**2        # ±15% c_v uncertainty
    )
    
    return {
        'P2_W': P2_W,
        'P2_u': P2_u,
        'dT_dt': dT_dt,
        'dT_dt_u': dT_dt_u,
        'T2_measured': T2,
        'P2_measured': P2_abs_kpa
    }


def compute_air_mass_g(t_initial, ambient_pressure=101.325, gas_constant=8.314, 
                       molar_mass=28.97, v_cylinder=0.0336):
    """
    Compute air mass using ideal gas law.
    
    Parameters:
    -----------
    t_initial : float
        Initial temperature (°C)
    ambient_pressure : float
        Ambient pressure (kPa), default 101.325
    gas_constant : float
        Gas constant (J/mol·K), default 8.314
    molar_mass : float
        Molar mass of air (g/mol), default 28.97
    v_cylinder : float
        Cylinder volume (m³), default 0.0336
        
    Returns:
    --------
    float : Air mass in grams
    """
    ambient_pressure_Pa = ambient_pressure * 1000
    t_initial_K = t_initial + 273.15
    air_mass = molar_mass * (ambient_pressure_Pa * v_cylinder) / (gas_constant * t_initial_K)
    return air_mass


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_temperature_and_heater_energy(data_key, df, sigma=10, cutoff_index=-1, 
                                       calc_start=None, calc_end=None, mass_g=None):
    """
    Plot temperature and heater energy with optional analysis region.
    
    Parameters:
    -----------
    data_key : str
        Identifier for the trial
    df : pd.DataFrame
        DataFrame containing the data
    sigma : float
        Gaussian smoothing sigma
    cutoff_index : int
        Index to cut off data
    calc_start : float
        Start time for analysis region
    calc_end : float
        End time for analysis region
    mass_g : float
        Mass of air in grams
    """
    # Apply Gaussian smoothing
    temp = gaussian_filter1d(np.array(df['T1(Deg C)'][:cutoff_index]), sigma=sigma)
    time = np.array(df['Time(s)'][:cutoff_index])
    heater_energy = gaussian_filter1d(np.array(df['Heater Energy (kJ)'][:cutoff_index]), sigma=sigma)
    
    initial_air_mass = compute_air_mass_g(temp[0])
    mass_added = mass_g
    total_mass = mass_g + initial_air_mass
    
    # Create subplots
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 8))
    
    # Temperature plot
    axs[0].plot(time, temp, label='T1 (°C)', linewidth=2)
    axs[0].set_ylabel('Temperature (°C)', fontsize=11)
    axs[0].legend(fontsize=10)
    axs[0].set_title(f'Temperature and Heater Energy Analysis - Trial {data_key[-1]}', 
                     fontsize=12, fontweight='bold')
    axs[0].grid(True, alpha=0.3)
    
    # Heater energy plot
    axs[1].plot(time, heater_energy, label='Heater Energy (kJ)', color='r', linewidth=2)
    axs[1].set_xlabel('Time (s)', fontsize=11)
    axs[1].set_ylabel('Heater Energy (kJ)', fontsize=11)
    axs[1].legend(fontsize=10)
    axs[1].grid(True, alpha=0.3)
    
    if calc_start is not None and calc_end is not None and mass_g is not None:
        # Find indices for analysis region
        start_index = np.searchsorted(time, calc_start)
        end_index = np.searchsorted(time, calc_end)
        
        init_temp = temp[start_index]
        end_temp = temp[end_index]
        init_energy = heater_energy[start_index]
        end_energy = heater_energy[end_index]
        
        # Mark analysis region
        axs[0].axhline(y=init_temp, color='green', linestyle='--', alpha=0.7)
        axs[0].axhline(y=end_temp, color='blue', linestyle='--', alpha=0.7)
        axs[1].axhline(y=init_energy, color='green', linestyle='--', alpha=0.7)
        axs[1].axhline(y=end_energy, color='blue', linestyle='--', alpha=0.7)
        
        axs[0].axvline(x=calc_start, color='gray', linestyle='--', alpha=0.7)
        axs[0].axvline(x=calc_end, color='gray', linestyle='--', alpha=0.7)
        axs[1].axvline(x=calc_start, color='gray', linestyle='--', alpha=0.7)
        axs[1].axvline(x=calc_end, color='gray', linestyle='--', alpha=0.7)
        
        # Calculate c_v
        delta_q = (end_energy - init_energy) * 1000  # Convert kJ to J
        delta_t = end_temp - init_temp
        c_v = delta_q / (total_mass * delta_t) if delta_t != 0 else float('inf')
        
        # Uncertainty calculation
        delta_energy_uncertainty = 0.005
        delta_temp_uncertainty = 0.05
        mass_uncertainty = 0.06 if data_key[-1] in ['a', 'c'] else 0.08
        
        delta_energy = abs(end_energy - init_energy) * 1000
        delta_temp = abs(end_temp - init_temp)
        relative_energy_uncertainty = delta_energy_uncertainty / delta_energy if delta_energy != 0 else 0
        relative_temp_uncertainty = delta_temp_uncertainty / delta_temp if delta_temp != 0 else 0
        relative_mass_uncertainty = mass_uncertainty / total_mass
        
        delta_c_v = c_v * np.sqrt(
            relative_energy_uncertainty**2 + relative_temp_uncertainty**2 + relative_mass_uncertainty**2
        )
        
        # Annotate c_v with readable background box and text stroke
        axs[0].annotate(f'c_v: {c_v:.2f} ± {delta_c_v:.2f} J/g°C', 
                        xy=(calc_end, end_temp), xytext=(calc_end + 5, end_temp - 3.5),
                        arrowprops=dict(arrowstyle='->', color='black', lw=2), 
                        fontsize=12, fontweight='bold', bbox=BBOX, path_effects=TEXT_STROKE, zorder=10)
        axs[0].annotate(f'(m_initial + m_added = {initial_air_mass:.2f} + {mass_added:.2f} = {total_mass:.2f}g)', 
                        xy=(calc_end + 5, end_temp - 5.3), fontsize=10, bbox=BBOX, zorder=10)
    
    plt.tight_layout()
    plt.savefig(f'graphs/part2_trial_{data_key[-1]}_temp_heater_energy.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: graphs/part2_trial_{data_key[-1]}_temp_heater_energy.png")
    plt.close()


def plot_heater_temp_maintain(trial, df, start_idx, end_idx, sigma=10):
    """
    Plot steady-state temperature maintenance with heat rate analysis.
    
    Parameters:
    -----------
    trial : str
        Trial identifier
    df : pd.DataFrame
        DataFrame containing the data
    start_idx : int
        Start index for analysis
    end_idx : int
        End index for analysis
    sigma : float
        Gaussian smoothing sigma
    """
    # Extract and smooth data
    temp = gaussian_filter1d(np.array(df['T1(Deg C)'].iloc[start_idx:end_idx]), sigma=sigma)
    time = np.array(df['Time(s)'].iloc[start_idx:end_idx])
    heater_energy = gaussian_filter1d(np.array(df['Heater Energy (kJ)'].iloc[start_idx:end_idx]), sigma=sigma)
    
    # Calculate heat rate
    time_diff = np.diff(time)
    energy_diff = np.diff(heater_energy * 1000)
    heat_rate = energy_diff / time_diff
    time_heat_rate = time[:-1] + time_diff / 2
    
    # Create subplots
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    
    # Temperature plot
    axs[0].plot(time, temp, label='Temperature (°C)', linewidth=2)
    axs[0].set_ylabel('Temperature (°C)', fontsize=11)
    axs[0].legend(fontsize=10)
    axs[0].set_title(f'Steady-State Temperature Maintenance - Trial {trial}', 
                     fontsize=12, fontweight='bold')
    axs[0].grid(True, alpha=0.3)
    
    # Heater energy plot
    axs[1].plot(time, heater_energy, label='Heater Energy (kJ)', color='r', linewidth=2)
    axs[1].set_ylabel('Heater Energy (kJ)', fontsize=11)
    axs[1].legend(fontsize=10)
    axs[1].grid(True, alpha=0.3)
    
    # Heat rate plot
    axs[2].plot(time_heat_rate, heat_rate, label='Heat Rate (W)', color='#27ae60', linewidth=2.5)
    axs[2].set_xlabel('Time (s)', fontsize=11)
    axs[2].set_ylabel('Heat Rate (W)', fontsize=11)
    axs[2].legend(fontsize=10, loc='upper right')
    axs[2].grid(True, alpha=0.3)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    
    # Average heat rate annotation with readable box
    avg_heat_rate = np.average(heat_rate)
    axs[2].axhline(y=avg_heat_rate, color='#888', linestyle='--', alpha=0.6, linewidth=2)
    axs[2].annotate(f'Average: {avg_heat_rate:.2f} W',
                    xy=(time_heat_rate[int(len(time_heat_rate)*0.5)], avg_heat_rate),
                    xytext=(time_heat_rate[int(len(time_heat_rate)*0.6)], avg_heat_rate + 8),
                    fontsize=11, fontweight='bold', bbox=BBOX, path_effects=TEXT_STROKE,
                    arrowprops=dict(arrowstyle='-|>', color='gray', lw=1.5), zorder=10)
    
    plt.tight_layout()
    plt.savefig(f'graphs/part2_trial_{trial}_heater_maintenance.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: graphs/part2_trial_{trial}_heater_maintenance.png")
    plt.close()


def plot_heat_loss_breakdown(trial, heat_split_results):
    """
    Plot heat loss breakdown (acrylic walls vs plates) as a bar chart.
    Enhanced with readable text labels and professional styling.
    
    Parameters:
    -----------
    trial : str
        Trial identifier
    heat_split_results : dict
        Results from analyze_heat_loss_split()
    """
    Q_total = heat_split_results['Qdot_total']
    Q_acrylic = heat_split_results['Qdot_acrylic']
    Q_plates = heat_split_results['Qdot_plates']
    Q_acrylic_u = heat_split_results['Qdot_acrylic_u']
    Q_plates_u = heat_split_results['Qdot_plates_u']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(9, 6.5), facecolor='white')
    
    # Bar positions
    x_pos = np.arange(2)
    values = [Q_acrylic, Q_plates]
    uncertainties = [Q_acrylic_u, Q_plates_u]
    labels = ['Acrylic Wall\nConduction', 'Plate Loss\n(Residual)']
    colors = ['#FF6B6B', '#4ECDC4']
    
    # Create bars with darker error bars for visibility
    bars = ax.bar(x_pos, values, yerr=uncertainties, capsize=12, color=colors, 
                   edgecolor='#333', linewidth=2, alpha=0.85, 
                   error_kw={'linewidth': 2.5, 'ecolor': '#333333'})
    
    # Add value labels on bars with auto-contrasting text color and outline
    max_height = max(values)
    for bar, val, unc, fill_color in zip(bars, values, uncertainties, colors):
        height = bar.get_height()
        txt_color = _auto_text_color(fill_color)
        label_text = f'{val:.2f} ± {unc:.2f} W'
        
        # Position label well above bar + uncertainty (at 12% of max_height above top)
        y_pos = height + unc + max_height * 0.03
        ax.text(bar.get_x() + bar.get_width()/2., y_pos,
                label_text,
                ha='center', va='bottom', fontsize=12, fontweight='bold',
                color=txt_color, path_effects=TEXT_STROKE, zorder=10)
    
    # Add total line with enhanced visibility
    ax.axhline(y=Q_total, color='#555', linestyle='--', linewidth=2.5, 
               label=f'Total: {Q_total:.2f} W', alpha=0.7, zorder=2)
    
    # Formatting
    ax.set_ylabel('Heat Loss Rate (W)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title(f'Heat Loss Breakdown - Trial {trial}', fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=11, fontweight='bold')
    ax.set_ylim(0, Q_total * 1.35)
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax.legend(fontsize=11, loc='upper right', framealpha=0.95)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'graphs/part2_trial_{trial}_loss_breakdown.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: graphs/part2_trial_{trial}_loss_breakdown.png")
    plt.close()
    plt.close()


# ============================================================================
# MAIN ANALYSIS ORCHESTRATION
# ============================================================================

def run_complete_analysis():
    """Run complete analysis for all trials."""
    
    # Create graphs directory if it doesn't exist
    os.makedirs('graphs', exist_ok=True)
    
    trials = ['a', 'b', 'c', 'd']
    part1_results = {}
    part2_results = {}
    part2_heat_split = {}
    part2_cv_results = {}
    part3_fan_results = {}
    
    # Lists to accumulate rows for CSV export
    p1_rows = []
    p2_heat_rows = []
    cv_rows = []
    fan_rows = []
    
    print("\n" + "="*70)
    print("CHE260 FIRST LAW OF THERMODYNAMICS LAB - COMPLETE ANALYSIS")
    print("="*70)
    
    # ========================================================================
    # PART 1: MASS ANALYSIS FOR ALL TRIALS
    # ========================================================================
    print("\n" + "="*70)
    print("PART 1: MASS BALANCE ANALYSIS")
    print("="*70)
    
    for trial in trials:
        try:
            results = analyze_part1(trial)
            part1_results[trial] = results
            print_part1_results(results)
            plot_part1_mass_flow(results)
            
            # Accumulate Part 1 results for CSV
            p1_rows.append({
                'Trial': trial,
                'P1_abs_initial_psia': results['P1_abs_initial'],
                'T1_initial_C': results['T1_initial'],
                'P1_abs_final_psia': results['P1_abs_final'],
                'T1_final_C': results['T1_final'],
                'm_added_g': results['m_added_grams'],
                'm_added_kg': results['m_added_kg'],
                'm_left_tank_kg': results['m_lefttank_kg'],
                'm_escaped_kg': results['m_escaped_kg'],
                'pressure_ratio': results['pressure_ratio']
            })
        except Exception as e:
            print(f"Error analyzing Part 1 Trial {trial}: {e}")
    
    # ========================================================================
    # PART 2: HEAT AND SPECIFIC HEAT ANALYSIS FOR ALL TRIALS
    # ========================================================================
    print("\n" + "="*70)
    print("PART 2: HEAT LOSS AND SPECIFIC HEAT ANALYSIS")
    print("="*70)
    
    for trial in trials:
        try:
            results = analyze_part2_basic(trial)
            part2_results[trial] = results
            print_part2_basic_results(results)
            
            # Analyze heat loss split (walls vs plates)
            split_results = analyze_heat_loss_split(results, trial)
            part2_heat_split[trial] = split_results
            
            print(f"\n  Heat Loss Breakdown (Trial {trial}):")
            print(f"    Total heat loss:     {split_results['Qdot_total']:.2f} ± {split_results['Qdot_total_u']:.2f} W")
            print(f"    Acrylic wall loss:   {split_results['Qdot_acrylic']:.2f} ± {split_results['Qdot_acrylic_u']:.2f} W")
            print(f"    Plate loss (residual): {split_results['Qdot_plates']:.2f} ± {split_results['Qdot_plates_u']:.2f} W")
            
            # Accumulate heat split results
            p2_heat_rows.append({
                'Trial': trial,
                'Q_total_W': split_results['Qdot_total'],
                'Q_total_u_W': split_results['Qdot_total_u'],
                'Q_acrylic_W': split_results['Qdot_acrylic'],
                'Q_acrylic_u_W': split_results['Qdot_acrylic_u'],
                'Q_plates_W': split_results['Qdot_plates'],
                'Q_plates_u_W': split_results['Qdot_plates_u'],
                'T_gas_ss_C': split_results['T_gas_ss']
            })
        except Exception as e:
            print(f"Error analyzing Part 2 Trial {trial}: {e}")
    
    # ========================================================================
    # PLOTTING AND VISUALIZATION
    # ========================================================================
    print("\n" + "="*70)
    print("GENERATING PLOTS")
    print("="*70)
    
    # Part 2 plots (customize with your specific times and masses)
    # calc_start and calc_end are row indices for c_v calculation (need large enough ΔT)
    plot_configs = {
        'a': {'calc_start': 60, 'calc_end': 300, 'mass_g': 23.57, 'cutoff': -1, 'steady_start': 2000, 'steady_end': 4000},
        'b': {'calc_start': 90, 'calc_end': 330, 'mass_g': 43.3, 'cutoff': -1, 'steady_start': 2000, 'steady_end': 4000},
        'c': {'calc_start': 30, 'calc_end': 270, 'mass_g': 23.17, 'cutoff': -1, 'steady_start': 2000, 'steady_end': 4000},
        'd': {'calc_start': 60, 'calc_end': 300, 'mass_g': 40.95, 'cutoff': -1, 'steady_start': 2000, 'steady_end': 4000},
    }
    
    for trial in trials:
        try:
            if trial in part2_results:
                config = plot_configs[trial]
                df = part2_results[trial]['df']
                
                # Temperature and heater energy plot
                plot_temperature_and_heater_energy(
                    f'pt_2_t_{trial}',
                    df,
                    sigma=10,
                    cutoff_index=config['cutoff'],
                    calc_start=config['calc_start'],
                    calc_end=config['calc_end'],
                    mass_g=config['mass_g']
                )
                
                # Steady-state maintenance plot
                plot_heater_temp_maintain(
                    trial,
                    df,
                    start_idx=config['steady_start'],
                    end_idx=config['steady_end'],
                    sigma=10
                )
                
                # Heat loss breakdown bar chart
                if trial in part2_heat_split:
                    plot_heat_loss_breakdown(trial, part2_heat_split[trial])
                
                # Extract c_v for this trial
                initial_air_mass = compute_air_mass_g(df['T1(Deg C)'].iloc[0]) / 1000.0  # convert g to kg
                total_mass_kg = config['mass_g'] / 1000.0 + initial_air_mass
                
                # Estimate c_v from energy balance using row indices
                calc_start = config['calc_start']
                calc_end = config['calc_end']
                
                # Ensure indices are within bounds
                if calc_start < len(df) and calc_end < len(df) and calc_start < calc_end:
                    dt_s = df.iloc[calc_end]['Time(s)'] - df.iloc[calc_start]['Time(s)']
                    Q_in_kJ = df.iloc[calc_end]['Heater Energy (kJ)'] - df.iloc[calc_start]['Heater Energy (kJ)']
                    Q_in_J = Q_in_kJ * 1000
                    T_delta = df.iloc[calc_end]['T1(Deg C)'] - df.iloc[calc_start]['T1(Deg C)']
                    
                    # Estimate heat loss during this period
                    Q_loss_W = part2_heat_split[trial]['Qdot_total']
                    Q_loss_J = Q_loss_W * dt_s
                    
                    c_v_estimate = estimate_cv(Q_in_J, Q_loss_W, dt_s, total_mass_kg, T_delta)
                    
                    # Only save if c_v is valid (not NaN)
                    if not np.isnan(c_v_estimate):
                        c_v_unc = 0.15 * c_v_estimate  # ~±15% uncertainty
                        
                        part2_cv_results[trial] = {
                            'trial': trial,
                            'c_v': c_v_estimate,
                            'c_v_u': c_v_unc,
                            'Q_in_J': Q_in_J,
                            'Q_loss_J': Q_loss_J,
                            'T_delta': T_delta,
                            'm_total_kg': total_mass_kg
                        }
                        
                        cv_rows.append({
                            'Trial': trial,
                            'c_v_J_per_kg_K': c_v_estimate,
                            'c_v_u_J_per_kg_K': c_v_unc,
                            'Q_in_J': Q_in_J,
                            'Q_loss_W': Q_loss_W,
                            'dt_s': dt_s,
                            'T_delta_K': T_delta,
                            'm_total_kg': total_mass_kg
                        })
                        
                        print(f"\n  c_v Estimate (Trial {trial}):")
                        print(f"    c_v = {c_v_estimate:.2f} ± {c_v_unc:.2f} J/kg·K")
                    else:
                        print(f"\n  c_v Estimate (Trial {trial}): Invalid (ΔT too small: {T_delta:.2f}K)")
                else:
                    print(f"\n  c_v Calculation (Trial {trial}): Index out of bounds")
                
        except Exception as e:
            print(f"Error plotting Trial {trial}: {e}")
    
    # ========================================================================
    # PART 3: PROPELLER WORK ANALYSIS
    # ========================================================================
    print("\n" + "="*70)
    print("PART 3: PROPELLER WORK ANALYSIS")
    print("="*70)
    
    for trial in trials:
        try:
            if trial in part2_results and trial in part2_cv_results:
                df = part2_results[trial]['df']
                ss = df.tail(int(len(df) * 0.2))
                m_total = part2_cv_results[trial]['m_total_kg']
                c_v_avg = part2_cv_results[trial]['c_v']
                
                fan_results = analyze_part3_fan(ss, m_total, c_v_avg)
                part3_fan_results[trial] = fan_results
                
                print(f"\nPart 3 Results (Trial {trial}):")
                print(f"  Propeller Power (P2): {fan_results['P2_W']:.4f} ± {fan_results['P2_u']:.4f} W")
                print(f"  Predicted dT/dt: {fan_results['dT_dt']:.4f} ± {fan_results['dT_dt_u']:.4f} K/s")
                
                fan_rows.append({
                    'Trial': trial,
                    'P2_W': fan_results['P2_W'],
                    'P2_u_W': fan_results['P2_u'],
                    'dT_dt_K_per_s': fan_results['dT_dt'],
                    'dT_dt_u_K_per_s': fan_results['dT_dt_u'],
                    'T2_measured_C': fan_results['T2_measured'],
                    'm_total_kg': m_total,
                    'c_v_used_J_per_kg_K': c_v_avg
                })
        except Exception as e:
            print(f"Error analyzing Part 3 Trial {trial}: {e}")
    
    # ========================================================================
    # EXPORT CSV TABLES
    # ========================================================================
    print("\n" + "="*70)
    print("EXPORTING CSV SUMMARY TABLES")
    print("="*70)
    
    try:
        pd.DataFrame(p1_rows).to_csv("part1_summary_all.csv", index=False)
        print("Saved: part1_summary_all.csv")
    except Exception as e:
        print(f"Error saving Part 1 CSV: {e}")
    
    try:
        pd.DataFrame(p2_heat_rows).to_csv("part2_heat_split_all.csv", index=False)
        print("Saved: part2_heat_split_all.csv")
    except Exception as e:
        print(f"Error saving Part 2 heat split CSV: {e}")
    
    try:
        pd.DataFrame(cv_rows).to_csv("part2_cv_all.csv", index=False)
        print("Saved: part2_cv_all.csv")
    except Exception as e:
        print(f"Error saving Part 2 c_v CSV: {e}")
    
    try:
        pd.DataFrame(fan_rows).to_csv("part3_fan_all.csv", index=False)
        print("Saved: part3_fan_all.csv")
    except Exception as e:
        print(f"Error saving Part 3 CSV: {e}")
    
    # ========================================================================
    # SUMMARY REPORT
    # ========================================================================
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  CSV Tables:")
    print("    - part1_summary_all.csv")
    print("    - part2_heat_split_all.csv")
    print("    - part2_cv_all.csv")
    print("    - part3_fan_all.csv")
    print("  Plots:")
    print("    - graphs/part1_trial_[a-d]_mass_flow.png")
    print("    - graphs/part2_trial_[a-d]_temp_heater_energy.png")
    print("    - graphs/part2_trial_[a-d]_heater_maintenance.png")
    print("    - graphs/part2_trial_[a-d]_loss_breakdown.png (wall vs plate split)")
    print("\nNext steps:")
    print("  1. Review CSV tables for data to populate LaTeX tables")
    print("  2. Review generated plots in the 'graphs' directory")
    print("  3. Update plot configurations in run_complete_analysis() as needed")
    print("="*70 + "\n")


if __name__ == "__main__":
    run_complete_analysis()
