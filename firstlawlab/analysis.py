#code to analyze data and generate graphs for lab2 of che260
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties.umath import log as ulog
from scipy.stats import linregress

def main():
    for trial in ["a","b","c","d"]:
        filepath1 = f"rawdata/Lab 2 - Part 1{trial}.txt"
        df1 = pd.read_csv(filepath1, sep="\t", skiprows=2, header=0)
    
        #For Part 1
        #Calculate mass inside tank for each trial
        #m_lefttank = m_addedtotank * (1+1/((P2*T1)/(P1*T2)-1))
        
        #get average initial + final values
        P1_gauge_initial = df1.head(10)['P1(PSI)'].mean()  # Average of first 10 tank pressures
        T1_initial = df1.head(10)['T1(Deg C)'].mean()  # Average of first 10 tank temperatures
        P1_gauge_final = df1.tail(10)['P1(PSI)'].mean()  # Average of last 10 tank pressures
        T1_final = df1.tail(10)['T1(Deg C)'].mean()  # Average of last 10 tank temperatures
        #convert to kelvin
        T1_initial_K = T1_initial + 273.15
        T1_final_K = T1_final + 273.15
        #adjust for ambient pressure
        ambient_pressure_inHg = 29.9
        ambient_pressure_psi = ambient_pressure_inHg * 0.491154  # 1 inHg = 0.491154 PSI
        P1_abs_initial = P1_gauge_initial + ambient_pressure_psi
        P1_abs_final = P1_gauge_final + ambient_pressure_psi
        
        
        # Mass flow rate is in g/min, time is in seconds
        time_steps = df1['Time(s)'].to_numpy()
        mass_flow_rates = df1['Mass Flowrate(g/min)'].to_numpy()
        time_minutes = time_steps / 60.0
        
        # Integrate mass flow rate over time using trapezoidal rule
        # Only consider positive flow rates (mass being added)
        positive_flow_mask = mass_flow_rates > 0
        if np.any(positive_flow_mask):
            positive_times = time_minutes[positive_flow_mask]
            positive_flows = mass_flow_rates[positive_flow_mask]
            m_added_grams = np.trapezoid(positive_flows, positive_times)  # grams
        else:
            m_added_grams = 0
            
        m_addedtotank = m_added_grams / 1000.0  # Convert to kg
        
        # Calculate mass left in tank using ideal gas law
        ratio = (P1_abs_final * T1_initial_K) / (P1_abs_initial * T1_final_K)
        
        # m_lefttank = m_addedtotank * (1+1/((P2*T1)/(P1*T2)-1))
        m_lefttank = m_addedtotank * (1 + 1/(ratio - 1))
        
        print(f"\nTrial {trial} Results:")
        print(f"  Initial conditions (avg of first 10): P1={P1_abs_initial:.3f} PSI, T1={T1_initial:.2f}°C ({T1_initial_K:.2f}K)")
        print(f"  Final conditions (avg of last 10): P1={P1_abs_final:.3f} PSI, T1={T1_final:.2f}°C ({T1_final_K:.2f}K)")
        print(f"  Mass added to tank: {m_addedtotank:.4f} kg ({m_added_grams:.2f} g)")
        print(f"  Mass left in tank: {m_lefttank:.4f} kg")
        print(f"  Mass that escaped: {m_addedtotank - m_lefttank:.4f} kg")
        print(f"  Pressure ratio (P2*T1)/(P1*T2): {ratio:.4f}")

        # Create integration plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # Plot mass flow rate
        ax1.plot(time_steps, mass_flow_rates, 'b-', linewidth=2)
        ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        ax1.set_ylabel('Mass Flow Rate (g/min)')
        ax1.set_title(f'Trial {trial}: Mass Flow Analysis')
        ax1.grid(True, alpha=0.3)
        # Calculate and plot cumulative mass
        cumulative_mass = np.cumsum(np.where(mass_flow_rates > 0, mass_flow_rates * np.gradient(time_minutes), 0.0))
        ax2.plot(time_steps, cumulative_mass, 'r-', linewidth=2)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Cumulative Mass Added (g)')
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'trial_{trial}_mass_flow.png', dpi=300, bbox_inches='tight')
        plt.show()

        #For Part 2
        filepath2 = f"rawdata/Lab 2 - Part 2{trial}.txt"
        df2 = pd.read_csv(filepath2, sep="\t", skiprows=2, header=0)
        
        # Convert Heater Energy from kJ to J for power calculation in Watts
        df2['Heater Energy (J)'] = df2['Heater Energy (kJ)'] * 1000
        # We assume the last 20% of the data represents the steady-state period
        steady_state_df = df2.tail(int(len(df2) * 0.2))
        # Use linear regression (polyfit) on the steady-state data for accuracy.
        # The slope of the Energy vs. Time graph is Power (dQ/dt).
        time_steady = steady_state_df['Time(s)'].to_numpy()
        energy_steady = steady_state_df['Heater Energy (J)'].to_numpy()
        #get the slope of energy vs time
        slope, intercept, r_value, p_value, std_err = linregress(time_steady, energy_steady)
        heat_loss_power = slope
        avg_steady_temp = steady_state_df['T1(Deg C)'].mean()

        print(f"\n--- Part 2 Analysis for Trial {trial} ---")
        print(f"Steady-state identified from t={time_steady[0]:.1f}s to t={time_steady[-1]:.1f}s.")
        print(f"Average temperature during steady state: {avg_steady_temp:.2f} °C")
        print(f"Calculated Heat Loss Rate (Input Power): {heat_loss_power:.2f} Watts (Std Err: {std_err:.3f})")

        # --- Heat Loss Breakdown with Uncertainty ---
        # dQ/dt = 2 * k * pi * L * (delta_T / ln(r2/r1))
        
        # Define physical constants with uncertainties.
        k_acrylic = ufloat(0.185, 0.015)      # Thermal conductivity of acrylic (W/m·K)
        L_tank = ufloat(11.25/39.3701, 0.001*/39.3701)      # Height of the tank in meters (converted from inches)
        r_outer = ufloat(8/39.3701, 0.045/39.3701)   # Inner radius in meters (converted from inches)
        thickness = ufloat((3/8)/39.3701, 0.15*(3/8)/39.3701) #thickness converted from inches to m
        #calculate r_inner
        r_inner = r_outer - thickness
        
        #Calculate the constant volume specific heat capacity of air by looking at the rate of temperature rise.   



        #For Part 3

        #Calculate the work done by the propeller using the above mentioned method. What would be the rate of rise in temperature due to work done by the propeller?


        #How significant is the work done on the system versus the heat added?
        
if __name__ == "__main__":
    main()