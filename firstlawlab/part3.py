#manufacturer data
D_1 = 2.5 * 0.0254 #converted from inches to cm
n_1 = 4200 #rpm
Q_1 = 25 #cfm
rho_1 = 1.2 #kg/m^3 standard conditions
P_1 = 0.73549875 #0.001 hp to w

#calculate our experimental data
# Q_2/Q_1 = n_2/n_1 * (D_2/D_2)^3
# Q_2 * n_1 = Q_1 * n_2 because there's only 1 impeller diamater
# assume shaft speed is the same as what the manufacturer said as we have no way of getting it, we dont have a sensor for circulating air
# rho_2 = final mass / volume = Final Pressure / (R * Final Temperature), from ideal gas assumption
rho_2 = 54.39550459999999 * 6.89475729 / (0.2870 * (40+273)) # for trial a, converting from psi to kpa, using r from textbook,using kelvin
P_2 = P_1 * rho_2/rho_1

print(f"Work done by propeller: {P_2}")

#calculate rate of rise in temperature due to work done by propeller
m = 0.051254666666666664#
c_v = 63000 #from experiment
print(P_2/(m * c_v))



 