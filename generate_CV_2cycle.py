import numpy as np
import matplotlib.pyplot as plt

# --- 1. Fixed Electrochemical Parameters (The System) ---
n = 1           # Number of electrons transferred
D = 1e-5        # Diffusion coefficient (cm^2/s)
C_bulk = 1.0E-7    # Bulk concentration (mol/cm^3)
C_red_init = 0  # Initial concentration of reduced species

# --- 2. Fixed Experimental Parameters (The Experiment) ---
E_i = 0.3       # Initial potential (V)
E_lambda = -0.3 # Switching potential (V)
E_f = 0.3       # Final potential (V)
v = 0.1         # Scan rate (V/s)
E_start = 0     # Formal potential E^0 (V)

# --- 3. Physical Constants ---
R_const = 8.314  # Gas constant
T = 298.15       # Temperature (K)
F = 96485        # Faraday constant
sigma = F / (R_const * T)

# --- 4. Numerical Parameters & TIME SETUP (MODIFIED) ---
L = 1000        # Spatial grid points
n_cycles = 2    # <--- CHANGE: Run for 2 cycles to stabilize the tail
# We will calculate M (time steps) dynamically below

# Calculate dimensions for ONE cycle first
voltage_range_1_cycle = abs(E_lambda - E_i) + abs(E_f - E_lambda)
T_1_cycle = voltage_range_1_cycle / v 

delta_x = 0.2 / L 
# We pick a delta_t that satisfies stability
delta_t = (0.45 * delta_x**2) / D  # Auto-calculate dt to be safe (alpha = 0.45)

# Calculate number of steps per cycle
steps_per_cycle = int(T_1_cycle / delta_t)

# Total M is steps per cycle * number of cycles
M = steps_per_cycle * n_cycles 

alpha = D * delta_t / (delta_x**2)

# --- 5. Potential Profile Calculation (MODIFIED) ---
# Create the potential array for ONE cycle
# 1. Forward Scan
steps_forward = int(steps_per_cycle / 2)
E_forward = np.linspace(E_i, E_lambda, steps_forward)

# 2. Reverse Scan
steps_backward = steps_per_cycle - steps_forward
E_backward = np.linspace(E_lambda, E_f, steps_backward)

# Combine into one cycle
E_cycle = np.concatenate((E_forward, E_backward))

# Repeat for n_cycles
E_vector = np.tile(E_cycle, n_cycles)

# Recalculate exact M in case of rounding differences
M = len(E_vector)
time_vector = np.linspace(0, M * delta_t, M)

# --- 6. Concentration Grid Initialization ---
C_O = np.ones((M, L)) * C_bulk  
C_R = np.ones((M, L)) * C_red_init 

# Initial conditions (t=0)
C_O[0, :] = C_bulk
C_R[0, :] = C_red_init
C_O[:, L - 1] = C_bulk
C_R[:, L - 1] = C_red_init

# --- 7. Main Iteration Loop ---
I_vector = np.zeros(M) 

for k in range(M - 1): 
    # --- A. Boundary Condition (Nernst) ---
    Lambda = np.exp(n * sigma * (E_vector[k] - E_start))
    C_R[k, 0] = C_bulk / (1 + Lambda)
    C_O[k, 0] = C_bulk - C_R[k, 0]
    
    # --- B. Calculate Current ---
    flux_O = (C_O[k, 1] - C_O[k, 0]) / delta_x
    I_vector[k] = -n * F * D * flux_O 

    # --- C. FDM Iteration ---
    # Update interior points
    C_O[k+1, 1:L-1] = C_O[k, 1:L-1] + alpha * (C_O[k, 2:L] - 2 * C_O[k, 1:L-1] + C_O[k, 0:L-2])
    C_R[k+1, 1:L-1] = C_R[k, 1:L-1] + alpha * (C_R[k, 2:L] - 2 * C_R[k, 1:L-1] + C_R[k, 0:L-2])

# --- 8. Final Current Calculation ---
Lambda_final = np.exp(n * sigma * (E_vector[M-1] - E_start))
C_R[M-1, 0] = C_bulk / (1 + Lambda_final)
C_O[M-1, 0] = C_bulk - C_R[M-1, 0]
flux_O_final = (C_O[M-1, 1] - C_O[M-1, 0]) / delta_x
I_vector[M-1] = -n * F * D * flux_O_final

# --- 9. Plotting (MODIFIED) ---
# We slice the data to plot ONLY the second cycle (the stable one)
start_idx = steps_per_cycle # Start plotting from end of cycle 1

plt.figure(figsize=(14, 10))

# Plot only from start_idx to end
plt.plot(E_vector[start_idx:], I_vector[start_idx:] * 1000000, color='blue', linewidth=3, label='Cycle 2')

# Optional: Plot Cycle 1 in dashed grey to see the difference
plt.plot(E_vector[:start_idx], I_vector[:start_idx] * 1000000, color='grey', linestyle='--', alpha=0.7, label='Cycle 1', linewidth=2)

plt.xlabel("Potential (V)", fontsize=14)
plt.ylabel("Current Density (μA/cm$^2$)", fontsize=14)
plt.title(f"Simulated Reversible CV (v = {v} V/s)", fontsize=16)
plt.grid(True, linestyle='--')
plt.gca().invert_xaxis() 
plt.legend(fontsize=14)
plt.show()