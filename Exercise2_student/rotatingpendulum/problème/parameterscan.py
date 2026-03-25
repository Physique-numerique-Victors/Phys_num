import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.special import ellipk

# Parameters
repertoire = '/Users/victorgrosjean/Documents/Cours/ba4/P_num/Phys_num/Exercise2_student/rotatingpendulum/problème'
executable = '/engine'
input_filename = repertoire + '/configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'tf': 10, # t final (overwritten if N >0)
    'N': 0, # number of excitation periods
    'nsteps': 2**8, # number of time steps per period (if N>0), number of timesteps total if N=0
    'r': 0.0,
    'kappa': 0.0,
    'm': 0.1,
    'L': 0.2,
    'g': 9.81,
    'Omega': 2,
    'theta0': 1E-8,
    'thetadot0': 0.,
    'sampling': 1
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

# Extraction des paramètres pour les calculs analytiques
tf = input_parameters['tf']
m = input_parameters['m']
L = input_parameters['L']
g = input_parameters['g']
thetadot0 = input_parameters['thetadot0']
nsteps = input_parameters['nsteps']

paramstr = 'theta0'

outstr = f"pendulum_kappa_{input_parameters['kappa']:.2g}_r_{input_parameters['r']:.2g}_Omega_{input_parameters['Omega']:.2g}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

theta0_array = np.linspace(1e-8, np.pi-1e-3, 100) # Array of theta0 to scan
nsimul = len(theta0_array)


#Exact period
T_exact = 4*np.sqrt(L/g)*ellipk(np.sin(theta0_array*0.5)**2)


#Find the zeros of a list
def find_zero(list, t):
    times = []
    zeros = []
    for i in range(len(list)-1):
        if list[i]*list[i+1] < 0 and list[i] > 0:
            zeros.append(i)
    for j in zeros:
        a = (list[j+1]-list[j])/(t[j+1]-t[j])
        b = (list[j]*t[j+1]-list[j+1]*t[j])/(t[j+1]-t[j])
        times.append(-b/a)

    return times 


outputs = []
totalsteps = []
theta_list = []
thetadot_list = []
Emec_list = []
Pnc_list = []
t_list = []
T_list = []



for i in range(nsimul):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = theta0_array[i]

    output_file = f"{outstr}_{paramstr}_{theta0_array[i]:.8f}.txt"
    output_path = os.path.join(outdir, output_file)
    outputs.append(output_path)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")

error = np.zeros(nsimul)

lw = 1.5
fs = 16

fig, axs = plt.subplots(1, 1)

for i in range(nsimul):
    data = np.loadtxt(outputs[i])
    t = data[:, 0]
    theta = data[:, 1]
    thetadot = data[:, 2]
    Emec = data[:, 3]
    Pnc = data[:, 4]
    theta_list.append(theta)
    thetadot_list.append(thetadot)
    Emec_list.append(Emec)
    Pnc_list.append(Pnc)
    t_list.append(t)

    zeros = find_zero(theta, t)

    if len(zeros) >= 2:
        T_list.append(abs(zeros[1] - zeros[0]))
    else:
        T_list.append(np.nan)
    ############################################

########################################################################
T_list = np.array(T_list)
T_err = np.abs(T_exact - T_list) / T_exact

plt.plot(theta0_array, T_list, label='Période numérique', linewidth = lw, color = 'red')
plt.plot(theta0_array, T_exact, 'k--', linewidth = 2, label='Période exacte')
plt.xlabel(r'$\theta_0$', fontsize=fs)
plt.ylabel(r'$T$', fontsize=fs)
plt.xlim(0, np.pi)
plt.legend(fontsize = 10, loc = 'upper left')
plt.grid(True)
plt.tight_layout()
figstr = f"Period_vs_theta0.png"
plt.savefig(os.path.join(outdir, figstr), dpi=300)
plt.figure()

plt.plot(theta0_array, T_err, label='Erreur relative', linewidth = lw, color = 'red')
plt.xlabel(r'$\theta_0$', fontsize=fs)
plt.ylabel(r'Erreur relative', fontsize=fs) 
plt.xlim(0, np.pi)
plt.yscale('log')
plt.legend(fontsize = 10, loc = 'upper left')
plt.grid(True)
plt.tight_layout()
figstr = f"Period_error_vs_theta0.png"
plt.savefig(os.path.join(outdir, figstr), dpi=300)

plt.show()