import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

# Parameters
repertoire = '/Users/joro/Documents/2025-2026/BA4/Physique\ numérique/Phys_num/Exercise2_student/rotatingpendulum/problème'
executable = '/engine'
input_filename = repertoire + '/configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'tf': 2, # t final (overwritten if N >0)
    'N': 0, # number of excitation periods
    'nsteps': 10000, # number of time steps per period (if N>0), number of timesteps total if N=0
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
theta0 = input_parameters['theta0']
thetadot0 = input_parameters['thetadot0']

paramstr = 'nsteps' # The parameter to scan, must be one of the keys in input_parameters

outstr = f"pendulum_kappa_{input_parameters['kappa']:.2g}_r_{input_parameters['r']:.2g}_Omega_{input_parameters['Omega']:.2g}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

nsteps_array = 2**np.arange(4, 10) # Array of nsteps to scan
nsimul = len(nsteps_array)


# ---- exact characteristic time ----
t_ref = np.linspace(0, tf, 200000)

#Exact energy
def Emec_calc(theta, thetadot, t) :
    return 0.5*m*np.power(L*thetadot,2) + m*g*L*(1-np.cos(theta))

def theta_calc(theta0, t) :
    w0 = np.sqrt(g/L)
    return theta0*np.cos(w0*t)

theta_exact = theta_calc(theta0, t_ref)
Emec_exact = Emec_calc(theta0, 0, 0)

outputs = []
totalsteps = []
theta_list = []
thetadot_list = []
Emec_list = []
Pnc_list = []


for i in range(len(nsteps_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = nsteps_array[i]

    output_file = f"{outstr}_{paramstr}_{nsteps_array[i]:.2g}.txt"
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

    dt = tf / nsteps_array[i]
    axs.plot(t, theta, label=f"dt={dt:.2e}", linewidth=lw, alpha=0.7)


#ax1.axhline(y=Emec_exact, color='k', linestyle='--', linewidth=2, label="Exact")
#ax1.set_xlabel(r'$\overline{t}$', fontsize=fs)
#ax1.set_ylabel(r'$\overline{E}$', fontsize=fs)
#ax1.set_xlim(0, tf)
#axs.set_ylim(0, Nf*1.2)
#plt.legend(fontsize=10)
#plt.grid(True)
#plt.tight_layout()
#figstr = "Energie_vs_t"
#plt.savefig(os.path.join(outdir, f"{figstr}_time.png"), dpi=300)

axs.plot(t_ref, theta_exact, 'k--', linewidth=2, label="Exact")
axs.set_xlabel(r't', fontsize=fs)
axs.set_ylabel(r'$\Theta$', fontsize=fs)
axs.set_xlim(0, tf)
plt.legend(fontsize=10)
plt.grid(True)
plt.tight_layout()
figstr = "Theta_vs_t"
plt.savefig(os.path.join(outdir, f"{figstr}_time.png"), dpi=300)

plt.show()