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
eps = 1E-2

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

#Exact period
T_exact = 2*np.pi*np.sqrt(L/g)

def theta_calc(theta0, t) :
    w0 = np.sqrt(g/L)
    return theta0*np.cos(w0*t)

theta_exact = theta_calc(theta0, t_ref)
Emec_exact = Emec_calc(theta0, 0, 0)

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

    theta1 = find_zero(theta, t)[0]
    theta2 = find_zero(theta, t)[1]
    T_list.append(np.abs(theta2-theta1))
    ############################################

    dt = tf / nsteps_array[i]
    plt.plot(t, theta, label=f"dt={dt:.2e}", linewidth=lw, alpha=0.7)



plt.plot(t_ref, theta_exact, 'k--', linewidth=2, label="Exacte")
plt.xlabel(r't', fontsize=fs)
plt.ylabel(r'$\theta$', fontsize=fs)
plt.xlim(0, tf)
plt.legend(fontsize=10, loc = 'best')
plt.grid(True)
plt.tight_layout()
figstr = "Theta_vs_t"
plt.savefig(os.path.join(outdir, f"{figstr}.png"), dpi=300)

########################################################################
plt.figure()

for i in range(nsimul):
    dt = tf / nsteps_array[i]
    plt.plot(t_list[i], Emec_list[i], label=f"dt={dt:.2e}", linewidth=lw, alpha=0.7)

plt.axhline(y=Emec_exact, color='k', linestyle='--', linewidth=2, label="Exacte")
plt.xlabel(r'$t$', fontsize=fs)
plt.ylabel(r'$E_{m}$', fontsize=fs)
plt.xlim(0, tf)
plt.legend(fontsize=10, loc = 'best')
plt.grid(True)
plt.tight_layout()
figstr = "Energie_vs_t"
plt.savefig(os.path.join(outdir, f"{figstr}.png"), dpi=300)

########################################################################
dtlist = tf / nsteps_array
T_list = np.array(T_list)
T_err = np.abs(T_exact - T_list) / T_exact

plt.figure()
plt.loglog(dtlist, T_err, 'r+-', label="numérique")
#plt.loglog(dtlist, dtlist, 'k--', label="O(dt)")
#plt.loglog(dtlist, dtlist**2, 'k-.', label="O(dt^2)")
plt.xlabel(r"$dt$")
plt.ylabel(r"Erreur relative sur $T$")
plt.legend()
plt.grid(True)
plt.tight_layout()
figstr = "Error_on_energy_vs_dt"
plt.savefig(os.path.join(outdir, f"{figstr}.png"), dpi=300)

###################################################################

plt.figure()
plt.plot(dtlist, T_list, 'r+-', label="numérique")
plt.axhline(T_exact, color='k', linestyle='--', label="Exacte")
plt.xlabel(r"$dt$")
plt.ylabel(r"$T$")
plt.xscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
figstr = "Energy_vs_dt"
plt.savefig(os.path.join(outdir, f"{figstr}_tau.png"), dpi=300)

plt.show()