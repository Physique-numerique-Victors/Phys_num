import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline # If you don't have this, you can use np.interp instead, but it may be less accurate
import os

# Parameters
repertoire =  '/Users/joro/Documents/2025-2026/BA4/Physique\ numérique/Phys_num/Exercise1_student/runaways/problème'
executable = '/engine2' # Change this if your executable has a different name or path, like last week
input_filename = repertoire + '/configuration.in.example'

tf = 32
N0 = 0.0
g = 0.5
d = 0.01

alpha = 0 # 1 explicit, 0 implicit, 0.5 semi-implicit

if alpha == 1:
    alphastr = "expl"
elif alpha == 0:
    alphastr = "impl"
else:
    alphastr = "semi_impl"

figstr = f"runaway_{alphastr}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Outputs_g_{g:.2g}_d_{d:.2g}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)
# -------------------------------------------------
dt = 1 #TODO: Adjust for your needs
#nsimul = len(dt)
# Exact solution #TODO: Fill
beta = np.sqrt(g*g+4*d)
Nfp = (g+beta)/2 # steady state solution at t=inf

tol = np.logspace(-1, -12, 15)

nsimul = len(tol)

def N_anal(t):
    exp_term = np.exp(-beta*t)
    return (beta+g)*(1-exp_term)/(2*(1+exp_term*(beta+g)/(beta-g)))
Nf = N_anal(tf)  # exact solution at tf

Nr = 0.2  # fraction of equilibrium defining characteristic time

# ---- exact characteristic time ----
t_ref = np.linspace(0, tf, 200000)

#TODO: calculate N_exact as function of time
N_exact = N_anal(t_ref)# exact solution as function of time

ratio_exact = N_exact / Nfp
#TODO: calculate tau_ref as the time when ratio_exact crosses Nr, using interpolation
tau_ref = np.interp(Nr, N_exact/Nfp, t_ref)
paramstr = 'tol'
param = tol

# Simulations
outputs = []
totalsteps = []
tau_list = []
N_list = []
error = np.zeros(nsimul)

for i in range(nsimul):
    dt_val = param[i]  # current dt

    output_file = f"{alphastr}_tol={dt_val:.15g}.out"
    output_path = os.path.join(outdir, output_file)
    outputs.append(output_path)
    # Almost all parameters are passed as command line arguments, but you can also use an input file if you prefer. Adjust the command below accordingly.
    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{paramstr}={dt_val:.15g} output={output_path}"
        f" alpha={alpha:.2g} tf={tf:.3f} N0={N0:.3f} g={g:.4f} d={d:.4f}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


error = np.zeros(nsimul)


lw = 1.5
fs = 16

fig, axs = plt.subplots(1, 1)


for i in range(nsimul):
    with open(outputs[i]) as f:
        lines = f.readlines()

        total_steps = int(lines[-1].split(":")[1])
        data = np.loadtxt(lines[:-1])
        t = data[:, 0]
        N = data[:, 1]

        NN = N[-1]
        N_list.append(NN)
        totalsteps.append(total_steps)

        #TODO: calculate ratio and tau using interpolation, and store in tau_list
        ratio = N/Nfp # ratio as function of time

        if ratio[0] <= Nr <= ratio[-1]: # Check if Nr is within the range of ratio for interpolation
            try:
                tau = np.interp(Nr, ratio, t) #TODO: interpolate to find tau where ratio crosses Nr
            except ValueError:
                tau = np.nan
        else:
            tau = np.nan

        tau_list.append(tau)

        error[i] = np.abs(NN-Nf)/Nf #TODO: calculate relative error on Nf and store in error[i]

    axs.plot(t, N, label=f"dt={param[i]:.2e}", linewidth=lw, alpha=0.7)


plt.plot(t_ref, N_exact, 'k--', linewidth=2, label="Exact")
axs.set_xlabel(r'tolerance', fontsize=fs)
axs.set_ylabel(r'$\overline{N}$', fontsize=fs)
axs.set_xlim(0, tf)
axs.set_ylim(0, Nf*1.2)
plt.legend(fontsize=10)
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_time.png"), dpi=300)

# Error vs dt
tollist = tol

plt.figure()
plt.loglog(tollist, error, 'r+-', label="numerical")
plt.loglog(tollist, tollist/1e6, 'k--', label="O(dt)")
plt.loglog(tollist, tollist**2/1e6, 'k-.', label="O(dt^2)")
plt.xlabel(r"tolerance")
plt.ylabel("Relative error on Nf")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_Nf_error.png"), dpi=300)

# Convergence plot
plt.figure()
plt.plot(tollist, N_list, 'r+-', label="numerical")
plt.axhline(Nf, color='k', linestyle='--', label="Exact")
plt.xlabel(r"tolerance")
plt.ylabel(r"Final $\overline{N}$")
plt.xscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_Nf_conv.png"), dpi=300)

tau_ref = 2.9943409481678906

plt.figure()
plt.plot(tollist, tau_list, 'r+-', label="numerical")
plt.axhline(tau_ref, color='k', linestyle='--', label="Exact")
plt.xlabel(r"tolerance")
plt.ylabel(r"Characteristic time $\overline{\tau}$")
plt.xscale('log')
#plt.ylim(0, tf/10)  # Set y-limits to focus on the relevant range
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau.png"), dpi=300)

tau_err = np.abs(1 - np.array(tau_list) / tau_ref)

plt.figure()
plt.loglog(tollist, tau_err, 'r+-', label="numerical")
plt.loglog(tollist, tollist, 'k--', label="O(dt)")
plt.loglog(tollist, tollist**2, 'k-.', label="O(dt^2)")
plt.xlabel(r"tolerance")
plt.ylabel(r"Relative error on $\overline{\tau}$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau_error.png"), dpi=300)

plt.figure()

plt.loglog(totalsteps, tau_err, 'r+-', label=f"{alphastr}")
plt.xlabel("Total steps")
plt.ylabel("Relative error on tau")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau_error_vs_steps.png"), dpi=300)

plt.show()