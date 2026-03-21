import numpy as np
import subprocess
import os

# Parameters
repertoire = '/Users/victorgrosjean/Documents/Cours/ba4/P_num/Phys_num/Exercise2_student/rotatingpendulum/problème'
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
    'theta0': 0.1,
    'thetadot0': 0.,
    'sampling': 1
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'nsteps' # The parameter to scan, must be one of the keys in input_parameters

outstr = f"pendulum_kappa_{input_parameters['kappa']:.2g}_r_{input_parameters['r']:.2g}_Omega_{input_parameters['Omega']:.2g}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

nsteps_array = 2**np.arange(2, 10) # Array of nsteps to scan, not used in the code but can be useful to have




for i in range(len(nsteps_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = nsteps_array[i]

    output_file = f"{outstr}_{paramstr}_{nsteps_array[i]:.2g}.txt"
    output_path = os.path.join(outdir, output_file)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")

