import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq

# Parameters
repertoire = '/Users/joro/Documents/2025-2026/BA4/Physique\ numérique/Phys_num/Exercise3_student/probleme'
executable = '/engine'
input_filename = repertoire + '/configuration.in.example'

input_parameters = {
    'rho0': 0., 
    'mA': 1,
    'mL': 7.3477* 1e22,
    'mT': 5.972 * 1e24,
    'G': 6.674e-11,
    'epsilon': 1e-3,
    'RT' : 6378.1 * 1e3,
    'd' : 384748* 1e3,
    'dt': 100,
    'tf': 365*86400, 
    'adaptatif': True,
    'xG' : 0, 'yG' : 0, 'vxG' : 1e5, 'vyG' : -4e3,
    'sampling': 1
}

# Paramêtres
epsilon = input_parameters['epsilon']
G = input_parameters['G']
mT = input_parameters['mT']
mA = input_parameters['mA']
mL = input_parameters['mL']
tf = input_parameters['tf']
RT = input_parameters['RT']
d = input_parameters['d']

# Positions et vitesses initiales du centre de masse G
xG = input_parameters['xG']
yG = input_parameters['yG']
vxG = input_parameters['vxG']
vyG = input_parameters['vyG']

# Positions et vitesses intiales de Terre, Lune et P0 (dans le ref absolu)
xT = -d*mL/(mT+mL)+xG; yT = yG
xL = d*mT/(mT+mL)+xG; yL = yG

xA = (xL+xT)/2 + xG; yA = d*np.sqrt(3)/2 + yG
xA_inG = xA - xG; yA_inG = yA - yG

omega = np.sqrt((mT+mL)*G/d**3)
r = np.sqrt(xA_inG**2 + yA_inG**2)
theta = np.arctan(yA_inG/xA_inG)

vxT = vxG; vyT = -mL*np.sqrt(G/(d*(mT+mL))) + vyG
vxL = vxG;  vyL = mT*np.sqrt(G/(d*(mT+mL))) + vyG
vxA = -omega*r*np.sin(theta) + vxG; vyA = omega*r*np.cos(theta) + vyG

input_parameters['xT'] = xT
input_parameters['yT'] = yT
input_parameters['vxT'] = vxT
input_parameters['vyT'] = vyT

input_parameters['xL'] = xL
input_parameters['yL'] = yL
input_parameters['vxL'] = vxL
input_parameters['vyL'] = vyL

input_parameters['xA'] = xA
input_parameters['yA'] = yA
input_parameters['vxA'] = vxA
input_parameters['vyA'] = vyA

# Energie mécanique totale


outstr = f"systeme_epsilon{input_parameters['epsilon']:.2g}"

outdir = f"Scan_{epsilon}_{mT}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

alpha_array = np.linspace(0, 2*np.pi, 5, endpoint = False)
nsimul = len(alpha_array)
outputs = []
dist = 1000000
for i in range(nsimul):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params['xA'] = xA + dist * np.cos(alpha_array[i])
    params['yA'] = yA + dist * np.sin(alpha_array[i])

    output_file = f"{outstr}_{'alpha'}_{i}.txt"
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

datasets = [np.loadtxt(i) for i in outputs]

# Terre
xT_inG_serie = []; yT_inG_serie = []
vxT_inG_serie = []; vyT_inG_serie = []

# Lune
xL_inG_serie = []; yL_inG_serie = []
vxL_inG_serie = []; vyL_inG_serie = []

dT_serie = []; dL_serie = []
err_dT_serie = []; err_dL_serie = []

xA_inG_serie = []
yA_inG_serie = []

for data in datasets:
    t  = data[:, 0]
    t_days = t/86400

    xT_inG = data[:, 1]; yT_inG = data[:, 2]; vxT_inG = data[:, 3]; vyT_inG = data[:, 4] # Terre
    xL_inG = data[:, 5]; yL_inG = data[:, 6]; vxL_inG = data[:, 7]; vyL_inG = data[:, 8] # Lune
    xA_inG = data[:, 9]; yA_inG = data[:, 10]; vxA_inG = data[:, 11]; vyA_inG = data[:, 12] # Artemis

    xT_inG_serie.append(xT_inG); yT_inG_serie.append(yT_inG); vxT_inG_serie.append(vxT_inG); vyT_inG_serie.append(vyT_inG)
    xL_inG_serie.append(xL_inG); yL_inG_serie.append(yL_inG); vxL_inG_serie.append(vxL_inG); vyL_inG_serie.append(vyL_inG)
    xA_inG_serie.append(xA_inG); yA_inG_serie.append(yA_inG)

    dt = data[:, 13]

    dT = np.sqrt((xT_inG-xA_inG)**2 + (yT_inG-yA_inG)**2)
    dT_serie.append(dT)
    err_dT_serie.append(np.abs((dT- d)/d))

    dL = np.sqrt((xL_inG-xA_inG)**2 + (yL_inG-yA_inG)**2)
    dL_serie.append(dL)
    err_dL_serie.append(np.abs((dL- d)/d))


lw = 1.5
fs = 16
ms = 8

#------------------------------------------------------------------
# Trajectoires dans le référentiel TOURNANT (synodique)
#------------------------------------------------------------------
plt.figure()

for i in range(nsimul):
    # On calcule l'angle de la Lune à chaque instant
    theta = np.arctan2(yL_inG_serie[i], xL_inG_serie[i])
    
    # On applique une rotation inverse pour rendre la Lune et la Terre fixes
    cos_t = np.cos(-theta)
    sin_t = np.sin(-theta)
    
    xA_rot = xA_inG_serie[i] * cos_t - yA_inG_serie[i] * sin_t
    yA_rot = xA_inG_serie[i] * sin_t + yA_inG_serie[i] * cos_t
    
    # On trace les orbites perturbées du troisième corps
    lbl_A = "Artemis (perturbée)" if i == 0 else None
    plt.plot(xA_rot, yA_rot, linewidth=1, label=f'alpha = {alpha_array[i]:.2f} rad')

plt.xlabel(r"$x'$ (dans R') [m]", fontsize=fs)
plt.ylabel(r"$y'$ (dans R') [m]", fontsize=fs)
plt.legend(fontsize=fs-5, loc='lower left')
plt.axis('equal')
plt.grid(True)
plt.show()