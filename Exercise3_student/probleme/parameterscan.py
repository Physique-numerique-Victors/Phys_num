import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq

# Parameters
repertoire = '/Users/victorgrosjean/Documents/Cours/ba4/P_num/Phys_num/Exercise3_student/probleme'
executable = '/engine'
input_filename = repertoire + '/configuration.in.example'

v0 = 1200
h_ref = 10000
mT = 5.972 * 1e24
G = 6.674e-11
RT = 6378.1e3
r0 = 314159 * 1e3
mu = G * mT
rp = RT + h_ref

v_ref = np.sqrt(v0**2 + 2*mu*(1/rp - 1/r0))
vt0 = rp/r0 * v_ref
vr0 = -np.sqrt(v0**2 - vt0**2)
theta_ref = np.arctan2(vt0, -vr0)

input_parameters = {
    'rho0': 0., 
    'mA': 8500,
    'mL': 7.3477 * 1e22,
    'mT': 5.972 * 1e24,
    'epsilon': 1e-3,
    'dt': 100,
    'tf': 1.728e5, # 2 jours
    'adaptatif': True,
    'xT': 0., 'yT': 0., 'xL': 384748 * 1e3, 'yL': 0, 'xA': r0, 'yA': 0,
    'vxT': 0., 'vyT': 0., 'vxL': 0., 'vyL': 0., 'vxA': vr0, 'vyA': vt0,
    'sampling': 1
}

epsilon = input_parameters['epsilon']
mT = input_parameters['mT']
tf = input_parameters['tf']


outstr = f"systeme_epsilon{input_parameters['epsilon']:.2g}"

outdir = f"Scan_{epsilon}_{mT}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)

def find_extremum(t, y, i_min, window=5, kind='min'):
    """
    Find a continuous extremum near discrete index i_min.
  
    Parameters
    ----------
    t, y     : full time and value arrays from RK4
    i_min    : index of the discrete extremum (from argmin/argmax)
    window   : number of points on each side to fit the spline on
    kind     : 'min' or 'max'
    """
    # Extract local window — enough points for a good spline, not too many
    lo = max(0, i_min - window)
    hi = min(len(t) - 1, i_min + window)
    t_loc = t[lo:hi+1]
    y_loc = y[lo:hi+1]

    # Fit spline (k=4 gives a differentiable cubic-like derivative, k=5 is smoother)
    # s=0 forces interpolation through all points (appropriate for RK4 output)
    spl = UnivariateSpline(t_loc, y_loc, k=5, s=0)
    dspl = spl.derivative()

    # Find root of derivative in the interval
    # brentq requires a sign change — check it exists
    da = dspl(t_loc[0])
    db = dspl(t_loc[-1])
    if da * db > 0:
        # No sign change found — fall back to denser window or return discrete min
        raise ValueError(f"No sign change in derivative over window [{t_loc[0]}, {t_loc[-1]}]")

    t_extremum = brentq(dspl, t_loc[0], t_loc[-1])
    y_extremum = spl(t_extremum)
    return float(t_extremum), float(y_extremum)

theta_array = np.linspace(0.189, 0.192, 50)
nsimul = len(theta_array)
outputs = []

for i in range(nsimul):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params['vxA'] = -v0 * np.cos(theta_array[i])
    params['vyA'] = v0 * np.sin(theta_array[i])

    output_file = f"{outstr}_{'epsilon'}_{theta_array[i]:.8f}.txt"
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

xA_serie = []
yA_serie = []
accA_serie = []
t_serie = []
h_serie = []
n_entries_list = []

accA_max = []
t_acc_max =[]
h_min = []
t_h_min = []
accA_nan =[]
t_acc_nan = []

h_atm = 100000 # Seuil ou l'on considère que Artémis traverse l'atmosphère terrestre

for data in datasets:
    t  = data[:, 0]
    xA = data[:, 1]
    yA = data[:, 2]
    vxA = data[:, 3]
    vyA = data[:, 4]
    dt = data[:, 5]
    accA = data[:, 6]
    Pt = data[:, 7]

    r = np.sqrt(xA**2 + yA**2)
    h = r - RT
    v = np.sqrt(vxA**2 + vyA**2)

    xA_serie.append(xA)
    yA_serie.append(yA)
    accA_serie.append(accA)
    t_serie.append(t)
    h_serie.append(h)

    i_accA_max = np.argmax(accA)
    i_h_min = np.argmin(h)

    try:
        t_a, a = find_extremum(t, accA, i_accA_max, window=5, kind='max')
        t_h, h_ = find_extremum(t, h, i_h_min, window=5, kind='min')
    except ValueError:
        t_a = t[i_accA_max]
        a = accA[i_accA_max]
        t_h = t[i_h_min]
        h_ = h[i_h_min]

    n_entries = 0

    for i in range(1,len(h)):
        if h[i] < h_atm and h[i-1] >= h_atm:
            n_entries += 1

    n_entries_list.append(n_entries)
    t_acc_max.append(t_a), accA_max.append(a)
    t_h_min.append(t_h), h_min.append(h_)

    if n_entries != 1:
        accA_nan.append(a)
        t_acc_nan.append(t_a)
        accA_max[-1] = np.nan # Invalide si Artémis traverse l'atmosphère terrestre
    else:
        accA_nan.append(np.nan)
        t_acc_nan.append(np.nan)


i_opt = np.nanargmin(accA_max)
theta_opt = theta_array[i_opt]
a_opt = accA_max[i_opt]

i_nan = np.argmin(accA_nan)
theta_nan = theta_array[i_nan]
a_nan = accA_nan[i_nan]

plt.figure()
plt.scatter(theta_array, accA_max, label='a_max', color='blue')
plt.scatter(theta_array, accA_nan, label='a_max non valide', color='red')
plt.xlabel('theta')
plt.ylabel('a_max')
plt.grid()
plt.legend()

plt.show()
    
print(f"Optimal theta: {theta_opt:.6f} rad, with a_max = {a_opt:.2f} m/s²")
print(f"First invalid theta: {theta_nan:.6f} rad, with a_max = {a_nan:.2f} m/s²")



