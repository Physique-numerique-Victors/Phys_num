import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"..."

plot_layout = {
    "theta_time": True,
    "phase_space": False,
    "energy": True,
    "real_space": False,
    "power": True,
    "energy_balance": True
}

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.txt")))

datasets = []
param_values = []
param_name = None

for f in files:

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".txt"

    parts = name.split("_")

    param_name = parts[-2]          # scanned parameter
    value = float(parts[-1])        # parameter value

    data = np.loadtxt(f)

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Helper: colored line
# ============================================================

def colored_line(x, y, t, ax, vmin=None, vmax=None):

    points = np.array([x, y]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='viridis')
    lc.set_array(t)

    if vmin is not None and vmax is not None:
        lc.set_clim(vmin, vmax)

    lc.set_linewidth(2)

    ax.add_collection(lc)
    ax.autoscale()

    return lc

# ============================================================
# Shared time range for color maps
# ============================================================

tmin = min(data[:,0].min() for data in datasets)
tmax = max(data[:,0].max() for data in datasets)

# ============================================================
# Axis layout helper
# ============================================================

def get_axes(plot_key, title):

    if plot_layout[plot_key]:

        fig, ax = plt.subplots()
        axes = [ax]*len(datasets)

    else:

        n = len(datasets)
        ncols = min(3, n)
        nrows = math.ceil(n/3)

        fig, axarr = plt.subplots(nrows, ncols,
                                  figsize=(5*ncols,4*nrows))

        axes = np.array(axarr).reshape(-1)

        for j in range(n, len(axes)):
            fig.delaxes(axes[j])

        axes = axes[:n]

    fig.suptitle(title)

    return fig, axes


cmap = plt.get_cmap("tab10")


# ============================================================
# Plot 1 : theta vs time
# ============================================================

fig, axes = get_axes("theta_time", "Angle vs time")

for i,data in enumerate(datasets):

    t = data[:,0]
    theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi

    color = cmap(i % 10)

    axes[i].plot(t, theta, color=color,
                 label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("theta")
    axes[i].grid()

    if not plot_layout["theta_time"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["theta_time"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"theta_vs_time_all.png"), dpi=300)


# ============================================================
# Plot 2 : phase space
# ============================================================

fig, axes = get_axes("phase_space", "Phase space")

lc_ref = None

for i,data in enumerate(datasets):

    t = data[:,0]
    theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi
    thetadot = data[:,2]

    lc = colored_line(theta, thetadot, t, axes[i], tmin, tmax)

    axes[i].set_xlabel("theta")
    axes[i].set_ylabel("thetadot")
    axes[i].grid()

    if not plot_layout["phase_space"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

    if lc_ref is None:
        lc_ref = lc

cbar = fig.colorbar(lc_ref, ax=axes)
cbar.set_label("time")

fig.savefig(os.path.join(fig_dir,"phase_space_all.png"), dpi=300)


# ============================================================
# Plot 3 : energy
# ============================================================

fig, axes = get_axes("energy", "Energy evolution")

for i,data in enumerate(datasets):

    t = data[:,0]
    E = data[:,3]

    color = cmap(i % 10)

    axes[i].plot(t, E, color=color,
                 label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Energy")
    axes[i].grid()

    if not plot_layout["energy"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["energy"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"energy_all.png"), dpi=300)


# ============================================================
# Plot 4 : real space trajectory
# ============================================================

L = 0.2

fig, axes = get_axes("real_space", "Real space trajectory")

lc_ref = None

for i,data in enumerate(datasets):

    t = data[:,0]
    theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi

    x = L*np.sin(theta)
    y = -L*np.cos(theta)

    lc = colored_line(x, y, t, axes[i], tmin, tmax)

    axes[i].set_xlabel("x")
    axes[i].set_ylabel("y")

    if not plot_layout["real_space"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

    if lc_ref is None:
        lc_ref = lc

cbar = fig.colorbar(lc_ref, ax=axes)
cbar.set_label("time")

fig.savefig(os.path.join(fig_dir,"real_space_all.png"), dpi=300)


# ============================================================
# Plot 5 : power
# ============================================================

fig, axes = get_axes("power", "Non-conservative power")

for i,data in enumerate(datasets):

    t = data[:,0]
    Pnc = data[:,4]

    color = cmap(i % 10)

    axes[i].plot(t, Pnc, color=color,
                 label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Pnc")
    axes[i].grid()

    if not plot_layout["power"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["power"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"power_all.png"), dpi=300)


# ============================================================
# Plot 6 : energy balance
# ============================================================

fig, axes = get_axes("energy_balance", "Energy balance")

for i,data in enumerate(datasets):

    color = cmap(i % 10)

    t = data[:,0]
    E = data[:,3]
    Pnc = data[:,4]

    dEdt = np.gradient(E,t)

    axes[i].plot(t, dEdt, linestyle='-', color=color, label="dE/dt")
    axes[i].plot(t, Pnc, linestyle='--', color=color, label="Pnc")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Power")
    axes[i].grid()

    if not plot_layout["energy_balance"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")
        axes[i].legend()

if plot_layout["energy_balance"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"energy_balance_all.png"), dpi=300)


plt.show()