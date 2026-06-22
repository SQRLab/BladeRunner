"""
Power sweep for N=11 and N=12, then plot all N=2..12 together.

Run this AFTER the notebook has been executed through the "save" cell,
which writes power_sweep_N2_10.pkl in the same directory.

If power_sweep_N11_12.pkl already exists the heavy computation is skipped
and the script goes straight to plotting.
"""

import os
import time
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants

import importlib
import tweezer_functions
importlib.reload(tweezer_functions)
from tweezer_functions import (
    run_optimal_mode_selection_tweezed_only_claude_large,
    run_optimal_mode_selection_untweezed_only_claude_large,
)
from IonChainTools import *

# ---------------------------------------------------------------------------
# Physics setup (must match the notebook)
# ---------------------------------------------------------------------------
pi = np.pi
m = 39.9626 * constants.atomic_mass
c = constants.c
hbar = constants.hbar

tweezer_wavelength = 532e-9
omega_tweezer = 2 * pi * c / tweezer_wavelength

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "S_P_only.csv"), sep=",", encoding="UTF-8")
lambdares = np.array(df["wavelength (nm)"]) * 1e-9
omega_res = 2 * pi * c / lambdares
linewidths = np.array(df["A_ki (s^-1)"])

f_rf_r = 3e6
f_rf_a = 125e3
w0 = 1e-6

ueq = {N: ion_spacing(N, pi * f_rf_a)[0] for N in range(2, 40)}

# ---------------------------------------------------------------------------
# Load N=2..10 results from notebook pickle
# ---------------------------------------------------------------------------
nb_pkl = os.path.join(HERE, "power_sweep_N2_10.pkl")
print(f"Loading {nb_pkl} ...")
with open(nb_pkl, "rb") as f:
    nb = pickle.load(f)

power_results_nb     = nb["power_results"]           # list[N-2] of list[P] of winners
P_list               = nb["P_list"]
N_values_nb          = nb["N_values"]                # np.arange(2, 11)
optimal_powers_W_nb  = nb["optimal_powers_W"]
untweezed_results         = nb["untweezed_results"]
results_high_N_untweezed  = nb["results_high_N_untweezed"]
results_large_N_untweezed = nb["results_large_N_untweezed"]
results_test_check        = nb["results_test_check"]
results_high_N_tweezed    = nb["results_high_N_tweezed"]
results_large_N_tweezed   = nb["results_large_N_tweezed"]

# ---------------------------------------------------------------------------
# Power sweep for N=11, 12  (skip if already done)
# ---------------------------------------------------------------------------
large_pkl = os.path.join(HERE, "power_sweep_N11_12.pkl")
N_large = [11, 12]

if os.path.exists(large_pkl):
    print(f"Found {large_pkl}, loading instead of re-running.")
    with open(large_pkl, "rb") as f:
        large_data = pickle.load(f)
    power_results_large = large_data["power_results_large"]
else:
    power_results_large = {}
    for N in N_large:
        print(f"\nRunning power sweep for N={N} ({len(P_list)} powers)...")
        results_for_N = []
        t_start = time.perf_counter()
        for i, P in enumerate(P_list):
            if i % 10 == 0:
                elapsed = time.perf_counter() - t_start
                print(f"  [{i+1:3d}/{len(P_list)}]  P={P*1e3:.2f} mW  elapsed={elapsed:.1f}s")
            winners = run_optimal_mode_selection_tweezed_only_claude_large(
                omega_tweezer, linewidths, omega_res, m, mode_calc_r,
                N, f_rf_r, ueq, P, w0, max_tweezed=1,
            )
            results_for_N.append(winners)
        power_results_large[N] = results_for_N
        print(f"  N={N} done in {time.perf_counter()-t_start:.1f}s")

    with open(large_pkl, "wb") as f:
        pickle.dump({"power_results_large": power_results_large, "P_list": P_list}, f)
    print(f"\nSaved {large_pkl}")

# ---------------------------------------------------------------------------
# Build combined lookup structures
# ---------------------------------------------------------------------------

# Merge power results: {N: list-of-winners-per-P}
all_power_results = {}
for N in N_values_nb:
    all_power_results[N] = power_results_nb[N - 2]
for N in N_large:
    all_power_results[N] = power_results_large[N]

# Untweezed reference scores (power-independent)
untweezed_ref_dict = {}
for N_val, winners in untweezed_results:
    if not winners:
        continue
    win = winners[0]
    if isinstance(win, (list, tuple)) and isinstance(win[0], (list, tuple)) and len(win[0]) >= 3:
        untweezed_ref_dict[N_val] = abs(win[0][2])
    elif isinstance(win, (list, tuple)) and len(win) >= 3:
        untweezed_ref_dict[N_val] = abs(win[2])
for N_val, winners in results_high_N_untweezed:
    if winners:
        untweezed_ref_dict[N_val] = abs(winners[0][2])
for N_val, winners in results_large_N_untweezed:
    if winners:
        untweezed_ref_dict[N_val] = abs(winners[0][2])

# Tweezed fixed-P scores (for tau plot)
tweezed_fixed_data = {}
for N_val, winners in results_test_check:
    if not winners:
        continue
    win = winners[0]
    tweezed_ion = win[0][0] if isinstance(win[0], (list, tuple)) else win[0]
    if tweezed_ion is None:
        continue
    score = win[0][2] if isinstance(win[0], (list, tuple)) and len(win[0]) >= 3 else (win[2] if len(win) >= 3 else None)
    if score is not None:
        tweezed_fixed_data[N_val] = abs(score)
for N_val, winners in results_high_N_tweezed:
    if winners:
        tweezed_fixed_data[N_val] = abs(winners[0][2])
for N_val, winners in results_large_N_tweezed:
    if winners:
        tweezed_fixed_data[N_val] = abs(winners[0][2])

# ---------------------------------------------------------------------------
# Extract optimal powers for all N=2..12
# ---------------------------------------------------------------------------
N_values_all = np.arange(2, 13)

optimal_powers_W  = []
optimal_powers_mW = []
optimal_ratios_all = []

fig, axes = plt.subplots(4, 3, figsize=(15, 20))
axes = axes.flatten()

for idx, N in enumerate(N_values_all):
    x_tw, y_tw = [], []
    for P_val, winners in zip(P_list, all_power_results[N]):
        if not winners:
            continue
        win = winners[0]
        if len(win) >= 3 and win[0] is not None:
            x_tw.append(P_val)
            y_tw.append(abs(win[2]))

    x_tw = np.array(x_tw)
    y_tw = np.array(y_tw)

    untweezed_ref = untweezed_ref_dict.get(N)
    if untweezed_ref is None or len(y_tw) == 0:
        axes[idx].set_visible(False)
        continue

    ratio = untweezed_ref / y_tw
    opt_idx = np.argmin(ratio)
    opt_P_W  = x_tw[opt_idx]
    opt_P_mW = opt_P_W * 1e3
    opt_ratio = ratio[opt_idx]

    optimal_powers_W.append(opt_P_W)
    optimal_powers_mW.append(opt_P_mW)
    optimal_ratios_all.append(opt_ratio)

    axes[idx].scatter(x_tw * 1e3, ratio, marker='o', color='blue', s=50)
    axes[idx].set_xlabel("Optical power P (mW)", fontsize=10)
    axes[idx].set_ylabel(r"$\tau$ Ratio", fontsize=10)
    axes[idx].set_title(f"N = {N}", fontsize=11, fontweight='bold')
    axes[idx].grid(True, alpha=0.3)
    axes[idx].scatter(opt_P_mW, opt_ratio, color='red', s=200, marker='*', zorder=5,
                      label=f'Opt: {opt_P_mW:.2f} mW')
    axes[idx].legend(fontsize=9)

for idx in range(len(N_values_all), len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.suptitle("Power Sweep: Tweezed/Untweezed Ratio (N=2..12)",
             fontsize=14, fontweight='bold', y=1.00)
plt.savefig(os.path.join(HERE, "power_sweep_all_N.pdf"), bbox_inches='tight')
plt.show()

optimal_powers_W  = np.array(optimal_powers_W)
optimal_powers_mW = np.array(optimal_powers_mW)

print("\nOptimal Powers Summary (N=2..12):")
print("=" * 60)
for N, P_W, P_mW, r in zip(N_values_all, optimal_powers_W, optimal_powers_mW, optimal_ratios_all):
    print(f"N = {N}: P_opt = {P_W:.6e} W = {P_mW:.4f} mW,  Ratio = {r:.4f}")
print("=" * 60)

# ---------------------------------------------------------------------------
# Fig: Tau vs N at fixed P (using pre-calculated single-P scores)
# ---------------------------------------------------------------------------
N_fixed = sorted(tweezed_fixed_data.keys())
N_fixed_common = [N for N in N_fixed if N in untweezed_ref_dict]
tau_tw_fixed = np.array([1.0 / tweezed_fixed_data[N] for N in N_fixed_common])
tau_un_fixed = np.array([1.0 / untweezed_ref_dict[N] for N in N_fixed_common])

fig, ax = plt.subplots(1, 1, figsize=(7, 5))
ax.scatter(N_fixed_common, tau_tw_fixed, marker='o', color='blue', label='Tweezed', zorder=3)
ax.scatter(N_fixed_common, tau_un_fixed, marker='s', color='red', label='Untweezed', zorder=3)
ax.set_xlabel('N (Chain Length)', fontsize=12)
ax.set_ylabel(r'$\tau$ (cooling time)', fontsize=12)
ax.set_title(r'Cooling Time $\tau = 1/|coupling|$ per N (fixed $P$)', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(HERE, "tau_vs_N_fixed_P.pdf"), bbox_inches='tight')
plt.show()

ratio_fixed = tau_tw_fixed / tau_un_fixed
print("\nN   |  tau_tweezed/tau_untweezed  (fixed P)")
print("-" * 45)
for N, r in zip(N_fixed_common, ratio_fixed):
    print(f"N={N}  |  {r:.4f}")

fig, ax = plt.subplots(1, 1, figsize=(7, 5))
ax.scatter(N_fixed_common, ratio_fixed, marker='o', color='green', zorder=3)
ax.set_xlabel('N', fontsize=12)
ax.set_ylabel(r'Normalized $1/\eta$', fontsize=12)
ax.set_ylim(0.55, 1)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(HERE, "ratio_vs_N_fixed_P.pdf"), bbox_inches='tight')
plt.show()

# ---------------------------------------------------------------------------
# Fig: Tau vs N at optimal power per N
# ---------------------------------------------------------------------------
tweezed_opt_data = {}
for i, N in enumerate(N_values_all):
    if i >= len(optimal_powers_W):
        continue
    opt_P = optimal_powers_W[i]
    P_idx = int(np.argmin(np.abs(P_list - opt_P)))
    winners = all_power_results[N][P_idx]
    if not winners:
        continue
    win = winners[0]
    if len(win) >= 3 and win[0] is not None:
        tweezed_opt_data[N] = abs(win[2])

N_opt = sorted(set(tweezed_opt_data.keys()) & set(untweezed_ref_dict.keys()))
tau_tw_opt = np.array([1.0 / tweezed_opt_data[N] for N in N_opt])
tau_un_opt = np.array([1.0 / untweezed_ref_dict[N] for N in N_opt])
ratio_opt  = tau_tw_opt / tau_un_opt

fig, ax = plt.subplots(1, 1, figsize=(7, 5))
ax.scatter(N_opt, tau_tw_opt, marker='o', color='blue', label='Tweezed (opt P)', zorder=3)
ax.scatter(N_opt, tau_un_opt, marker='s', color='red', label='Untweezed', zorder=3)
ax.set_xlabel('N (Chain Length)', fontsize=12)
ax.set_ylabel(r'$\tau$ (cooling time)', fontsize=12)
ax.set_title(r'Cooling Time $\tau = 1/|coupling|$ at Optimal Power', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(HERE, "tau_vs_N_opt_P.pdf"), bbox_inches='tight')
plt.show()

print("\nN   |  tau_tweezed/tau_untweezed  (at optimal power)")
print("-" * 52)
for N, r in zip(N_opt, ratio_opt):
    print(f"N={N}  |  {r:.4f}")

fig, ax = plt.subplots(1, 1, figsize=(7, 5))
ax.scatter(N_opt, ratio_opt, marker='o', color='green', zorder=3)
ax.set_xlabel('N', fontsize=12)
ax.set_ylabel(r'Normalized $1/\eta$', fontsize=12)
ax.set_ylim(0.55, 1)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(HERE, "ratio_vs_N_opt_P.pdf"), bbox_inches='tight')
plt.show()
