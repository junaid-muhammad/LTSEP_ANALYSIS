#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-03-13 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid III <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Import relevant packages
import uproot
import uproot as up
import numpy as np
import pandas as pd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
from array import array as arr
import csv
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString, TF1
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from ROOT import TGraphErrors
from functools import reduce
import math as ma
import uncertainties as u
from ctypes import c_double
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from copy import deepcopy
import random

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=4:
    print("!!!!! ERROR !!!!!\n Expected 4 arguments\n Usage is with - CSV file suffix \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
ITERATION = sys.argv[2]
XSECT_LEPS_DAT = sys.argv[3]
XSECT_HEPS_DAT = sys.argv[4]

################################################################################################################################################

'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__), "Plot_ProdCoin")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
XSECT_PATH = "%s/LTSEP_ANALYSIS/src/xsects" % (REPLAYPATH)
XSECT_OUTPATH = "%s/LTSEP_ANALYSIS/src/output" % (REPLAYPATH)
XSECT_PARAMPATH = "%s/LTSEP_ANALYSIS/src/fit_params" % (REPLAYPATH)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, lt_analysis path assumed as %s" % (USER, HOST, REPLAYPATH))
#Pion_Analysis_Distributions = "%s/%s_ProdCoin_Pion_Analysis_xsect_Distributions.pdf" % (XSECT_PATH, PHY_SETTING)

# Input file location and variables taking
xsect_loweps_dat = "%s/%s.dat" % (XSECT_PATH, XSECT_LEPS_DAT)
xsect_higheps_dat = "%s/%s.dat" % (XSECT_PATH, XSECT_HEPS_DAT)

###############################################################################################################################################

# Section for average yield calculation
print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("="*40)

# Define the output path for loweps cross-sections
xsect_loweps_output = "%s/%s_pion_physics_loweps_xsect_iter%s.csv" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)
# List your actual column names here, e.g.:
xsect_loweps_headers = ["Xsection_real", "Xsection_real_error", "Xsection_model", "epsilon", "theta_cm", "phi_central", "t_central", "W", "Q2"]
# Read the file with provided headers
xsect_loweps_df = pd.read_csv(xsect_loweps_dat, sep=r'\s+', names=xsect_loweps_headers, engine='python')
# Save as CSV
xsect_loweps_df.to_csv(xsect_loweps_output, index=False)
print(f"Arranged loweps cross-sections in CSV and saved to the path: {xsect_loweps_output}")


# Define the output path for higheps cross-sections
xsect_higheps_output = "%s/%s_pion_physics_higheps_xsect_iter%s.csv" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)
# List your actual column names here, e.g.:
xsect_higheps_headers = ["Xsection_real", "Xsection_real_error", "Xsection_model", "epsilon", "theta_cm", "phi_central", "t_central", "W", "Q2"]
# Read the file with provided headers
xsect_higheps_df = pd.read_csv(xsect_higheps_dat, sep=r'\s+', names=xsect_higheps_headers, engine='python')
# Save as CSV
xsect_higheps_df.to_csv(xsect_higheps_output, index=False)
print(f"Arranged higheps cross-sections in CSV and saved to the path: {xsect_higheps_output}")

print("-"*40)

#=========================================================================================================================================================================================================================================================

# Plotting Section
# Load dataframes
xsect_loweps_df = pd.read_csv(xsect_loweps_output)
xsect_higheps_df = pd.read_csv(xsect_higheps_output)

# Ensure error bars are positive
xsect_loweps_df["Xsection_real_error"] = xsect_loweps_df["Xsection_real_error"].abs()
xsect_higheps_df["Xsection_real_error"] = xsect_higheps_df["Xsection_real_error"].abs()

# Convert cross section and error from MeV to GeV (μb/MeV² → μb/GeV²)
xsect_loweps_df["Xsection_real"] = xsect_loweps_df["Xsection_real"] * 1e6
xsect_loweps_df["Xsection_real_error"] = xsect_loweps_df["Xsection_real_error"] * 1e6
xsect_higheps_df["Xsection_real"] = xsect_higheps_df["Xsection_real"] * 1e6
xsect_higheps_df["Xsection_real_error"] = xsect_higheps_df["Xsection_real_error"] * 1e6

# Use ALL unique t_central values (from either dataset), sorted
t_central_values = np.union1d(
    xsect_loweps_df["t_central"].unique(),
    xsect_higheps_df["t_central"].unique()
)
t_central_values = np.sort(t_central_values)

# Output file path
UnSep_Xsection_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_ProdCoin_Pion_Analysis_UnSep_xsection_iter%s_Distributions.pdf" % (REPLAYPATH, PHY_SETTING, ITERATION)

# 2x3 grid, first subplot is for annotation/text
fig_1, axs_1 = plt.subplots(2, 3, figsize=(20, 10))
axs_1 = axs_1.flatten()
axs_1[0].axis('off')
axs_1[0].text(
    0.5, 0.5, f"Physics Setting:\n{PHY_SETTING}\n Iteration: {ITERATION}",
    ha='center', va='center', fontsize=18, fontweight='bold',
    transform=axs_1[0].transAxes
)

for idx, t_central in enumerate(t_central_values[:5]):
    ax = axs_1[idx + 1]
    df_t_low = xsect_loweps_df[xsect_loweps_df["t_central"] == t_central]
    df_t_high = xsect_higheps_df[xsect_higheps_df["t_central"] == t_central]
    # Plot loweps (blue) if exists
    if not df_t_low.empty:
        ax.errorbar(
            df_t_low["phi_central"],
            df_t_low["Xsection_real"],
            yerr=df_t_low["Xsection_real_error"],
            fmt='o',
            color='blue',
            ecolor='blue',
            capsize=2,
            markersize=4,
            label='low ε'
        )
    # Plot higheps (red) if exists
    if not df_t_high.empty:
        ax.errorbar(
            df_t_high["phi_central"],
            df_t_high["Xsection_real"],
            yerr=df_t_high["Xsection_real_error"],
            fmt='s',
            color='red',
            ecolor='red',
            capsize=2,
            markersize=4,
            label='high ε'
        )
    ax.set_title(
        f"t central = {t_central:.3f}",
        fontsize=18, fontweight='bold'
    )
    ax.set_xlabel("phi (deg)", fontsize=16, fontweight='bold')
    ax.set_xlim(0, 360)
    ax.set_ylabel(r"$d^2\sigma/dt\,d\phi$ ($\mu$b/GeV$^2$)", fontsize=16, fontweight='bold')
    ax.set_ylim(0, max(df_t_low["Xsection_real"].max(), df_t_high["Xsection_real"].max()) * 1.1)  # Adjust y-limits for visibility
    ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig(UnSep_Xsection_pdf)
plt.close(fig_1)

print(f"Saved Unsep cross-section plots to {UnSep_Xsection_pdf}")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Full Rosenbluth function for fitting
def rosenbluth_fit(x, sigma_L, sigma_T, sigma_LT, sigma_TT):
    phi, eps = x
    phi_rad = np.deg2rad(phi)  # Convert phi to radians for np.cos
    A = eps * sigma_L + sigma_T
    B = np.sqrt(2 * eps * (eps + 1)) * sigma_LT
    C = eps * sigma_TT
    return ((1 / (2 * np.pi)) * (A + B * np.cos(phi_rad) + C * np.cos(2 * phi_rad)))

# Output file path
RS_Sep_Xsection_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_ProdCoin_Pion_Analysis_Rosenbluth_Sep_xsection_iter%s_Distributions.pdf" % (REPLAYPATH, PHY_SETTING, ITERATION)

# --- Prepare lists to collect fit results ---
sigma_L_list, sigma_L_err_list = [], []
sigma_T_list, sigma_T_err_list = [], []
sigma_LT_list, sigma_LT_err_list = [], []
sigma_TT_list, sigma_TT_err_list = [], []
abs_t_cent_list = []

# 2x3 grid, first subplot is for annotation/text
fig_2, axs_2 = plt.subplots(2, 3, figsize=(20, 10))
axs_2 = axs_2.flatten()
axs_2[0].axis('off')
axs_2[0].text(
    0.5, 0.5, f"Physics Setting:\n{PHY_SETTING}\n Iteration: {ITERATION}",
    ha='center', va='center', fontsize=18, fontweight='bold',
    transform=axs_2[0].transAxes
)

# --- Create 2x3 grid for main plots ---
for idx, t_central in enumerate(t_central_values[:5]):
    ax = axs_2[idx + 1]
    df_t_low = xsect_loweps_df[xsect_loweps_df["t_central"] == t_central]
    df_t_high = xsect_higheps_df[xsect_higheps_df["t_central"] == t_central]
    # Plot loweps (blue) if exists
    if not df_t_low.empty:
        ax.errorbar(
            df_t_low["phi_central"],
            df_t_low["Xsection_real"],
            yerr=df_t_low["Xsection_real_error"],
            fmt='o',
            color='blue',
            ecolor='blue',
            capsize=2,
            markersize=4,
            label='low ε'
        )
    # Plot higheps (red) if exists
    if not df_t_high.empty:
        ax.errorbar(
            df_t_high["phi_central"],
            df_t_high["Xsection_real"],
            yerr=df_t_high["Xsection_real_error"],
            fmt='s',
            color='red',
            ecolor='red',
            capsize=2,
            markersize=4,
            label='high ε'
        )
    # --- Fit Rosenbluth to both datasets ---
    # Combine data
    df_combined = pd.concat([df_t_low, df_t_high], ignore_index=True)
    if not df_combined.empty:
        phi_vals = df_combined["phi_central"].values
        eps_vals = df_combined["epsilon"].values
        x_vals = np.vstack((phi_vals, eps_vals))
        y_vals = df_combined["Xsection_real"].values 
        y_errs = df_combined["Xsection_real_error"].values

        # Initial guess for [sigma_L, sigma_T, sigma_LT, sigma_TT]
        p0 = [y_vals.mean()/2, y_vals.mean(), 0, 0]
        try:
            popt, pcov = curve_fit(
                rosenbluth_fit, x_vals, y_vals, sigma=y_errs, p0=p0, absolute_sigma=True, maxfev=50000
            )
            sigma_L, sigma_T, sigma_LT, sigma_TT = popt
            perr = np.sqrt(np.diag(pcov))
            # Store all fit results for plotting
            abs_t_cent_list.append(t_central)
            sigma_L_list.append(sigma_L)
            sigma_L_err_list.append(perr[0])
            sigma_T_list.append(sigma_T)
            sigma_T_err_list.append(perr[1])
            sigma_LT_list.append(sigma_LT)
            sigma_LT_err_list.append(perr[2])
            sigma_TT_list.append(sigma_TT)
            sigma_TT_err_list.append(perr[3])
            # Plot fit curve
            phi_fit = np.linspace(0, 360, 200)
            # Use average eps for each dataset to plot fit curves for low/high eps
            for eps_val, color, label in [
                (df_t_low["epsilon"].mean() if not df_t_low.empty else None, 'blue', 'Fit (low ε)'),
                (df_t_high["epsilon"].mean() if not df_t_high.empty else None, 'red', 'Fit (high ε)')
            ]:
                if eps_val is not None:
                    y_fit = rosenbluth_fit((phi_fit, np.full_like(phi_fit, eps_val)), sigma_L, sigma_T, sigma_LT, sigma_TT)
                    ax.plot(phi_fit, y_fit, color=color, linestyle='--', label=label)
            # Print cross-section fit results and residuals
            print(f"t_central = {t_central:.3f}")
            print(f"  sigma_L  = {sigma_L:.4f} ± {perr[0]:.4f}")
            print(f"  sigma_T  = {sigma_T:.4f} ± {perr[1]:.4f}")
            print(f"  sigma_LT = {sigma_LT:.4f} ± {perr[2]:.4f}")
            print(f"  sigma_TT = {sigma_TT:.4f} ± {perr[3]:.4f}")
            # Calculate chi^2/ndof for the fit
            y_fit_all = rosenbluth_fit((phi_vals, eps_vals), sigma_L, sigma_T, sigma_LT, sigma_TT)
            residuals = y_vals - y_fit_all
            chi2 = np.sum((residuals / y_errs) ** 2)
            ndof = len(y_vals) - len(popt)
            print(f"  chi2     = {chi2:.2f}")
            print(f"  ndof     = {ndof}")
            print(f"  chi2/ndof= {chi2/ndof if ndof > 0 else float('nan'):.2f}")
#            print("  Residuals (data - fit):")
#            for i, (phi, data, fit, err) in enumerate(zip(phi_vals, y_vals, y_fit_all, y_errs)):
#                print(f"    phi={phi:.2f}  data={data:.4f}  fit={fit:.4f}  err={err:.4f}  resid={data-fit:.4f}")
            # Annotate fit results
            fit_text = (
                r"$\sigma_L$ = %.3f $\pm$ %.3f" % (sigma_L, perr[0]) + "\n"
                r"$\sigma_T$ = %.3f $\pm$ %.3f" % (sigma_T, perr[1]) + "\n"
                r"$\sigma_{LT}$ = %.3f $\pm$ %.3f" % (sigma_LT, perr[2]) + "\n"
                r"$\sigma_{TT}$ = %.3f $\pm$ %.3f" % (sigma_TT, perr[3])
            )
#            fit_text += "\n$\chi^2$/ndof = %.2f/%d = %.2f" % (chi2, ndof, chi2/ndof if ndof > 0 else float('nan'))
            ax.text(
                0.98, 0.02, fit_text, ha='right', va='bottom', fontsize=14,
                transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray')
            )
        except Exception as e:
            ax.text(0.5, 0.5, "Fit failed", ha='center', va='center', color='red', fontsize=14, transform=ax.transAxes)
            print(f"Fit failed for t_central={t_central}: {e}")

    ax.set_title(
        f"t central = {t_central:.3f}",
        fontsize=18, fontweight='bold'
    )
    ax.set_xlabel("phi (deg)", fontsize=16, fontweight='bold')
    ax.set_xlim(0, 360)
    ax.set_ylabel(r"$d^2\sigma/dt\,d\phi$ ($\mu$b/GeV$^2$)", fontsize=16, fontweight='bold')
    ymax = max(df_t_low["Xsection_real"].max(), df_t_high["Xsection_real"].max())
    ax.set_ylim(-0.1 * ymax, 1.2 * ymax)
    ax.legend(fontsize=12)
plt.tight_layout()

# Save the figure with Rosenbluth fits
fig_2.savefig(RS_Sep_Xsection_pdf)
plt.close(fig_2)

print(f"Saved Rosenbluth Sep cross-section plots to {RS_Sep_Xsection_pdf}")

# --- Save separated cross-sections (sigma_L, sigma_T, sigma_LT, sigma_TT) vs t, theta_cm, Q2 to CSV ---
sep_xsection_csv = "%s/%s_pion_physics_sep_xsect_iter%s.csv" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)
sep_xsection_df = pd.DataFrame({
    't_central': abs_t_cent_list,
    'theta_cm': [xsect_loweps_df[xsect_loweps_df['t_central'] == t]['theta_cm'].iloc[0] if not xsect_loweps_df[xsect_loweps_df['t_central'] == t].empty else (xsect_higheps_df[xsect_higheps_df['t_central'] == t]['theta_cm'].iloc[0] if not xsect_higheps_df[xsect_higheps_df['t_central'] == t].empty else None) for t in abs_t_cent_list],
    'Q2': [xsect_loweps_df[xsect_loweps_df['t_central'] == t]['Q2'].iloc[0] if not xsect_loweps_df[xsect_loweps_df['t_central'] == t].empty else (xsect_higheps_df[xsect_higheps_df['t_central'] == t]['Q2'].iloc[0] if not xsect_higheps_df[xsect_higheps_df['t_central'] == t].empty else None) for t in abs_t_cent_list],
    'W': [xsect_loweps_df[xsect_loweps_df['t_central'] == t]['W'].iloc[0] if not xsect_loweps_df[xsect_loweps_df['t_central'] == t].empty else (xsect_higheps_df[xsect_higheps_df['t_central'] == t]['W'].iloc[0] if not xsect_higheps_df[xsect_higheps_df['t_central'] == t].empty else None) for t in abs_t_cent_list],
    'sigma_L': sigma_L_list,
    'sigma_L_err': sigma_L_err_list,
    'sigma_T': sigma_T_list,
    'sigma_T_err': sigma_T_err_list,
    'sigma_LT': sigma_LT_list,
    'sigma_LT_err': sigma_LT_err_list,
    'sigma_TT': sigma_TT_list,
    'sigma_TT_err': sigma_TT_err_list
})
sep_xsection_df.to_csv(sep_xsection_csv, index=False)
print(f"Saved separated cross-sections to {sep_xsection_csv}")

#============================================================================================================================================

# Extract Q2 value from PHY_SETTING argument
Q2par = 0
Q2str = PHY_SETTING.split('_')[0]
if Q2str.startswith('Q'):
    Q2par = Q2str[1:].replace('p', '')  # '3p85' -> '385' (string, no decimal)

# Read fit parameters from file (first column only)
param_file = "%s/iter%s/par.pl_%s" % (XSECT_PARAMPATH, ITERATION, Q2par)
params = []
with open(param_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 1:
            params.append(float(parts[0]))

# Create parameter dictionary for cleaner access
p_vals = {f'p{i+1}': params[i] for i in range(len(params))}

# Assign to individual variables (16 parameters expected)
p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = [p_vals.get(f'p{i}', 0.0) for i in range(1, 17)]
#p10, p11, p12, p13, p14 = -100000, 0, 0, 0, 0
#p10 = -9
print("Fit parameters loaded from file:")
print(f"p1={p1}, p2={p2}, p3={p3}, p4={p4}, p5={p5}, p6={p6}, p7={p7}, p8={p8}, p9={p9}, p10={p10}, p11={p11}, p12={p12}, p13={p13}, p14={p14}, p15={p15}, p16={p16}")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print("==========ROOT fit section==========")
# ROOT version of fitting functional form to data
LT_Sep_Xsection_root_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_ProdCoin_Pion_Analysis_LT_Sep_xsection_iter%s_Distributions_ROOT.pdf" % (REPLAYPATH, PHY_SETTING, ITERATION)

# Prepare arrays for -t, values, and errors (from fit results lists)
t_vals = sep_xsection_df['t_central'].values
w_vals = sep_xsection_df['W'].values
Q2_vals = np.array([sep_xsection_df[sep_xsection_df["t_central"] == var_tcm]["Q2"].iloc[0] for var_tcm in t_vals])
theta_vals = np.array([sep_xsection_df[sep_xsection_df["t_central"] == var_tcm]["theta_cm"].iloc[0] for var_tcm in t_vals])
m_p = 0.93827231
m_pi = 0.13957018
abs_t_cent_var = t_vals.astype(np.float64)
Q2_var = Q2_vals.astype(np.float64)
W_var = w_vals.astype(np.float64)
theta_cm_var = np.deg2rad(theta_vals.astype(np.float64))
sigma_T_var = sep_xsection_df['sigma_T'].values.astype(np.float64)
sigma_T_err = sep_xsection_df['sigma_T_err'].values.astype(np.float64)
sigma_L_var = sep_xsection_df['sigma_L'].values.astype(np.float64)
sigma_L_err = sep_xsection_df['sigma_L_err'].values.astype(np.float64)
sigma_LT_var = sep_xsection_df['sigma_LT'].values.astype(np.float64)
sigma_LT_err = sep_xsection_df['sigma_LT_err'].values.astype(np.float64)
sigma_TT_var = sep_xsection_df['sigma_TT'].values.astype(np.float64)
sigma_TT_err = sep_xsection_df['sigma_TT_err'].values.astype(np.float64)
npoints = len(abs_t_cent_var)
zeros_arr_f64 = np.zeros(npoints, dtype=np.float64)

graphs = {
    'T': ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_T_var, zeros_arr_f64, sigma_T_err),
    'L': ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_L_var, zeros_arr_f64, sigma_L_err),
    'LT': ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_LT_var, zeros_arr_f64, sigma_LT_err),
    'TT': ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_TT_var, zeros_arr_f64, sigma_TT_err)
}
gr_T = graphs['T']
gr_L = graphs['L']
gr_LT = graphs['LT']
gr_TT = graphs['TT']

# --- Chi-square quality thresholds (configurable) ---
CHI2_GOOD_THRESHOLD = 2.0      # Reduced chi2 < 2.0 is considered GOOD
CHI2_ACCEPTABLE_THRESHOLD = 5.0  # Reduced chi2 < 5.0 is considered ACCEPTABLE

# --- NEW: Parameter search configuration ---
MAX_SEARCH_ITERATIONS = 10000  # Maximum parameter search attempts per cross-section
PARAM_VARIATION_FACTOR = 1e6  # Parameter variation range: (1-factor, 1+factor) = (0.5, 1.5) for ±50%
RANDOM_SEARCH_ATTEMPTS = 5  # Number of random initial guesses to try

# --- Fit control flags: Set to True to fit, False to fix ---
FIT_SIGMA_L = False   # Fit σL parameters (p4, p5, p6, p7, p8)  
FIT_SIGMA_T = False  # Fit σT parameters (p1, p2, p3)
FIT_SIGMA_LT = False  # Fit σLT parameters (p9, p10, p11, p12)
FIT_SIGMA_TT = True  # Fit σTT parameters (p13, p14, p15, p16)

print("\nFit control settings:")
print(f"FIT_SIGMA_T  = {FIT_SIGMA_T}")
print(f"FIT_SIGMA_L  = {FIT_SIGMA_L}")
print(f"FIT_SIGMA_LT = {FIT_SIGMA_LT}")
print(f"FIT_SIGMA_TT = {FIT_SIGMA_TT}")

# --- NEW: Function to generate parameter variations ---
def generate_param_variations(base_params, n_variations, fit_enabled=True):
    # If fitting is disabled, only return the base parameters (no variations)
    if not fit_enabled:
        return [base_params]
    variations = [base_params]  # Include original as first variation
    # Generate n_variations additional parameter sets with random variations
    for _ in range(n_variations):
        varied = []
        for p in base_params:
            if p == 0:
                # For zero parameters, use small random values
                varied.append(random.uniform(-1, 1))
            else:
                # Use PARAM_VARIATION_FACTOR to determine variation range
                factor = random.uniform(1 - PARAM_VARIATION_FACTOR, 1 + PARAM_VARIATION_FACTOR)
                varied.append(p * factor)
        variations.append(varied)
    return variations

# --- Modified: Function to try fitting with direct parameter control ---
def fit_with_search(graph, fit_func, param_names, initial_params, function_fit_enabled, chi2_threshold=CHI2_GOOD_THRESHOLD):
    print(f"\n{'='*60}")
    print(f"Starting parameter search for {fit_func.GetName()}")
    
    # If function-level fitting is disabled, skip parameter search entirely
    if not function_fit_enabled:
        print("Function-level fitting DISABLED - Using fixed parameters")
        print(f"{'='*60}")
        
        # Perform single evaluation with fixed parameters (already set in main code)
        try:
            fit_result = graph.Fit(fit_func, "RMES")
            if fit_result and fit_result.IsValid() and fit_result.Ndf() > 0:
                chi2_reduced = fit_result.Chi2() / fit_result.Ndf()
                print(f"Fixed parameter evaluation: χ²/NDF = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {chi2_reduced:.3f}")
                return fit_result, initial_params[:], chi2_reduced, 1
            else:
                print("WARNING: Fixed parameter evaluation failed")
                chi2_reduced = 999999.0
                dummy_result = type('obj', (object,), {
                    'Chi2': lambda: 999999.0,
                    'Ndf': lambda: 1,
                    'Status': lambda: -1,
                    'IsValid': lambda: False
                })()
                return dummy_result, initial_params[:], chi2_reduced, 1
        except Exception as e:
            print(f"ERROR: Fixed parameter evaluation exception: {str(e)}")
            chi2_reduced = 999999.0
            dummy_result = type('obj', (object,), {
                'Chi2': lambda: 999999.0,
                'Ndf': lambda: 1,
                'Status': lambda: -1,
                'IsValid': lambda: False
            })()
            return dummy_result, initial_params[:], chi2_reduced, 1
    
    print(f"Function-level fitting ENABLED - Searching for optimal parameters")
    print(f"Target: Reduced χ² < {chi2_threshold}")
    print(f"{'='*60}")
    
    best_chi2_reduced = float('inf')
    best_fit_result = None
    best_params = initial_params[:]
    successful_attempts = 0
    
    # Generate parameter variations
    param_variations = generate_param_variations(initial_params, RANDOM_SEARCH_ATTEMPTS, fit_enabled=True)
    
    for attempt, params in enumerate(param_variations, 1):
        if attempt > MAX_SEARCH_ITERATIONS:
            print(f"Reached maximum iterations ({MAX_SEARCH_ITERATIONS})")
            break
            
        # Set parameters (parameter fixing is already handled in main code)
        for i, p in enumerate(params):
            fit_func.SetParameter(i, p)
        
        # Perform fit
        try:
            fit_result = graph.Fit(fit_func, "RMES")  # Q for quiet mode
            
            if fit_result and fit_result.IsValid() and fit_result.Status() == 0 and fit_result.Ndf() > 0:
                successful_attempts += 1
                chi2_reduced = fit_result.Chi2() / fit_result.Ndf()
                
                print(f"Attempt {attempt:3d}: χ²/NDF = {fit_result.Chi2():.3f}/{fit_result.Ndf()} = {chi2_reduced:.3f}", end='')
                
                # Update best result
                if chi2_reduced < best_chi2_reduced:
                    best_chi2_reduced = chi2_reduced
                    best_fit_result = fit_result
                    best_params = [fit_func.GetParameter(i) for i in range(len(params))]
                    print(" [NEW BEST]", end='')
                
                # Check if we achieved good fit
                if chi2_reduced < chi2_threshold:
                    print(" [SUCCESS - GOOD FIT]")
                    print(f"Final parameters: {[f'{p:.4g}' for p in best_params]}")
                    return best_fit_result, best_params, chi2_reduced, attempt
                
                print()
            else:
                status = fit_result.Status() if fit_result else "None"
                print(f"Attempt {attempt:3d}: Fit failed (status={status})")
        except Exception as e:
            print(f"Attempt {attempt:3d}: Exception during fit: {str(e)}")
    
    # Handle case where no successful fits occurred
    if best_fit_result is None or successful_attempts == 0:
        print(f"\nWARNING: All {attempt} fit attempts failed!")
        print("Performing one final fit with original parameters...")
        for i, p in enumerate(initial_params):
            fit_func.SetParameter(i, p)
        
        try:
            best_fit_result = graph.Fit(fit_func, "RQMES")  # No Q flag for visibility
            if best_fit_result and best_fit_result.IsValid() and best_fit_result.Ndf() > 0:
                best_chi2_reduced = best_fit_result.Chi2() / best_fit_result.Ndf()
                best_params = [fit_func.GetParameter(i) for i in range(len(initial_params))]
                print(f"Final fallback fit: χ² reduced = {best_chi2_reduced:.3f}")
            else:
                print("ERROR: Even fallback fit failed. Using original parameters with dummy chi2.")
                best_chi2_reduced = 999999.0
                best_params = initial_params[:]
                best_fit_result = type('obj', (object,), {
                    'Chi2': lambda: 999999.0,
                    'Ndf': lambda: 1,
                    'Status': lambda: -1,
                    'IsValid': lambda: False
                })()
        except Exception as e:
            print(f"ERROR: Fallback fit exception: {str(e)}")
            best_chi2_reduced = 999999.0
            best_params = initial_params[:]
            best_fit_result = type('obj', (object,), {
                'Chi2': lambda: 999999.0,
                'Ndf': lambda: 1,
                'Status': lambda: -1,
                'IsValid': lambda: False
            })()
        
        return best_fit_result, best_params, best_chi2_reduced, attempt
    
    # Return best result even if threshold not met
    print(f"\nBest result after {successful_attempts} successful attempts (out of {attempt} total): χ² reduced = {best_chi2_reduced:.3f}")
    if best_chi2_reduced < CHI2_ACCEPTABLE_THRESHOLD:
        print("Status: ACCEPTABLE (not optimal but usable)")
    else:
        print("Status: POOR (consider manual parameter adjustment)")
    
    # Restore best parameters
    for i, p in enumerate(best_params):
        fit_func.SetParameter(i, p)
    
    return best_fit_result, best_params, best_chi2_reduced, attempt

# --- Chi-square tracking and fitting for each component (T, L, LT, TT) ---
chi2_results = {
    'T': {'chi2': None, 'ndf': None, 'reduced_chi2': None},
    'L': {'chi2': None, 'ndf': None, 'reduced_chi2': None},
    'LT': {'chi2': None, 'ndf': None, 'reduced_chi2': None},
    'TT': {'chi2': None, 'ndf': None, 'reduced_chi2': None}
}

# Fit T vs -t
def sigma_T_root_func(x, par):
    var_tcm = x[0]
    Q2 = np.interp(var_tcm, abs_t_cent_var, Q2_var)
    W = np.interp(var_tcm, abs_t_cent_var, W_var)
    W_factor = 1.0 / (W**2 - m_p**2)**2
#    return (W_factor * ((par[0] / Q2) + (par[1] / (Q2**2))))
#    return (W_factor * (((par[0] / Q2) + (par[1] / (Q2**2))) * (np.exp(par[2] * np.abs(var_tcm)))))
    return(W_factor * ((par[0] / Q2) * (np.exp(par[1] * Q2**2)) * (np.exp(par[2] * np.abs(var_tcm)))))
sigma_T_root = ROOT.TF1("sigma_T_func", sigma_T_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 3)
#p3 = -0.5
sigma_T_root.SetParameters(p1, p2, p3)  # Set par[0] = p1, par[1] = p2, par[2] = p3
# Control parameter fixing for σT
if not FIT_SIGMA_T:
    sigma_T_root.FixParameter(0, p1)  # Fix p1
    sigma_T_root.FixParameter(1, p2)  # Fix p2
    sigma_T_root.FixParameter(2, p3)  # Fix p3
    print("σ_T parameters FIXED")
else:
    sigma_T_root.FixParameter(1, p2)  # Fix p2
    sigma_T_root.FixParameter(2, p3)  # Fix p3
    print("σ_T parameters will be FITTED")
sigma_T_root.SetLineColor(ROOT.kRed)
sigma_T_root.SetLineStyle(2)  # Dashed line

# --- NEW: Use parameter search for σT ---
initial_params_T = [p1, p2, p3]
fit_result_T, best_params_T, best_chi2_T, attempts_T = fit_with_search(
    gr_T, sigma_T_root, ['p1', 'p2', 'p3'], initial_params_T, FIT_SIGMA_T
)

# Safely extract chi2 results
if fit_result_T and hasattr(fit_result_T, 'Chi2'):
    chi2_results['T']['chi2'] = fit_result_T.Chi2()
    chi2_results['T']['ndf'] = fit_result_T.Ndf()
    chi2_results['T']['reduced_chi2'] = best_chi2_T
else:
    chi2_results['T']['chi2'] = 999999.0
    chi2_results['T']['ndf'] = 1
    chi2_results['T']['reduced_chi2'] = 999999.0
print(f"σ_T final: Chi2/NDF: {chi2_results['T']['chi2']:.3f}/{chi2_results['T']['ndf']} = {best_chi2_T:.3f} (attempts: {attempts_T})")

# Fit L vs -t
def sigma_L_root_func(x, par):
    var_tcm = x[0]
    Q2 = np.interp(var_tcm, abs_t_cent_var, Q2_var)
    W = np.interp(var_tcm, abs_t_cent_var, W_var)
    W_factor = 1.0 / (W**2 - m_p**2)**2
#    return (par[0] * Q2 * np.exp((par[1] - par[2] * np.log(Q2)) * np.abs(var_tcm))) / (1 + par[3] * Q2 + par[4] * (Q2**2))**2
#    return (W_factor * ((par[0] + par[1] / Q2) * np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2 * (np.exp(par[2] * np.abs(var_tcm))) *(Q2/((1 + par[3] * Q2 + par[4] * (Q2**2))**2))))
    return (W_factor * ((par[0] + par[1] / Q2) * (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2) * (np.exp(par[2] * np.abs(var_tcm))) *(Q2/((1 + par[3] * Q2 + par[4] * (Q2**2))**2))))
sigma_L_root = ROOT.TF1("sigma_L_func", sigma_L_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 5)
#p6 = 0.08
sigma_L_root.SetParameters(p4, p5, p6, p7, p8)  # All parameters are fitted
# Control parameter fixing for σL
if not FIT_SIGMA_L:
    sigma_L_root.FixParameter(0, p4)  # Fix p4
    sigma_L_root.FixParameter(1, p5)  # Fix p5
    sigma_L_root.FixParameter(2, p6)  # Fix p6
    sigma_L_root.FixParameter(3, p7)  # Fix p7
    sigma_L_root.FixParameter(4, p8)  # Fix p8
    print("σ_L parameters FIXED")
else:
    sigma_L_root.FixParameter(3, p7)  # Fix p7
    sigma_L_root.FixParameter(4, p8)  # Fix p8
    print("σ_L parameters will be FITTED")
sigma_L_root.SetLineColor(ROOT.kRed)
sigma_L_root.SetLineStyle(2)  # Dashed line

# --- NEW: Use parameter search for σL ---
initial_params_L = [p4, p5, p6, p7, p8]
fit_result_L, best_params_L, best_chi2_L, attempts_L = fit_with_search(
    gr_L, sigma_L_root, ['p4', 'p5', 'p6', 'p7', 'p8'], initial_params_L, FIT_SIGMA_L
)

# Safely extract chi2 results
if fit_result_L and hasattr(fit_result_L, 'Chi2'):
    chi2_results['L']['chi2'] = fit_result_L.Chi2()
    chi2_results['L']['ndf'] = fit_result_L.Ndf()
    chi2_results['L']['reduced_chi2'] = best_chi2_L
else:
    chi2_results['L']['chi2'] = 999999.0
    chi2_results['L']['ndf'] = 1
    chi2_results['L']['reduced_chi2'] = 999999.0
print(f"σ_L final: Chi2/NDF: {chi2_results['L']['chi2']:.3f}/{chi2_results['L']['ndf']} = {best_chi2_L:.3f} (attempts: {attempts_L})")

# Fit LT vs -t
def sigma_LT_root_func(x, par):
    var_tcm = x[0]
    Q2 = np.interp(var_tcm, abs_t_cent_var, Q2_var)
    theta = np.interp(var_tcm, abs_t_cent_var, theta_cm_var)
    W = np.interp(var_tcm, abs_t_cent_var, W_var)
    W_factor = 1.0 / (W**2 - m_p**2)**2
#    return (W_factor * (np.exp(par[0] + (par[1] * np.abs(var_tcm) / np.sqrt(Q2))) + par[2] + (par[3] / (Q2**2))) * np.sin(theta))
#    return (W_factor * ((((par[0] / (1 + Q2)) * np.exp(par[1] * np.abs(var_tcm))) + (par[2] / np.abs(var_tcm)**2)) * np.sin(theta)))
#    return (W_factor * ((((par[0] / (1 + (par[1] * Q2))) * (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2)) + (par[2] / np.abs(var_tcm)**2)) * np.sin(theta)))
#    return (W_factor * ((((par[0] / (Q2)) + np.exp(par[1] * np.abs(var_tcm))) * (par[2] / np.abs(var_tcm)**2)) * np.sin(theta)))
#    return (W_factor * (((par[0] / (Q2)) + (np.exp(par[1] * np.abs(var_tcm))) * (par[2] / (par[3] + np.abs(var_tcm))**2)) * np.sin(theta)))
#    return (W_factor * (((par[0] / (Q2)) + (np.exp(par[1] * np.abs(var_tcm))) * (par[2] / (np.abs(var_tcm))**par[3])) * np.sin(theta)))
#    return (W_factor * ((par[0] / (Q2)) * (np.exp(par[1] * np.abs(var_tcm))) * (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2) * np.sin(theta)))
#    return (W_factor * (((par[0] / (Q2)) + (par[1] / np.abs(var_tcm))) * np.sin(theta)))
#    return (W_factor * (((par[0] / (Q2)) + ((np.exp(-par[1] * np.abs(var_tcm))) * (par[2] / (par[3] + np.abs(var_tcm))))) * np.sin(theta)))
#    return (W_factor * (((par[0] / (Q2)) + ((np.exp(-par[1] * np.abs(var_tcm))) * (par[2] * np.abs(var_tcm)**par[3]))) * np.sin(theta)))
    return (W_factor * (((par[0] / (Q2)) + (par[1] / (np.abs(var_tcm))) + ((np.exp(par[2] * np.abs(var_tcm))) * (par[3] / (np.abs(var_tcm))**2))) * np.sin(theta)))
sigma_LT_root = ROOT.TF1("sigma_LT_func", sigma_LT_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 4)
#p11 = -15
#p11 = -0.003
sigma_LT_root.SetParameters(p9, p10, p11, p12)  # All parameters are fitted
# Control parameter fixing for σLT
if not FIT_SIGMA_LT:
    sigma_LT_root.FixParameter(0, p9)  # Fix p9
    sigma_LT_root.FixParameter(1, p10)  # Fix p10
    sigma_LT_root.FixParameter(2, p11)  # Fix p11
    sigma_LT_root.FixParameter(3, p12)  # Fix p12
    print("σ_LT parameters FIXED")
else:
    sigma_LT_root.FixParameter(2, p11)  # Fix p11
    print("σ_LT parameters will be FITTED")
sigma_LT_root.SetLineColor(ROOT.kRed)
sigma_LT_root.SetLineStyle(2)  # Dashed line

# --- NEW: Use parameter search for σLT ---
initial_params_LT = [p9, p10, p11, p12]
fit_result_LT, best_params_LT, best_chi2_LT, attempts_LT = fit_with_search(
    gr_LT, sigma_LT_root, ['p9', 'p10', 'p11', 'p12'], initial_params_LT, FIT_SIGMA_LT
)

# Safely extract chi2 results
if fit_result_LT and hasattr(fit_result_LT, 'Chi2'):
    chi2_results['LT']['chi2'] = fit_result_LT.Chi2()
    chi2_results['LT']['ndf'] = fit_result_LT.Ndf()
    chi2_results['LT']['reduced_chi2'] = best_chi2_LT
else:
    chi2_results['LT']['chi2'] = 999999.0
    chi2_results['LT']['ndf'] = 1
    chi2_results['LT']['reduced_chi2'] = 999999.0
print(f"σ_LT final: Chi2/NDF: {chi2_results['LT']['chi2']:.3f}/{chi2_results['LT']['ndf']} = {best_chi2_LT:.3f} (attempts: {attempts_LT})")

# Fit TT vs -t
def sigma_TT_root_func(x, par):
    var_tcm = x[0]
    Q2 = np.interp(var_tcm, abs_t_cent_var, Q2_var)
    theta = np.interp(var_tcm, abs_t_cent_var, theta_cm_var)
    W = np.interp(var_tcm, abs_t_cent_var, W_var)
    W_factor = 1.0 / (W**2 - m_p**2)**2
#    return (W_factor * ((par[0] / (Q2**2)) * (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2) * np.sin(theta)**2))
#    return (W_factor * ((par[0] / (Q2**2)) * (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2) * (np.exp(par[1] * np.abs(var_tcm))) * np.sin(theta)**2))
#    return (W_factor * ((((par[0] /( 1 + Q2)) * np.exp(par[1] * np.abs(var_tcm))) + (par[2] / np.abs(var_tcm)**3)) * np.sin(theta)**2))
#    return (W_factor * (((par[0] /(Q2)) + (np.abs(var_tcm) / (np.abs(var_tcm) + m_pi**2)**2) * (par[1] / np.abs(var_tcm)**3)) * np.sin(theta)**2))
#    return (W_factor * (((par[0] /(Q2))  + (par[1] / np.abs(var_tcm)**3)) * np.sin(theta)**2))
#    return (W_factor * (((par[0] /(Q2)) + (np.exp(par[1] * np.abs(var_tcm))) * (par[2] / (par[3] + np.abs(var_tcm))**3)) * np.sin(theta)**2))
#    return (W_factor * (((par[0] / (Q2)) + (par[1] / (np.abs(var_tcm))**2) + ((np.exp(par[2] * np.abs(var_tcm))) * (par[3] / (np.abs(var_tcm))**3))) * np.sin(theta)**2))
#    return (W_factor * (((par[0] /(Q2)) + (np.exp(par[1] * np.abs(var_tcm))) * (par[2] / (par[3] + np.abs(var_tcm))**3)) * np.sin(theta)**2))
#    return (W_factor * (((par[0] / (Q2)) + (par[1] / (np.abs(var_tcm))) + ((np.exp(par[2] / np.abs(var_tcm))) * (par[3] / (np.abs(var_tcm))**3))) * np.sin(2 * theta)))
#    return (W_factor * (((par[0] /(Q2)) + (np.exp(par[1] / np.abs(var_tcm))) * (par[2] / (par[3] + np.abs(var_tcm))**3)) * np.sin(theta)**2))
    return (W_factor * (((par[0] / (Q2)) + ((np.exp(-par[1] * np.abs(var_tcm))) * (par[2] / (np.abs(var_tcm))**par[3]))) * np.sin(theta)**2))
sigma_TT_root = ROOT.TF1("sigma_TT_func", sigma_TT_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 4)
#p13= 1
p14 = -11
p16 = 3
sigma_TT_root.SetParameters(p13, p14, p15, p16)  # All parameters are fitted
# Control parameter fixing for σTT
if not FIT_SIGMA_TT:
    sigma_TT_root.FixParameter(0, p13)  # Fix p13
    sigma_TT_root.FixParameter(1, p14)  # Fix p14
    sigma_TT_root.FixParameter(2, p15)  # Fix p15
    sigma_TT_root.FixParameter(3, p16)  # Fix p16
    print("σ_TT parameters FIXED")
else:
#    sigma_TT_root.FixParameter(0, p13)  # Fix p13
    sigma_TT_root.FixParameter(1, p14)  # Fix p14
#    sigma_TT_root.FixParameter(2, p15)  # Fix p15
    sigma_TT_root.FixParameter(3, p16)  # Fix p16
    print("σ_TT parameters will be FITTED")
sigma_TT_root.SetLineColor(ROOT.kRed)
sigma_TT_root.SetLineStyle(2)  # Dashed line

# --- NEW: Use parameter search for σTT ---
initial_params_TT = [p13, p14, p15, p16]
fit_result_TT, best_params_TT, best_chi2_TT, attempts_TT = fit_with_search(
    gr_TT, sigma_TT_root, ['p13', 'p14', 'p15', 'p16'], initial_params_TT, FIT_SIGMA_TT
)

# Safely extract chi2 results
if fit_result_TT and hasattr(fit_result_TT, 'Chi2'):
    chi2_results['TT']['chi2'] = fit_result_TT.Chi2()
    chi2_results['TT']['ndf'] = fit_result_TT.Ndf()
    chi2_results['TT']['reduced_chi2'] = best_chi2_TT
else:
    chi2_results['TT']['chi2'] = 999999.0
    chi2_results['TT']['ndf'] = 1
    chi2_results['TT']['reduced_chi2'] = 999999.0
print(f"σ_TT final: Chi2/NDF: {chi2_results['TT']['chi2']:.3f}/{chi2_results['TT']['ndf']} = {best_chi2_TT:.3f} (attempts: {attempts_TT})")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Print summary of chi-square results ---
print("\n" + "="*60)
print("CHI-SQUARE FIT QUALITY SUMMARY (AFTER PARAMETER SEARCH)")
print(f"Thresholds: GOOD < {CHI2_GOOD_THRESHOLD}, ACCEPTABLE < {CHI2_ACCEPTABLE_THRESHOLD}")
print("="*60)
for name in ['T', 'L', 'LT', 'TT']:
    result = chi2_results[name]
    status = "GOOD" if result['reduced_chi2'] < CHI2_GOOD_THRESHOLD else ("ACCEPTABLE" if result['reduced_chi2'] < CHI2_ACCEPTABLE_THRESHOLD else "POOR")
    attempts = {'T': attempts_T, 'L': attempts_L, 'LT': attempts_LT, 'TT': attempts_TT}[name]
    print(f"σ_{name:2s}: χ²/NDF = {result['chi2']:.3f}/{result['ndf']} = {result['reduced_chi2']:.3f} [{status}] (attempts: {attempts})")

# Identify best and worst fits
valid_fits = {k: v for k, v in chi2_results.items() if v['reduced_chi2'] != float('inf')}
if valid_fits:
    best_fit = min(valid_fits, key=lambda k: valid_fits[k]['reduced_chi2'])
    worst_fit = max(valid_fits, key=lambda k: valid_fits[k]['reduced_chi2'])
    print(f"\nBest fit:  σ_{best_fit} with reduced χ² = {chi2_results[best_fit]['reduced_chi2']:.3f}")
    print(f"Worst fit: σ_{worst_fit} with reduced χ² = {chi2_results[worst_fit]['reduced_chi2']:.3f}")
print("="*60 + "\n")

# Print fit results
LT_fit_rt = ROOT.TCanvas("LT_fit_rt", "Fitted Cross Sections vs -t", 1400, 900)
LT_fit_rt.Divide(2,2)

# sigma_T plot (top right - position 2)
LT_fit_rt.cd(2)
gr_T.SetMarkerStyle(20)
gr_T.SetMarkerSize(1.4)
gr_T.SetTitle(r'#sigma_{T} vs t; t (GeV^{2}); #sigma_{T} (#mub/GeV^{2})')
gr_T.Draw("AP")
gr_T.GetXaxis().SetTitle("t (GeV^{2})")
gr_T.GetYaxis().SetTitle("#sigma_{T} (#mub/GeV^{2})")
gr_T.GetXaxis().SetTitleSize(0.05)
gr_T.GetYaxis().SetTitleSize(0.05)
gr_T.GetXaxis().CenterTitle()
gr_T.GetYaxis().CenterTitle()
ROOT.gPad.SetTitle("")
sigma_T_root.Draw("SAME")

# sigma_L plot (top left - position 1)
LT_fit_rt.cd(1)
gr_L.SetMarkerStyle(20)
gr_L.SetMarkerSize(1.4)
gr_L.SetTitle(r'#sigma_{L} vs t; t (GeV^{2}); #sigma_{L} (#mub/GeV^{2})')
gr_L.Draw("AP")
gr_L.GetXaxis().SetTitle("t (GeV^{2})")
gr_L.GetYaxis().SetTitle("#sigma_{L} (#mub/GeV^{2})")
gr_L.GetXaxis().SetTitleSize(0.05)
gr_L.GetYaxis().SetTitleSize(0.05)
gr_L.GetXaxis().CenterTitle()
gr_L.GetYaxis().CenterTitle()
ROOT.gPad.SetTitle("")
sigma_L_root.Draw("SAME")

# sigma_LT plot (bottom left - position 3)
LT_fit_rt.cd(3)
gr_LT.SetMarkerStyle(20)
gr_LT.SetMarkerSize(1.4)
gr_LT.SetTitle(r'#sigma_{LT} vs t; t (GeV^{2}); #sigma_{LT} (#mub/GeV^{2})')
gr_LT.Draw("AP")
gr_LT.GetXaxis().SetTitle("t (GeV^{2})")
gr_LT.GetYaxis().SetTitle("#sigma_{LT} (#mub/GeV^{2})")
gr_LT.GetXaxis().SetTitleSize(0.05)
gr_LT.GetYaxis().SetTitleSize(0.05)
gr_LT.GetXaxis().CenterTitle()
gr_LT.GetYaxis().CenterTitle()
ROOT.gPad.SetTitle("")
sigma_LT_root.Draw("SAME")

# sigma_TT plot (bottom right - position 4)
LT_fit_rt.cd(4)
gr_TT.SetMarkerStyle(20)
gr_TT.SetMarkerSize(1.4)
gr_TT.SetTitle(r'#sigma_{TT} vs t; t (GeV^{2}); #sigma_{TT} (#mub/GeV^{2})')
gr_TT.Draw("AP")
gr_TT.GetXaxis().SetTitle("t (GeV^{2})")
gr_TT.GetYaxis().SetTitle("#sigma_{TT} (#mub/GeV^{2})")
gr_TT.GetXaxis().SetTitleSize(0.05)
gr_TT.GetYaxis().SetTitleSize(0.05)
gr_TT.GetXaxis().CenterTitle()
gr_TT.GetYaxis().CenterTitle()
ROOT.gPad.SetTitle("")
sigma_TT_root.Draw("SAME")
LT_fit_rt.Print(LT_Sep_Xsection_root_pdf + '(')

# --- Add a new page with 4 canvases, extended t range and fit lines ---
x_range = 0.07
extended_t_min = min(abs_t_cent_var) - x_range
extended_t_max = max(abs_t_cent_var) + x_range
y_min_L = -1
y_max_L = 1
y_min_T = 0
y_max_T = 1
y_min_LT = -0.5
y_max_LT = 0.5
y_min_TT = -0.5
y_max_TT = 0.5

canvas_ext = ROOT.TCanvas("canvas_ext", "Extended Range Cross Sections vs -t", 1400, 900)
canvas_ext.Divide(2,2)

# sigma_T extended plot (top right - position 2)
canvas_ext.cd(2)
gr_T.Draw("AP")
gr_T.GetXaxis().SetLimits(extended_t_min, extended_t_max)
gr_T.GetYaxis().SetRangeUser(y_min_T, y_max_T)
sigma_T_root.SetLineColor(ROOT.kBlue)
sigma_T_root.SetLineStyle(2)
ext_sigma_T = ROOT.TF1("ext_sigma_T_func", sigma_T_root_func, extended_t_min, extended_t_max, 3)
ext_sigma_T.SetParameters(*[sigma_T_root.GetParameter(i) for i in range(3)])
ext_sigma_T.SetLineColor(ROOT.kBlue)
ext_sigma_T.SetLineStyle(1)
ext_sigma_T.Draw("SAME")

# sigma_L extended plot (top left - position 1)
canvas_ext.cd(1)
gr_L.Draw("AP")
gr_L.GetXaxis().SetLimits(extended_t_min, extended_t_max)
gr_L.GetYaxis().SetRangeUser(y_min_L, y_max_L)
sigma_L_root.SetLineColor(ROOT.kBlue)
sigma_L_root.SetLineStyle(2)
ext_sigma_L = ROOT.TF1("ext_sigma_L_func", sigma_L_root_func, extended_t_min, extended_t_max, 5)
ext_sigma_L.SetParameters(*[sigma_L_root.GetParameter(i) for i in range(5)])
ext_sigma_L.SetLineColor(ROOT.kBlue)
ext_sigma_L.SetLineStyle(1)
ext_sigma_L.Draw("SAME")

# sigma_LT extended plot (bottom left - position 3)
canvas_ext.cd(3)
gr_LT.Draw("AP")
gr_LT.GetXaxis().SetLimits(extended_t_min, extended_t_max)
gr_LT.GetYaxis().SetRangeUser(y_min_LT, y_max_LT)
sigma_LT_root.SetLineColor(ROOT.kBlue)
sigma_LT_root.SetLineStyle(2)
ext_sigma_LT = ROOT.TF1("ext_sigma_LT_func", sigma_LT_root_func, extended_t_min, extended_t_max, 4)
ext_sigma_LT.SetParameters(*[sigma_LT_root.GetParameter(i) for i in range(4)])
#ext_sigma_LT.SetParameter(0, sigma_LT_root.GetParameter(0))
ext_sigma_LT.SetLineColor(ROOT.kBlue)
ext_sigma_LT.SetLineStyle(1)
ext_sigma_LT.Draw("SAME")

# sigma_TT extended plot (bottom right - position 4)
canvas_ext.cd(4)
gr_TT.Draw("AP")
gr_TT.GetXaxis().SetLimits(extended_t_min, extended_t_max)
gr_TT.GetYaxis().SetRangeUser(y_min_TT, y_max_TT)
sigma_TT_root.SetLineColor(ROOT.kBlue)
sigma_TT_root.SetLineStyle(2)
ext_sigma_TT = ROOT.TF1("ext_sigma_TT_func", sigma_TT_root_func, extended_t_min, extended_t_max, 4)
ext_sigma_TT.SetParameters(*[sigma_TT_root.GetParameter(i) for i in range(4)])
ext_sigma_TT.SetLineColor(ROOT.kBlue)
ext_sigma_TT.SetLineStyle(1)
ext_sigma_TT.Draw("SAME")
canvas_ext.Print(LT_Sep_Xsection_root_pdf + ')')
print(f"Saved ROOT fits to {LT_Sep_Xsection_root_pdf}")

# Collect all parameters and errors in order, reflecting which are fitted/fixed
params = [
    (best_params_T[0], sigma_T_root.GetParError(0), 'fitted' if FIT_SIGMA_T else 'fixed'),   # p1
    (best_params_T[1], sigma_T_root.GetParError(1), 'fitted' if FIT_SIGMA_T else 'fixed'),   # p2
    (best_params_T[2], sigma_T_root.GetParError(2), 'fitted' if FIT_SIGMA_T else 'fixed'),   # p3

    (best_params_L[0], sigma_L_root.GetParError(0), 'fitted' if FIT_SIGMA_L else 'fixed'),   # p4
    (best_params_L[1], sigma_L_root.GetParError(1), 'fitted' if FIT_SIGMA_L else 'fixed'),   # p5
    (best_params_L[2], sigma_L_root.GetParError(2), 'fitted' if FIT_SIGMA_L else 'fixed'),   # p6
    (best_params_L[3], sigma_L_root.GetParError(3), 'fitted' if FIT_SIGMA_L else 'fixed'),   # p7
    (best_params_L[4], sigma_L_root.GetParError(4), 'fitted' if FIT_SIGMA_L else 'fixed'),   # p8

    (best_params_LT[0], sigma_LT_root.GetParError(0), 'fitted' if FIT_SIGMA_LT else 'fixed'), # p9
    (best_params_LT[1], sigma_LT_root.GetParError(1), 'fitted' if FIT_SIGMA_LT else 'fixed'), # p10
    (best_params_LT[2], sigma_LT_root.GetParError(2), 'fitted' if FIT_SIGMA_LT else 'fixed'), # p11
    (best_params_LT[3], sigma_LT_root.GetParError(3), 'fitted' if FIT_SIGMA_LT else 'fixed'), # p12

    (best_params_TT[0], sigma_TT_root.GetParError(0), 'fitted' if FIT_SIGMA_TT else 'fixed'), # p13
    (best_params_TT[1], sigma_TT_root.GetParError(1), 'fitted' if FIT_SIGMA_TT else 'fixed'), # p14
    (best_params_TT[2], sigma_TT_root.GetParError(2), 'fitted' if FIT_SIGMA_TT else 'fixed'), # p15
    (best_params_TT[3], sigma_TT_root.GetParError(3), 'fitted' if FIT_SIGMA_TT else 'fixed')  # p16
]

# Check for any problematic fits and warn user
fit_warnings = []
if not fit_result_T or not hasattr(fit_result_T, 'Status') or fit_result_T.Status() != 0:
    fit_warnings.append(f"σ_T fit status: {fit_result_T.Status() if fit_result_T and hasattr(fit_result_T, 'Status') else 'FAILED'}")
if not fit_result_L or not hasattr(fit_result_L, 'Status') or fit_result_L.Status() != 0:
    fit_warnings.append(f"σ_L fit status: {fit_result_L.Status() if fit_result_L and hasattr(fit_result_L, 'Status') else 'FAILED'}")
if not fit_result_LT or not hasattr(fit_result_LT, 'Status') or fit_result_LT.Status() != 0:
    fit_warnings.append(f"σ_LT fit status: {fit_result_LT.Status() if fit_result_LT and hasattr(fit_result_LT, 'Status') else 'FAILED'}")
if not fit_result_TT or not hasattr(fit_result_TT, 'Status') or fit_result_TT.Status() != 0:
    fit_warnings.append(f"σ_TT fit status: {fit_result_TT.Status() if fit_result_TT and hasattr(fit_result_TT, 'Status') else 'FAILED'}")

if fit_warnings:
    print("\n" + "!"*60)
    print("FIT WARNINGS:")
    for warning in fit_warnings:
        print(f"  - {warning}")
    print("!"*60 + "\n")

# Print fit results
print("Fitted and fixed parameters (AFTER PARAMETER SEARCH):")
for i, (val, err, status) in enumerate(params, start=1):
    print(f"p{i} = {val:.6f} ± {err:.6f} [{status}]")

# --- Save new fit parameters and errors to a text file in alphabetical order ---
param_names = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15', 'p16']
NEW_ITERATION = str(int(ITERATION) + 1).zfill(2)
param_outfile = f"{XSECT_OUTPATH}/new_fitparams_iter{NEW_ITERATION}_par.pl_{Q2par}"
with open(param_outfile, 'w') as f:
    for idx, (val, err, status) in enumerate(params, 1):
        f.write(f"{val:12.8g} {err:12.4g} {idx:3d}\n")
print(f"Saved new fit parameters to {param_outfile}")

print("Processing Complete")