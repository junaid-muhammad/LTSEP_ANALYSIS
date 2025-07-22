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

np.bool = bool
np.float = float

import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
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
def rosenbluth_fit(phi, sigma_T, sigma_L, sigma_LT, sigma_TT, eps):
    phi_rad = np.deg2rad(phi)  # Convert phi to radians for np.cos
    A = eps * sigma_L + sigma_T
    B = np.sqrt(2 * eps * (eps + 1)) * sigma_LT
    C = eps * sigma_TT
    return (1 / (2 * np.pi)) * (A + B * np.cos(phi_rad) + C * np.cos(2 * phi_rad))
    
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

        # Initial guess for [sigma_T, sigma_L, sigma_LT, sigma_TT]
        p0 = [y_vals.mean(), y_vals.mean()/2, 0, 0]
        try:
            popt, pcov = curve_fit(
                rosenbluth_fit, x_vals, y_vals, sigma=y_errs, p0=p0, absolute_sigma=True, maxfev=50000
            )
            sigma_T, sigma_L, sigma_LT, sigma_TT = popt
            perr = np.sqrt(np.diag(pcov))
            # Store all fit results for plotting
            abs_t_cent_list.append(t_central)
            sigma_L_list.append(sigma_L)
            sigma_L_err_list.append(perr[1])
            sigma_T_list.append(sigma_T)
            sigma_T_err_list.append(perr[0])
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
                    y_fit = rosenbluth_fit((phi_fit, np.full_like(phi_fit, eps_val)), sigma_T, sigma_L, sigma_LT, sigma_TT)
                    ax.plot(phi_fit, y_fit, color=color, linestyle='--', label=label)
            # Print cross-section fit results and residuals
            print(f"t_central = {t_central:.3f}")
            print(f"  sigma_L  = {sigma_L:.4f} ± {perr[1]:.4f}")
            print(f"  sigma_T  = {sigma_T:.4f} ± {perr[0]:.4f}")
            print(f"  sigma_LT = {sigma_LT:.4f} ± {perr[2]:.4f}")
            print(f"  sigma_TT = {sigma_TT:.4f} ± {perr[3]:.4f}")
            # Calculate chi^2/ndof for the fit
            y_fit_all = rosenbluth_fit((phi_vals, eps_vals), sigma_T, sigma_L, sigma_LT, sigma_TT)
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
                r"$\sigma_L$ = %.3f $\pm$ %.3f" % (sigma_L, perr[1]) + "\n"
                r"$\sigma_T$ = %.3f $\pm$ %.3f" % (sigma_T, perr[0]) + "\n"
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
    ax.set_ylim(0, max(df_t_low["Xsection_real"].max(), df_t_high["Xsection_real"].max()) * 1.1)
#    ax.set_ylim(-0.10, 0.10)
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

# --- fit functions for cross sections ---
def sigma_T_func(Q2_var, p1, p2):
    return (p1 / Q2_var) + (p2 / (Q2_var**2))
# p3 and p4 are unassigned parameters in the code so far.

def sigma_L_func(abs_t_cent_var, Q2_var, p5, p6, p7, p8, p9):
    return (p5 * Q2_var * np.exp((p6 - p7 * np.log(Q2_var)) * np.abs(abs_t_cent_var))) / (1 + p8 * Q2_var + p9 * (Q2_var**2))**2

def sigma_LT_func(abs_t_cent_var, Q2_var, p10, p11, p12, p13, theta_cm_var):
    return (np.exp(p10 + (p11 * np.abs(abs_t_cent_var) / np.sqrt(Q2_var))) + p12 + (p13 / (Q2_var**2))) * np.sin(theta_cm_var)

def sigma_TT_func(abs_t_cent_var, Q2_var, p14, m_pi, theta_cm_var):
    return ((p14 / (Q2_var**2)) * (np.abs(abs_t_cent_var) / (np.abs(abs_t_cent_var) + m_pi**2)**2) * np.sin(theta_cm_var)**2)
# p15 and p16 are unassigned parameters in the code so far.

# Extract Q2 value from PHY_SETTING argument
Q2par = 0
Q2str = PHY_SETTING.split('_')[0]
if Q2str.startswith('Q'):
    Q2par = Q2str[1:].replace('p', '')  # '3p85' -> '385' (string, no decimal)

# Read fit parameters from file and assign to variables (first column only)
param_file = "%s/iter%s/par.pl_%s" % (XSECT_PARAMPATH, ITERATION, Q2par)
params = []
with open(param_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 1:
            params.append(float(parts[0]))

# Assign parameters to variables (first value of each row)
p1 = params[0]  #1
p2 = params[1]  #2
p3 = params[2]  #3 - unassigned
p4 = params[3]  #4 - unassigned
p5 = params[4]  #5
p6 = params[5]  #6
p7 = params[6]  #7
p8 = params[7]  #8
p9 = params[8]  #9
p10 = params[9]  #10
p11 = params[10] #11
p12 = params[11] #12
p13 = params[12] #13
p14 = params[13] #14
p15 = params[14] #15 - unassigned
p16 = params[15] #16 - unassigned
print ("Fit parameters loaded from file:")
print(f"p1={p1}, p2={p2}, p3={p3}, p4={p4}, p5={p5}, p6={p6}, p7={p7}, p8={p8}, p9={p9}, p10={p10}, p11={p11}, p12={p12}, p13={p13}, p14={p14}, p15={p15}, p16={p16}")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ROOT version of fitting functional form to data
LT_Sep_Xsection_root_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_ProdCoin_Pion_Analysis_LT_Sep_xsection_iter%s_Distributions_ROOT.pdf" % (REPLAYPATH, PHY_SETTING, ITERATION)

# Prepare arrays for -t, values, and errors (from fit results lists)
Q2_var = sep_xsection_df['Q2'].values
theta_cm_var = np.deg2rad(sep_xsection_df['theta_cm'].values)
m_pi = 0.13957018
abs_t_cent_var = np.array(abs_t_cent_list)
sigma_L_var = np.array(sigma_L_list)
sigma_L_err_var = np.array(sigma_L_err_list)
sigma_T_var = np.array(sigma_T_list)
sigma_T_err_var = np.array(sigma_T_err_list)
sigma_LT_var = np.array(sigma_LT_list)
sigma_LT_err_var = np.array(sigma_LT_err_list)
sigma_TT_var = np.array(sigma_TT_list)
sigma_TT_err_var = np.array(sigma_TT_err_list)

npoints = len(abs_t_cent_var)
zeros_arr = np.zeros(npoints)

# Create TGraphErrors objects with error bars (using pre-converted arrays)
gr_T = ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_T_var, zeros_arr, sigma_T_err_var)
gr_L = ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_L_var, zeros_arr, sigma_L_err_var)
gr_LT = ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_LT_var, zeros_arr, sigma_LT_err_var)
gr_TT = ROOT.TGraphErrors(npoints, abs_t_cent_var, sigma_TT_var, zeros_arr, sigma_TT_err_var)

# Fit T vs -t
def sigma_T_root_func(x, par):
    var_f = x[0]
    Q2 = np.interp(var_f, abs_t_cent_var, Q2_var)
    return (par[0] / Q2) + (par[1] / (Q2**2))
sigma_T_root = ROOT.TF1("sigma_T_func", sigma_T_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 2)
sigma_T_root.SetParameters(p1, p2)  # Set par[0] = p1, par[1] = p2
sigma_T_root.SetLineColor(ROOT.kRed)
sigma_T_root.SetLineStyle(2)  # Dashed line
fit_result_T = gr_T.Fit(sigma_T_root, "RS")  # S option for better statistics
print(f"σ_T fit status: {fit_result_T.Status()}, Chi2/NDF: {fit_result_T.Chi2():.3f}/{fit_result_T.Ndf()}")

# Fit L vs -t
def sigma_L_root_func(x, par):
    var_f = x[0]
    Q2 = np.interp(var_f, abs_t_cent_var, Q2_var)
    return (par[0] * Q2 * np.exp((par[1] - par[2] * np.log(Q2)) * np.abs(var_f))) / (1 + par[3] * Q2 + par[4] * (Q2**2))**2
sigma_L_root = ROOT.TF1("sigma_L_func", sigma_L_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 5)
sigma_L_root.SetParameters(p5, p6, p7, p8, p9)  # Set par[0] = p5, par[1] = p6, par[2] = p7, par[3] = p8, par[4] = p9
sigma_L_root.SetLineColor(ROOT.kRed)
sigma_L_root.SetLineStyle(2)  # Dashed line
fit_result_L = gr_L.Fit(sigma_L_root, "RS")  # S option for better statistics
print(f"σ_L fit status: {fit_result_L.Status()}, Chi2/NDF: {fit_result_L.Chi2():.3f}/{fit_result_L.Ndf()}")

# Fit LT vs -t
def sigma_LT_root_func(x, par):
    var_f = x[0]
    Q2 = np.interp(var_f, abs_t_cent_var, Q2_var)
    theta = np.interp(var_f, abs_t_cent_var, theta_cm_var)
    return (np.exp(par[0] + (par[1] * np.abs(var_f) / np.sqrt(Q2))) + par[2] + (par[3] / (Q2**2))) * np.sin(theta)
sigma_LT_root = ROOT.TF1("sigma_LT_func", sigma_LT_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 4)
sigma_LT_root.SetParameters(p10, p11, p12, p13)  # Set par[0] = p10, par[1] = p11, par[2] = p12, par[3] = p13
sigma_LT_root.SetLineColor(ROOT.kRed)
sigma_LT_root.SetLineStyle(2)  # Dashed line
fit_result_LT = gr_LT.Fit(sigma_LT_root, "RS")  # S option for better statistics
print(f"σ_LT fit status: {fit_result_LT.Status()}, Chi2/NDF: {fit_result_LT.Chi2():.3f}/{fit_result_LT.Ndf()}")

# Fit TT vs -t
def sigma_TT_root_func(x, par):
    var_f = x[0]
    Q2 = np.interp(var_f, abs_t_cent_var, Q2_var)
    theta = np.interp(var_f, abs_t_cent_var, theta_cm_var)
    return ((par[0] / (Q2**2)) * (np.abs(var_f) / (np.abs(var_f) + m_pi**2)**2) * np.sin(theta)**2)
sigma_TT_root = ROOT.TF1("sigma_TT_func", sigma_TT_root_func, min(abs_t_cent_var), max(abs_t_cent_var), 1)
sigma_TT_root.SetParameter(0, p14)  # Set par[0] = p14
sigma_TT_root.SetLineColor(ROOT.kRed)
sigma_TT_root.SetLineStyle(2)  # Dashed line
fit_result_TT = gr_TT.Fit(sigma_TT_root, "RS")  # S option for better statistics
print(f"σ_TT fit status: {fit_result_TT.Status()}, Chi2/NDF: {fit_result_TT.Chi2():.3f}/{fit_result_TT.Ndf()}")

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

LT_fit_rt.Update()
LT_fit_rt.SaveAs(LT_Sep_Xsection_root_pdf)
print(f"Saved ROOT fits to {LT_Sep_Xsection_root_pdf}")

# Collect all parameters and errors in order
params = [
    (sigma_T_root.GetParameter(0), sigma_T_root.GetParError(0)),   # p1
    (sigma_T_root.GetParameter(1), sigma_T_root.GetParError(1)),   # p2
    # p3 and p4 are not fitted in this code, so keep original values
    (p3, 0.0),               # p3
    (p4, 0.0),               # p4
    (sigma_L_root.GetParameter(0), sigma_L_root.GetParError(0)),   # p5
    (sigma_L_root.GetParameter(1), sigma_L_root.GetParError(1)),   # p6
    (sigma_L_root.GetParameter(2), sigma_L_root.GetParError(2)),   # p7
    (sigma_L_root.GetParameter(3), sigma_L_root.GetParError(3)),   # p8
    (sigma_L_root.GetParameter(4), sigma_L_root.GetParError(4)),   # p9
    (sigma_LT_root.GetParameter(0), sigma_LT_root.GetParError(0)), # p10
    (sigma_LT_root.GetParameter(1), sigma_LT_root.GetParError(1)), # p11
    (sigma_LT_root.GetParameter(2), sigma_LT_root.GetParError(2)), # p12
    (sigma_LT_root.GetParameter(3), sigma_LT_root.GetParError(3)), # p13
    (sigma_TT_root.GetParameter(0), sigma_TT_root.GetParError(0)), # p14
    (p15, 0.0),              # p15
    (p16, 0.0)               # p16
]

# Check for any problematic fits and warn user
fit_warnings = []
if fit_result_T.Status() != 0:
    fit_warnings.append(f"σ_T fit status: {fit_result_T.Status()}")
if fit_result_L.Status() != 0:
    fit_warnings.append(f"σ_L fit status: {fit_result_L.Status()}")
if fit_result_LT.Status() != 0:
    fit_warnings.append(f"σ_LT fit status: {fit_result_LT.Status()}")
if fit_result_TT.Status() != 0:
    fit_warnings.append(f"σ_TT fit status: {fit_result_TT.Status()}")

print("Fitted parameters and uncertainties:")
for i, (val, err) in enumerate(params, start=1):
    print(f"p{i} = {val:.6f} ± {err:.6f}")

# --- Save new fit parameters and errors to a text file in alphabetical order ---
param_names = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15', 'p16']
param_values = [
    sigma_T_root.GetParameter(0), sigma_T_root.GetParameter(1),  # p1, p2
    p3, p4,  # Use loaded values for p3, p4 (not fitted)
    sigma_L_root.GetParameter(0), sigma_L_root.GetParameter(1), sigma_L_root.GetParameter(2), sigma_L_root.GetParameter(3), sigma_L_root.GetParameter(4),  # p5, p6, p7, p8, p9
    sigma_LT_root.GetParameter(0), sigma_LT_root.GetParameter(1), sigma_LT_root.GetParameter(2), sigma_LT_root.GetParameter(3),  # p10, p11, p12, p13
    sigma_TT_root.GetParameter(0),  # p14
    p15, p16   # Use loaded values for p15, p16 (not fitted)
]
param_errors = [
    sigma_T_root.GetParError(0), sigma_T_root.GetParError(1),  # p1, p2 errors
    0, 0,  # p3, p4 not fitted
    sigma_L_root.GetParError(0), sigma_L_root.GetParError(1), sigma_L_root.GetParError(2), sigma_L_root.GetParError(3), sigma_L_root.GetParError(4),  # p5, p6, p7, p8, p9 errors
    sigma_LT_root.GetParError(0), sigma_LT_root.GetParError(1), sigma_LT_root.GetParError(2), sigma_LT_root.GetParError(3),  # p10, p11, p12, p13 errors
    sigma_TT_root.GetParError(0),  # p14 error
    0, 0   # p15, p16 not fitted
]
# Write to file in numerical order
NEW_ITERATION = str(int(ITERATION) + 1).zfill(2)
param_outfile = f"{XSECT_OUTPATH}/new_fitparams_iter{NEW_ITERATION}_par.pl_{Q2par}"
with open(param_outfile, 'w') as f:
    for idx, (name, val, err) in enumerate(zip(param_names, param_values, param_errors), 1):
        f.write(f"{val:12.8g}{' '*7}{err:12.4g}{' '*5}{idx}\n")
print(f"Saved new fit parameters to {param_outfile}")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Python version of fitting functional form to data
LT_Sep_Xsection_py_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_ProdCoin_Pion_Analysis_LT_Sep_xsection_iter%s_Distributions_py.pdf" % (REPLAYPATH, PHY_SETTING, ITERATION)

# sigma_L fit
p0_L = [p5, p6, p7, p8, p9]  # Fit all parameters
def sigma_L_fitfunc(x_tuple, p5, p6, p7, p8, p9):
    t, Q2 = x_tuple
    return sigma_L_func(t, Q2, p5, p6, p7, p8, p9)
x_input_L = np.vstack((abs_t_cent_var, Q2_var))
popt_L, pcov_L = curve_fit(
    sigma_L_fitfunc, x_input_L, sigma_L, sigma=sigma_L_err, p0=p0_L,
    bounds=(-1e8, 1e8), absolute_sigma=True, maxfev=100000
)

# sigma_T fit
p0_T = [p1, p2]
def sigma_T_fitfunc(Q2, p1, p2):
    return sigma_T_func(Q2, p1, p2)
popt_T, pcov_T = curve_fit(
    sigma_T_fitfunc, Q2_var, sigma_T, sigma=sigma_T_err, p0=p0_T,
    absolute_sigma=True, maxfev=100000
)

# sigma_LT fit
p0_LT = [p10, p11, p12, p13]
def sigma_LT_fitfunc(x_tuple, p10, p11, p12, p13):
    t, Q2, theta_cm = x_tuple
    return sigma_LT_func(t, Q2, p10, p11, p12, p13, theta_cm)
x_input_LT = np.vstack((abs_t_cent_var, Q2_var, theta_cm_var))
popt_LT, pcov_LT = curve_fit(
    sigma_LT_fitfunc, x_input_LT, sigma_LT, sigma=sigma_LT_err, p0=p0_LT,
    absolute_sigma=True, maxfev=100000
)

# sigma_TT fit
p0_TT = [p14]
def sigma_TT_fitfunc(x_tuple, p14):
    t, Q2, theta_cm = x_tuple
    return sigma_TT_func(t, Q2, p14, m_pi, theta_cm)
x_input_TT = np.vstack((abs_t_cent_var, Q2_var, theta_cm_var))
popt_TT, pcov_TT = curve_fit(
    sigma_TT_fitfunc, x_input_TT, sigma_TT, sigma=sigma_TT_err, p0=p0_TT,
    absolute_sigma=True, maxfev=100000
)

# --- zPlot all four sigma's vs t ---
fig_3, axs_3 = plt.subplots(2, 2, figsize=(14, 10))
axs_3 = axs_3.flatten()

# σ_L vs t
axs_3[0].errorbar(abs_t_cent_var, sigma_L, yerr=sigma_L_err, fmt='o', color='blue', label='Data')
axs_3[0].set_title(r'$\sigma_L$ vs $t$', fontsize=16, fontweight='bold')
axs_3[0].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs_3[0].set_ylabel(r'$\sigma_L$ ($\mu$b/GeV$^2$)', fontsize=14)
# Functional fit for sigma_L
sigma_L_fit = sigma_L_func(abs_t_cent_var, Q2_var, *popt_L)
axs_3[0].plot(abs_t_cent_var, sigma_L_fit, 'r--', label='Fit')
axs_3[0].legend()

# σ_T vs t
axs_3[1].errorbar(abs_t_cent_var, sigma_T, yerr=sigma_T_err, fmt='s', color='red', label='Data')
axs_3[1].set_title(r'$\sigma_T$ vs $t$', fontsize=16, fontweight='bold')
axs_3[1].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs_3[1].set_ylabel(r'$\sigma_T$ ($\mu$b/GeV$^2$)', fontsize=14)
# Functional fit for sigma_T
sigma_T_fit = sigma_T_func(Q2_var, *popt_T)
axs_3[1].plot(abs_t_cent_var, sigma_T_fit, 'r--', label='Fit')
axs_3[1].legend()

# σ_LT vs t
axs_3[2].errorbar(abs_t_cent_var, sigma_LT, yerr=sigma_LT_err, fmt='^', color='green', label='Data')
axs_3[2].set_title(r'$\sigma_{LT}$ vs $t$', fontsize=16, fontweight='bold')
axs_3[2].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs_3[2].set_ylabel(r'$\sigma_{LT}$ ($\mu$b/GeV$^2$)', fontsize=14)
# Functional fit for sigma_LT
sigma_LT_fit = sigma_LT_func(abs_t_cent_var, Q2_var, *popt_LT, theta_cm_var)
axs_3[2].plot(abs_t_cent_var, sigma_LT_fit, 'r--', label='Fit')
axs_3[2].legend()

# σ_TT vs t
axs_3[3].errorbar(abs_t_cent_var, sigma_TT, yerr=sigma_TT_err, fmt='v', color='purple', label='Data')
axs_3[3].set_title(r'$\sigma_{TT}$ vs $t$', fontsize=16, fontweight='bold')
axs_3[3].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs_3[3].set_ylabel(r'$\sigma_{TT}$ ($\mu$b/GeV$^2$)', fontsize=14)
# Functional fit for sigma_TT
sigma_TT_fit = sigma_TT_func(abs_t_cent_var, Q2_var, *popt_TT, m_pi, theta_cm_var)
axs_3[3].plot(abs_t_cent_var, sigma_TT_fit, 'r--', label='Fit')
axs_3[3].legend()

plt.tight_layout()

# --- Save both figures to the same PDF ---
fig_3.savefig(LT_Sep_Xsection_py_pdf)
plt.close(fig_3)
print(f"Saved Sep cross-section plots to {LT_Sep_Xsection_py_pdf}")

print ("Processing Complete")