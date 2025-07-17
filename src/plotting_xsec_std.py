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
import array
import csv
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
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
XSECT_PATH = "%s/scripts/ltsep_analysis/src/xsects" % (UTILPATH)
XSECT_OUTPATH = "%s/scripts/ltsep_analysis/src/output" % (UTILPATH)

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
UnSep_Xsection_pdf = f"%s/%s_ProdCoin_Pion_Analysis_UnSep_xsection_iter%s_Distributions.pdf" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)

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
            capsize=4,
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
            capsize=4,
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
def rosenbluth_fit(phi_eps, sigma_T, sigma_L, sigma_LT, sigma_TT):
    phi, eps = phi_eps
    phi_rad = np.deg2rad(phi)  # Convert phi to radians for np.cos
    A = eps * sigma_L + sigma_T
    B = np.sqrt(2 * eps * (eps + 1)) * sigma_LT
    C = eps * sigma_TT
    return (1 / (2 * np.pi)) * (A + B * np.cos(phi_rad) + C * np.cos(2 * phi_rad))
    
# Output file path
Sep_Xsection_pdf = f"%s/%s_ProdCoin_Pion_Analysis_Sep_xsection_iter%s_Distributions.pdf" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)

# --- Prepare lists to collect fit results ---
sigma_L_list, sigma_L_err_list = [], []
sigma_T_list, sigma_T_err_list = [], []
sigma_LT_list, sigma_LT_err_list = [], []
sigma_TT_list, sigma_TT_err_list = [], []
t_plot_list = []

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
            capsize=4,
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
            capsize=4,
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
                rosenbluth_fit, x_vals, y_vals, sigma=y_errs, p0=p0, absolute_sigma=True, maxfev=10000
            )
            sigma_T, sigma_L, sigma_LT, sigma_TT = popt
            perr = np.sqrt(np.diag(pcov))
            # Store all fit results for plotting
            t_plot_list.append(t_central)
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
    ax.legend(fontsize=12)
plt.tight_layout()

# --- Save separated cross-sections (sigma_L, sigma_T, sigma_LT, sigma_TT) vs t, theta_cm, Q2 to CSV ---
sep_xsection_csv = "%s/%s_pion_physics_sep_xsect_iter%s.csv" % (XSECT_OUTPATH, PHY_SETTING, ITERATION)
sep_xsection_df = pd.DataFrame({
    't_central': t_plot_list,
    'theta_cm': [xsect_loweps_df[xsect_loweps_df['t_central'] == t]['theta_cm'].iloc[0] if not xsect_loweps_df[xsect_loweps_df['t_central'] == t].empty else (xsect_higheps_df[xsect_higheps_df['t_central'] == t]['theta_cm'].iloc[0] if not xsect_higheps_df[xsect_higheps_df['t_central'] == t].empty else None) for t in t_plot_list],
    'Q2': [xsect_loweps_df[xsect_loweps_df['t_central'] == t]['Q2'].iloc[0] if not xsect_loweps_df[xsect_loweps_df['t_central'] == t].empty else (xsect_higheps_df[xsect_higheps_df['t_central'] == t]['Q2'].iloc[0] if not xsect_higheps_df[xsect_higheps_df['t_central'] == t].empty else None) for t in t_plot_list],
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

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Plot all four sigma's vs t ---
fig4, axs4 = plt.subplots(2, 2, figsize=(14, 10))
axs4 = axs4.flatten()

# σ_L vs t
axs4[0].errorbar(t_plot_list, sigma_L_list, yerr=sigma_L_err_list, fmt='o', color='blue')
axs4[0].set_title(r'$\sigma_L$ vs $t$', fontsize=16, fontweight='bold')
axs4[0].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs4[0].set_ylabel(r'$\sigma_L$ ($\mu$b/GeV$^2$)', fontsize=14)

# σ_T vs t
axs4[1].errorbar(t_plot_list, sigma_T_list, yerr=sigma_T_err_list, fmt='s', color='red')
axs4[1].set_title(r'$\sigma_T$ vs $t$', fontsize=16, fontweight='bold')
axs4[1].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs4[1].set_ylabel(r'$\sigma_T$ ($\mu$b/GeV$^2$)', fontsize=14)

# σ_LT vs t
axs4[2].errorbar(t_plot_list, sigma_LT_list, yerr=sigma_LT_err_list, fmt='^', color='green')
axs4[2].set_title(r'$\sigma_{LT}$ vs $t$', fontsize=16, fontweight='bold')
axs4[2].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs4[2].set_ylabel(r'$\sigma_{LT}$ ($\mu$b/GeV$^2$)', fontsize=14)

# σ_TT vs t
axs4[3].errorbar(t_plot_list, sigma_TT_list, yerr=sigma_TT_err_list, fmt='v', color='purple')
axs4[3].set_title(r'$\sigma_{TT}$ vs $t$', fontsize=16, fontweight='bold')
axs4[3].set_xlabel(r'$t$ (GeV$^2$)', fontsize=14)
axs4[3].set_ylabel(r'$\sigma_{TT}$ ($\mu$b/GeV$^2$)', fontsize=14)

plt.tight_layout()

# --- Save both figures to the same PDF ---
with PdfPages(Sep_Xsection_pdf) as pdf:
    pdf.savefig(fig_2)
    pdf.savefig(fig4)
plt.close(fig_2)
plt.close(fig4)

print(f"Saved Sep cross-section plots to {Sep_Xsection_pdf}")

print ("Processing Complete")