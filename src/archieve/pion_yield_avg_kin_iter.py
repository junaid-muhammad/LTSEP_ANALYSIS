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

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=7:
    print("!!!!! ERROR !!!!!\n Expected 7 arguments\n Usage is with - CSV file suffix \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
DATA_YIELD_CSV = sys.argv[3]
DATA_AVG_KIN_CSV = sys.argv[4]
SIMC_YIELD_CSV = sys.argv[5]
SIMC_AVG_KIN_CSV = sys.argv[6]
Iteration = sys.argv[7]

################################################################################################################################################

# Naming Scheme for the creating input file for ltsep analysis

if PHY_SETTING == "Q3p85_W2p62_t0p21":
    loweps = "0.292"
    higheps = "0.779"
    pol = +1
    theta_c = +0.000
else:
    print("!!!!! Please declare low/mid/high epsilon values of you physics setting !!!!!")

# Section fo theta_cm calculation
# Define constants (masses in GeV, pi)
mp = 0.93827231    # Proton mass
mp2 = 0.88035493   # Proton mass squared
mpi = 0.13956995   # Charged pion mass
mpi2 = 0.01947977  # Charged pion mass squared
mn = 0.93956563    # Neutron mass
mn2 = 0.88278357   # Neutron mass squared
pi = 3.14159265    # Pi

# Assume m3 and m32 for pion by default
m3 = mpi
m32 = mpi2

# Set npol: e.g., npol = +1 for proton, npol = -1 for neutron
# npol = ... (set this value before using the logic below)

if pol > 0:
    m2 = mp
    m22 = mp2
    m4 = mn
    m42 = mn2
else:
    m2 = mn
    m22 = mn2
    m4 = mp
    m42 = mp2

##################################################################################################################################################

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
VOLATILEPATH=lt.VOLATILEPATH

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_iter%s" % (setting_name, Iteration)
DATA_CSV_PATH = "%s/LTSEP_ANALYSIS/LTSep_CSVs/physics_yields_csv/%s_std" % (REPLAYPATH, setting_name)
SIMC_CSV_PATH = "%s/LTSEP_ANALYSIS/LTSep_CSVs/simc_yields_csv/%s" % (REPLAYPATH, physet_dir_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

# Input file location and variables taking
data_yield_csv = "%s/%s_%s.csv" % (DATA_CSV_PATH, PHY_SETTING, DATA_YIELD_CSV)
data_avgkin_csv = "%s/%s_%s.csv" % (DATA_CSV_PATH, PHY_SETTING, DATA_AVG_KIN_CSV)
simc_yield_csv = "%s/%s_%s.csv" % (SIMC_CSV_PATH, PHY_SETTING, SIMC_YIELD_CSV)
simc_avgkin_csv = "%s/%s_%s.csv" % (SIMC_CSV_PATH, PHY_SETTING, SIMC_AVG_KIN_CSV)

###############################################################################################################################################

# Section for average yield calculation
print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("="*40)

# Define the output path
combined_weighted_loweps_yields = "%s/LTSEP_ANALYSIS/LTSep_CSVs/avg_kinematics_csv/%s/%s_pion_physics_loweps_avg_yields.csv" % (REPLAYPATH, physet_dir_name, PHY_SETTING)

# Load the CSV file
data_yield_df = pd.read_csv(data_yield_csv)
data_yield_df.columns = data_yield_df.columns.str.strip()  # clean up headers
simc_yield_df = pd.read_csv(simc_yield_csv)
simc_yield_df.columns = simc_yield_df.columns.str.strip()  # clean up headers

# Filter the DataFrame for specific Physics_Setting values
filtered_data_yield_loweps_center_df = data_yield_df[data_yield_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_center"]
filtered_data_yield_loweps_left_df = data_yield_df[data_yield_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_left"]
filtered_simc_yield_loweps_center_df = simc_yield_df[simc_yield_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_center"]
filtered_simc_yield_loweps_left_df = simc_yield_df[simc_yield_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_left"]

# Define binning keys
bin_yield_keys = ["tbin_number", "t_min", "t_max", "phibin_number", "phi_min", "phi_max"]
key_yield_cols = ["Physics_Setting"] + bin_yield_keys

# Merge filtered DataFrames on binning keys
merged_loweps_center_df = pd.merge(filtered_data_yield_loweps_center_df, filtered_simc_yield_loweps_center_df, on=key_yield_cols, suffixes=('_loweps_center_data', '_loweps_center_simc'))
merged_loweps_left_df = pd.merge(filtered_data_yield_loweps_left_df, filtered_simc_yield_loweps_left_df, on=key_yield_cols, suffixes=('_loweps_left_data', '_loweps_left_simc'))

# Calculate Physics_Yield ratio for low-ε center
merged_loweps_center_df["physics_yield"] = merged_loweps_center_df["physics_yield"].clip(lower=0)
merged_loweps_center_df["simc_yield"] = merged_loweps_center_df["simc_yield"].clip(lower=0)
merged_loweps_center_df["ratio_loweps_center"] = merged_loweps_center_df["physics_yield"] / merged_loweps_center_df["simc_yield"]
merged_loweps_center_df["ratio_loweps_center"] = merged_loweps_center_df["ratio_loweps_center"].replace(np.nan, 0)
merged_loweps_center_df["ratio_error_loweps_center"] = merged_loweps_center_df["ratio_loweps_center"] * np.sqrt(
    (merged_loweps_center_df["physics_yield_error"] / merged_loweps_center_df["physics_yield"].replace(0, np.nan))**2 +
    (merged_loweps_center_df["simc_yield_error"] / merged_loweps_center_df["simc_yield"].replace(0, np.nan))**2
)
merged_loweps_center_df["ratio_loweps_center"] = merged_loweps_center_df["ratio_loweps_center"].replace([np.inf, -np.inf], 0)
merged_loweps_center_df["ratio_error_loweps_center"] = merged_loweps_center_df["ratio_error_loweps_center"].replace(np.nan, 1000)
merged_loweps_center_df["ratio_error_loweps_center"] = np.abs(merged_loweps_center_df["ratio_error_loweps_center"])
# Print the results for low-ε center
#print("\n--- Low-ε Center Ratios ---")
#for _, row in merged_loweps_center_df.iterrows():
#    print(f"{row['Physics_Setting']} | t=({row['t_min']}, {row['t_max']}) | phi=({row['phi_min']}, {row['phi_max']}) | "
#          f"ratio = {row['ratio_loweps_center']:.4f} | error = {row['ratio_error_loweps_center']:.4f}")


# Calculate Physics_Yield ratio for low-ε left
merged_loweps_left_df["physics_yield"] = merged_loweps_left_df["physics_yield"].clip(lower=0)
merged_loweps_left_df["simc_yield"] = merged_loweps_left_df["simc_yield"].clip(lower=0)
merged_loweps_left_df["ratio_loweps_left"] = merged_loweps_left_df["physics_yield"] / merged_loweps_left_df["simc_yield"]
merged_loweps_left_df["ratio_loweps_left"] = merged_loweps_left_df["ratio_loweps_left"].replace(np.nan, 0)
merged_loweps_left_df["ratio_loweps_left"] = merged_loweps_left_df["ratio_loweps_left"].replace([np.inf, -np.inf], 0)
merged_loweps_left_df["ratio_error_loweps_left"] = merged_loweps_left_df["ratio_loweps_left"] * np.sqrt(
    (merged_loweps_left_df["physics_yield_error"] / merged_loweps_left_df["physics_yield"].replace(0, np.nan))**2 +
    (merged_loweps_left_df["simc_yield_error"] / merged_loweps_left_df["simc_yield"].replace(0, np.nan))**2
)
merged_loweps_left_df["ratio_error_loweps_left"] = merged_loweps_left_df["ratio_error_loweps_left"].replace(np.nan, 1000)
merged_loweps_left_df["ratio_error_loweps_left"] = np.abs(merged_loweps_left_df["ratio_error_loweps_left"])
# Print the results for low-ε left
#print("\n--- Low-ε Left Ratios ---")
#for _, row in merged_loweps_left_df.iterrows():
#    print(f"{row['Physics_Setting']} | t=({row['t_min']}, {row['t_max']}) | phi=({row['phi_min']}, {row['phi_max']}) | "
#          f"ratio = {row['ratio_loweps_left']:.4f} | error = {row['ratio_error_loweps_left']:.4f}")
    

# Merge center and left results on bin keys
merged_avg_loweps_df = pd.merge(merged_loweps_center_df, merged_loweps_left_df, on=bin_yield_keys, suffixes=('_loweps_center', '_loweps_left'))

# Drop duplicated Physics_Setting columns if they exist
for col in merged_avg_loweps_df.columns:
    if "Physics_Setting" in col:merged_avg_loweps_df.drop(columns=col, inplace=True)
# Add unified Physics_Setting label
merged_avg_loweps_df.insert(0, "Physics_Setting", f"{PHY_SETTING}_loweps")

# Perform error-weighted average for merged center and left
def error_weighted_ratio_loweps(row):
    y1, e1 = row["ratio_loweps_center"], row["ratio_error_loweps_center"]
    y2, e2 = row["ratio_loweps_left"], row["ratio_error_loweps_left"]
    
    # Handle zero or very small errors to avoid division by zero
    min_error = 1e-10  # Minimum error threshold
    e1 = max(e1, min_error)
    e2 = max(e2, min_error)
    
    # Compute weights (inverse variance)
    weights = 1 / np.array([e1**2, e2**2])
    ratios = np.array([y1, y2])
    weighted_data_ratio = np.sum(ratios * weights) / np.sum(weights)
    weighted_data_ratio_err = np.sqrt(1 / np.sum(weights))
    # Weighted average and propagated error
    return pd.Series({
        "avg_ratio": weighted_data_ratio,
        "avg_ratio_error": weighted_data_ratio_err
    })

# Apply row-wise
avg_loweps_results = merged_avg_loweps_df.apply(error_weighted_ratio_loweps, axis=1)
# Concatenate results
merged_avg_loweps_df = pd.concat([merged_avg_loweps_df, avg_loweps_results], axis=1)

# Select final output columns
final_ratio_loweps_df = merged_avg_loweps_df[[
    "Physics_Setting"] + bin_yield_keys +
    ["avg_ratio", "avg_ratio_error"
]]

# Save the final DataFrame to CSV
final_ratio_loweps_df.to_csv(combined_weighted_loweps_yields, index=False)

print(f"Average weighted loweps yields saved to the path: {combined_weighted_loweps_yields}")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define the output path
combined_weighted_higheps_yields = "%s/LTSEP_ANALYSIS/LTSep_CSVs/avg_kinematics_csv/%s/%s_pion_physics_higheps_avg_yields.csv" % (REPLAYPATH, physet_dir_name, PHY_SETTING)

# Filter the DataFrame for specific Physics_Setting values
filtered_data_yield_higheps_right_df = data_yield_df[data_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_right"]
filtered_data_yield_higheps_center_df = data_yield_df[data_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_center"]
filtered_data_yield_higheps_left_df = data_yield_df[data_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_left"]
filtered_simc_yield_higheps_right_df = simc_yield_df[simc_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_right"]
filtered_simc_yield_higheps_center_df = simc_yield_df[simc_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_center"]
filtered_simc_yield_higheps_left_df = simc_yield_df[simc_yield_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_left"]

# Merge filtered DataFrames on binning keys
merged_higheps_right_df = pd.merge(filtered_data_yield_higheps_right_df, filtered_simc_yield_higheps_right_df, on=key_yield_cols, suffixes=('_higheps_right_data', '_higheps_right_simc'))
merged_higheps_center_df = pd.merge(filtered_data_yield_higheps_center_df, filtered_simc_yield_higheps_center_df, on=key_yield_cols, suffixes=('_higheps_center_data', '_higheps_center_simc'))
merged_higheps_left_df = pd.merge(filtered_data_yield_higheps_left_df, filtered_simc_yield_higheps_left_df, on=key_yield_cols, suffixes=('_higheps_left_data', '_higheps_left_simc'))

# Calculate Physics_Yield ratio for high-ε center
# RIGHT
merged_higheps_right_df["physics_yield"] = merged_higheps_right_df["physics_yield"].clip(lower=0)
merged_higheps_right_df["simc_yield"] = merged_higheps_right_df["simc_yield"].clip(lower=0)
merged_higheps_right_df["ratio_higheps_right"] = merged_higheps_right_df["physics_yield"] / merged_higheps_right_df["simc_yield"]
merged_higheps_right_df["ratio_higheps_right"] = merged_higheps_right_df["ratio_higheps_right"].replace(np.nan, 0)
merged_higheps_right_df["ratio_higheps_right"] = merged_higheps_right_df["ratio_higheps_right"].replace([np.inf, -np.inf], 0)
merged_higheps_right_df["ratio_error_higheps_right"] = merged_higheps_right_df["ratio_higheps_right"] * np.sqrt(
    (merged_higheps_right_df["physics_yield_error"] / merged_higheps_right_df["physics_yield"].replace(0, np.nan))**2 +
    (merged_higheps_right_df["simc_yield_error"] / merged_higheps_right_df["simc_yield"].replace(0, np.nan))**2
)
merged_higheps_right_df["ratio_error_higheps_right"] = merged_higheps_right_df["ratio_error_higheps_right"].replace(np.nan, 1000)
merged_higheps_right_df["ratio_error_higheps_right"] = np.abs(merged_higheps_right_df["ratio_error_higheps_right"])

# CENTER
merged_higheps_center_df["physics_yield"] = merged_higheps_center_df["physics_yield"].clip(lower=0)
merged_higheps_center_df["simc_yield"] = merged_higheps_center_df["simc_yield"].clip(lower=0)
merged_higheps_center_df["ratio_higheps_center"] = merged_higheps_center_df["physics_yield"] / merged_higheps_center_df["simc_yield"]
merged_higheps_center_df["ratio_higheps_center"] = merged_higheps_center_df["ratio_higheps_center"].replace(np.nan, 0)
merged_higheps_center_df["ratio_higheps_center"] = merged_higheps_center_df["ratio_higheps_center"].replace([np.inf, -np.inf], 0)
merged_higheps_center_df["ratio_error_higheps_center"] = merged_higheps_center_df["ratio_higheps_center"] * np.sqrt(
    (merged_higheps_center_df["physics_yield_error"] / merged_higheps_center_df["physics_yield"].replace(0, np.nan))**2 +
    (merged_higheps_center_df["simc_yield_error"] / merged_higheps_center_df["simc_yield"].replace(0, np.nan))**2
)
merged_higheps_center_df["ratio_error_higheps_center"] = merged_higheps_center_df["ratio_error_higheps_center"].replace(np.nan, 1000)
merged_higheps_center_df["ratio_error_higheps_center"] = np.abs(merged_higheps_center_df["ratio_error_higheps_center"])

# LEFT
merged_higheps_left_df["physics_yield"] = merged_higheps_left_df["physics_yield"].clip(lower=0)
merged_higheps_left_df["simc_yield"] = merged_higheps_left_df["simc_yield"].clip(lower=0)
merged_higheps_left_df["ratio_higheps_left"] = merged_higheps_left_df["physics_yield"] / merged_higheps_left_df["simc_yield"]
merged_higheps_left_df["ratio_higheps_left"] = merged_higheps_left_df["ratio_higheps_left"].replace(np.nan, 0)
merged_higheps_left_df["ratio_higheps_left"] = merged_higheps_left_df["ratio_higheps_left"].replace([np.inf, -np.inf], 0)
merged_higheps_left_df["ratio_error_higheps_left"] = merged_higheps_left_df["ratio_higheps_left"] * np.sqrt(
    (merged_higheps_left_df["physics_yield_error"] / merged_higheps_left_df["physics_yield"].replace(0, np.nan))**2 +
    (merged_higheps_left_df["simc_yield_error"] / merged_higheps_left_df["simc_yield"].replace(0, np.nan))**2
)
merged_higheps_left_df["ratio_error_higheps_left"] = merged_higheps_left_df["ratio_error_higheps_left"].replace(np.nan, 1000)
merged_higheps_left_df["ratio_error_higheps_left"] = np.abs(merged_higheps_left_df["ratio_error_higheps_left"])

# Merge all three (right, center, left) on bin keys only (t, phi bins)
merged_temp_higheps_df = pd.merge(merged_higheps_right_df, merged_higheps_center_df, on=bin_yield_keys, suffixes=('_higheps_right', '_higheps_center'))
merged_avg_higheps_df = pd.merge(merged_temp_higheps_df, merged_higheps_left_df, on=bin_yield_keys, suffixes=('', '_higheps_left'))

# Drop duplicated Physics_Setting columns if they exist
for col in merged_avg_higheps_df.columns:
    if "Physics_Setting" in col:
        merged_avg_higheps_df.drop(columns=col, inplace=True)
merged_avg_higheps_df.insert(0, "Physics_Setting", f"{PHY_SETTING}_higheps")

# Perform error-weighted average for merged center and left
def error_weighted_ratio_higheps(row):
    y1, e1 = row["ratio_higheps_right"], row["ratio_error_higheps_right"]
    y2, e2 = row["ratio_higheps_center"], row["ratio_error_higheps_center"]
    y3, e3 = row["ratio_higheps_left"], row["ratio_error_higheps_left"]
    
    # Handle zero or very small errors to avoid division by zero
    min_error = 1e-10  # Minimum error threshold
    e1 = max(e1, min_error)
    e2 = max(e2, min_error)
    e3 = max(e3, min_error)
    
    # Compute weights (inverse variance)
    weights = 1 / np.array([e1**2, e2**2, e3**2])
    ratios = np.array([y1, y2, y3])
    weighted_data_ratio = np.sum(ratios * weights) / np.sum(weights)
    weighted_data_ratio_err = np.sqrt(1 / np.sum(weights))
    # Weighted average and propagated error
    return pd.Series({
        "avg_ratio": weighted_data_ratio,
        "avg_ratio_error": weighted_data_ratio_err
    })

# Apply row-wise
avg_higheps_results = merged_avg_higheps_df.apply(error_weighted_ratio_higheps, axis=1)
# Concatenate results
merged_avg_higheps_df = pd.concat([merged_avg_higheps_df, avg_higheps_results], axis=1)

# Select final output columns
final_ratio_higheps_df = merged_avg_higheps_df[
    ["Physics_Setting"] + bin_yield_keys +
    ["avg_ratio", "avg_ratio_error"
]]

# Save the final DataFrame to CSV
final_ratio_higheps_df.to_csv(combined_weighted_higheps_yields, index=False)

print(f"Average weighted higheps yields saved to the path: {combined_weighted_higheps_yields}")
print("-"*40)

#===========================================================================================================================================================================================================

# Section for average kinematic calculation
# Define the output path
combined_weighted_kinematics = "%s/LTSEP_ANALYSIS/LTSep_CSVs/avg_kinematics_csv/%s/%s_pion_physics_avg_kinematics.csv" % (REPLAYPATH, physet_dir_name, PHY_SETTING)

# Load the CSV file
data_avgkin_df = pd.read_csv(data_avgkin_csv)
data_avgkin_df.columns = data_avgkin_df.columns.str.strip()  # clean up headers
simc_avgkin_df = pd.read_csv(simc_avgkin_csv)
simc_avgkin_df.columns = simc_avgkin_df.columns.str.strip()  # clean up headers

# Filter the DataFrame for specific Physics_Setting values from Data
filtered_data_avgkin_loweps_center_df = data_avgkin_df[data_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_center"]
filtered_data_avgkin_loweps_left_df = data_avgkin_df[data_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_left"]
filtered_data_avgkin_higheps_right_df = data_avgkin_df[data_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_right"]
filtered_data_avgkin_higheps_center_df = data_avgkin_df[data_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_center"]
filtered_data_avgkin_higheps_left_df = data_avgkin_df[data_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_left"]
# Filter the DataFrame for specific Physics_Setting values from SimC
filtered_simc_avgkin_loweps_center_df = simc_avgkin_df[simc_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_center"]
filtered_simc_avgkin_loweps_left_df = simc_avgkin_df[simc_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_loweps_left"]
filtered_simc_avgkin_higheps_right_df = simc_avgkin_df[simc_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_right"]
filtered_simc_avgkin_higheps_center_df = simc_avgkin_df[simc_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_center"]
filtered_simc_avgkin_higheps_left_df = simc_avgkin_df[simc_avgkin_df["Physics_Setting"] == f"{PHY_SETTING}_higheps_left"]

# Define binning keys
bin_avgkin_keys = ["total_tbins", "tbin_number", "t_min", "t_max", "t_central", "total_phibins", "phi_central"]

# Rename columns for each setting (except for bin keys)
def data_add_suffix(df, suffix):
    return df.rename(columns={col: f"{col}_{suffix}" for col in df.columns if col not in bin_avgkin_keys})
loweps_center  = data_add_suffix(filtered_data_avgkin_loweps_center_df,  "loweps_center")
loweps_left    = data_add_suffix(filtered_data_avgkin_loweps_left_df,    "loweps_left")
higheps_right  = data_add_suffix(filtered_data_avgkin_higheps_right_df,  "higheps_right")
higheps_center = data_add_suffix(filtered_data_avgkin_higheps_center_df, "higheps_center")
higheps_left   = data_add_suffix(filtered_data_avgkin_higheps_left_df,   "higheps_left")
# Merge all on bin keys
merged_data_all = loweps_center.merge(loweps_left, on=bin_avgkin_keys, how='outer') \
    .merge(higheps_right,  on=bin_avgkin_keys, how='outer') \
    .merge(higheps_center, on=bin_avgkin_keys, how='outer') \
    .merge(higheps_left,   on=bin_avgkin_keys, how='outer')
#print(list(merged_data_all.columns))

# Perform error-weighted average for merged data
def error_weighted_avgkin_data(row):
    vars = ["avg_Q2", "avg_W"]
    suffixes = ["loweps_center", "loweps_left", "higheps_right", "higheps_center", "higheps_left"]
    result = {}
    for var in vars:
        values = [row[f"{var}_{sfx}"] for sfx in suffixes]
        errors = [row[f"{var}_err_{sfx}"] for sfx in suffixes]
        
        # Handle zero or very small errors to avoid division by zero
        errors = np.array(errors)
        min_error = 1e-10  # Minimum error threshold
        errors = np.where(errors <= min_error, min_error, errors)
        
        weights = 1 / np.square(errors)
        weighted_avg = np.sum(np.array(values) * weights) / np.sum(weights)
        weighted_avg_err = np.sqrt(1 / np.sum(weights))
        result[f"data_{var}"] = weighted_avg
        result[f"data_{var}_err"] = weighted_avg_err
    # theta_cm Calculation
    Eb_loweps = row["Beam_Energy_loweps_center"]
    Eb_higheps = row["Beam_Energy_higheps_center"]
    tm = row["t_central"]
    q2 = float(result["data_avg_Q2"])
    w = float(result["data_avg_W"])
    s = w * w
    m12 = -q2
    omega = (s + q2 - m22) / (2 * m2)
    q = math.sqrt(q2 + omega**2)
    e1cm = (s + m12 - m22) / (2 * w)
    e3cm = (s + m32 - m42) / (2 * w)
    p1lab = q
    p1cm = p1lab * m2 / w
    p3cm = math.sqrt(e3cm * e3cm - m32)
    tmin = -((e1cm - e3cm)**2 - (p1cm - p3cm)**2)
    thetacm = 2 * math.asin(math.sqrt((tm - tmin) / (4 * p1cm * p3cm)))
    thetacm_deg = thetacm * 180 / math.pi  # Convert to degrees
    low_eps = 1/(1 + 2 * (q2 + omega**2)/(4 * Eb_loweps * (Eb_loweps-omega) - q2))
    high_eps = 1/(1 + 2 * (q2 + omega**2)/(4 * Eb_higheps * (Eb_higheps-omega) - q2))
    result["theta_cm"] = thetacm_deg
    result["epsilon_low"] = low_eps
    result["epsilon_high"] = high_eps
    result["Beam_Energy_loweps"] = row["Beam_Energy_loweps_center"]
    result["Beam_Energy_higheps"] = row["Beam_Energy_higheps_center"]
    return pd.Series(result)

# Rename columns for each setting (except for bin keys)
def simc_add_suffix(df, suffix):
    return df.rename(columns={col: f"{col}_{suffix}" for col in df.columns if col not in bin_avgkin_keys})
loweps_center  = simc_add_suffix(filtered_simc_avgkin_loweps_center_df,  "loweps_center")
loweps_left    = simc_add_suffix(filtered_simc_avgkin_loweps_left_df,    "loweps_left")
higheps_right  = simc_add_suffix(filtered_simc_avgkin_higheps_right_df,  "higheps_right")
higheps_center = simc_add_suffix(filtered_simc_avgkin_higheps_center_df, "higheps_center")
higheps_left   = simc_add_suffix(filtered_simc_avgkin_higheps_left_df,   "higheps_left")
# Merge all on bin keys
merged_simc_all = loweps_center.merge(loweps_left, on=bin_avgkin_keys, how='outer') \
    .merge(higheps_right,  on=bin_avgkin_keys, how='outer') \
    .merge(higheps_center, on=bin_avgkin_keys, how='outer') \
    .merge(higheps_left,   on=bin_avgkin_keys, how='outer')
#print(list(merged_data_all.columns))

def error_weighted_avgkin_simc(row):
    vars = ["avg_Q2", "avg_W"]
    suffixes = ["loweps_center", "loweps_left", "higheps_right", "higheps_center", "higheps_left"]
    result = {}
    for var in vars:
        values = [row[f"{var}_{sfx}"] for sfx in suffixes]
        errors = [row[f"{var}_err_{sfx}"] for sfx in suffixes]
        
        # Handle zero or very small errors to avoid division by zero
        errors = np.array(errors)
        min_error = 1e-10  # Minimum error threshold
        errors = np.where(errors <= min_error, min_error, errors)
        
        weights = 1 / np.square(errors)
        weighted = np.sum(np.array(values) * weights) / np.sum(weights)
        weighted_err = np.sqrt(1 / np.sum(weights))
        result[f"simc_{var}"] = weighted
        result[f"simc_{var}_err"] = weighted_err
    return pd.Series(result)

# Assuming merged_all is your merged DataFrame
result_data_df = merged_data_all.apply(error_weighted_avgkin_data, axis=1)
result_simc_df = merged_simc_all.apply(error_weighted_avgkin_simc, axis=1)

# Combine all into a single DataFrame
final_data_df = pd.concat([merged_data_all[bin_avgkin_keys].reset_index(drop=True), result_data_df.reset_index(drop=True)], axis=1)
final_simc_df = pd.concat([merged_simc_all[bin_avgkin_keys].reset_index(drop=True), result_simc_df.reset_index(drop=True)], axis=1)

# 8. Combine both (drop duplicate bin columns from simc first!)
final_simc_nobin = final_simc_df.drop(columns=bin_avgkin_keys)
final_combined_df = pd.concat([final_data_df, final_simc_nobin], axis=1)

# Add Physics_Setting column and compute ratio
final_combined_df.insert(0, "Physics_Setting", f"{PHY_SETTING}")

# Reorder columns for final output
final_avgkin_df = final_combined_df[[
    "Physics_Setting", "Beam_Energy_loweps", "Beam_Energy_higheps", "epsilon_low", "epsilon_high", "total_tbins",
    "tbin_number", "t_min", "t_max", "t_central", "total_phibins", "phi_central", "theta_cm",
    "data_avg_Q2", "data_avg_Q2_err",
    "data_avg_W", "data_avg_W_err",
    "simc_avg_Q2", "simc_avg_Q2_err",
    "simc_avg_W", "simc_avg_W_err",
]]

final_avgkin_df.to_csv(combined_weighted_kinematics, index=False)

print(f"Average weighted loweps kinematics saved to the path: {combined_weighted_kinematics}")
print("-"*40)

#=========================================================================================================================================================================================================================================================
# Plotting Section fot loweps
loweps_df = pd.read_csv(combined_weighted_loweps_yields)
loweps_df["phi_avg"] = (loweps_df["phi_min"] + loweps_df["phi_max"]) / 2
loweps_tbins = loweps_df[["t_min", "t_max"]].drop_duplicates().reset_index(drop=True)

# Setup plot grid
fig_loweps, axs_loweps = plt.subplots(2, 3, figsize=(20, 10), sharey=False)  # sharey=False to show individual y-ticks
axs_loweps = axs_loweps.flatten()

# First subplot: physics setting text
axs_loweps[0].axis('off')
axs_loweps[0].text(0.5, 0.5, f"Physics Setting:\n{PHY_SETTING}_loweps",
                   ha='center', va='center', fontsize=18, fontweight='bold',
                   transform=axs_loweps[0].transAxes)

# Plot the 5 t-bins
for t_idx in range(min(5, len(loweps_tbins))):
    t_min = loweps_tbins.loc[t_idx, "t_min"]
    t_max = loweps_tbins.loc[t_idx, "t_max"]
    t_df = loweps_df[(loweps_df["t_min"] == t_min) & (loweps_df["t_max"] == t_max)]

    ax = axs_loweps[t_idx + 1]
    ax.errorbar(t_df["phi_avg"], t_df["avg_ratio"], yerr=t_df["avg_ratio_error"],
                fmt='o', color='blue', capsize=4, label='Data/SIMC')
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2)
    ax.set_title(f"$t$: [{t_min:.2f} - {t_max:.2f}]")
    ax.set_xlabel("ϕ (deg)")
    ax.set_xlim(0, 360)
#    ax.set_ylim(0.0, max(1.5, t_df["avg_ratio"].max() + 0.4))
    ax.set_ylim(0.5, 1.5)
    ax.set_ylabel("Avg Data/SIMC Ratio")
    ax.tick_params(axis='y', which='both', labelleft=True)  # ensure y-tick labels are shown

# Adjust layout and save
plt.tight_layout()

# Plotting Section fot higheps
higheps_df = pd.read_csv(combined_weighted_higheps_yields)
higheps_df["phi_avg"] = (higheps_df["phi_min"] + higheps_df["phi_max"]) / 2
higheps_tbins = higheps_df[["t_min", "t_max"]].drop_duplicates().reset_index(drop=True)

# Setup plot grid
fig_higheps, axs_higheps = plt.subplots(2, 3, figsize=(20, 10), sharey=False)  # sharey=False to show individual y-ticks
axs_higheps = axs_higheps.flatten()

# First subplot: physics setting text
axs_higheps[0].axis('off')
axs_higheps[0].text(0.5, 0.5, f"Physics Setting:\n{PHY_SETTING}_higheps",
                   ha='center', va='center', fontsize=18, fontweight='bold',
                   transform=axs_higheps[0].transAxes)

# Plot the 5 t-bins
for t_idx in range(min(5, len(higheps_tbins))):
    t_min = higheps_tbins.loc[t_idx, "t_min"]
    t_max = higheps_tbins.loc[t_idx, "t_max"]
    t_df = higheps_df[(higheps_df["t_min"] == t_min) & (higheps_df["t_max"] == t_max)]

    ax = axs_higheps[t_idx + 1]
    ax.errorbar(t_df["phi_avg"], t_df["avg_ratio"], yerr=t_df["avg_ratio_error"],
                fmt='o', color='black', capsize=4, label='Data/SIMC')
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2)
    ax.set_title(f"$t$: [{t_min:.2f} - {t_max:.2f}]")
    ax.set_xlabel("ϕ (deg)")
    ax.set_xlim(0, 360)
#    ax.set_ylim(0.0, max(1.5, t_df["avg_ratio"].max() + 0.4))
    ax.set_ylim(0.5, 1.5)
    ax.set_ylabel("Avg Data/SIMC Ratio")
    ax.tick_params(axis='y', which='both', labelleft=True)  # ensure y-tick labels are shown

# Adjust layout and save
plt.tight_layout()

ratio_pdf = "%s/LTSEP_ANALYSIS/src/plots/%s_%s_ProdCoin_Pion_Analysis_avgkin_avgyield_Distributions_iter%s.pdf" % (REPLAYPATH, PHY_SETTING, MaxEvent, Iteration)
# Save both plots into one PDF
with PdfPages(ratio_pdf) as pdf:
    # Save loweps figure
    pdf.savefig(fig_loweps)
    plt.close(fig_loweps)  # Close after saving to free memory
    # Save higheps figure
    pdf.savefig(fig_higheps)
    plt.close(fig_higheps)  # Close after saving

#==========================================================================================================================================================================================================================================================

# Section for creating input files for cross-section and LTSep analysis
# Creating naming scheme for input files
q_name_str = PHY_SETTING.split('_')[0]
q_name_val = int(q_name_str[1:].replace('p', ''))  # Remove 'Q', replace 'p' with '', then convert to int
loweps_val = int(loweps.replace('.', '')[1:])
higheps_val = int(higheps.replace('.', '')[1:]) 
#print(f"Q value: {q_name_val}, Loweps value: {loweps_val}, Higheps value: {higheps_val}")

# Input file for loweps data
avgkin_dat = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_avek_%s.dat" % (REPLAYPATH, physet_dir_name, PHY_SETTING, q_name_val)
with open(avgkin_dat, "w") as f:
    for idx, row in final_avgkin_df.iterrows():
        # Build the line as: Q2   (3 spaces) Q2err   (3 spaces) W   (3) Werr   (3) theta   (3) tbin_number
        line = (
            f"{row['data_avg_W']:<10.6f}"     # W
            + " " * 3 +
            f"{row['data_avg_W_err']:<10.6f}"  # W error
            + " " * 3 +
            f"{row['data_avg_Q2']:<10.6f}"   # Q2
            + " " * 3 +
            f"{row['data_avg_Q2_err']:<10.6f}"  # Q2 error
            + " " * 3 +
            f"{row['theta_cm']:<10.6f}"  # theta
            + " " * 3 +
            f"{int(row['tbin_number'])}\n"    # tbin_number as integer (no decimals)
        )
        f.write(line)
print(f"Saved input avgkin_dat file to: {avgkin_dat}")
print("-"*40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# input file with ratios for loweps
avgratio_loweps_dat = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_aver_pi_%s_%s.dat" % (REPLAYPATH, physet_dir_name, PHY_SETTING, q_name_val, loweps_val)
# Read the CSV file into a DataFrame
with open(avgratio_loweps_dat, "w") as f:
    for idx, row in final_ratio_loweps_df.iterrows():
        line = (
            f"{row['avg_ratio']:<10.6f}"
            + " " * 3 +
            f"{row['avg_ratio_error']:<10.6f}"
            + " " * 3 +
            f"{int(row['tbin_number'])}"
            + " " * 3 +
            f"{int(row['phibin_number'])}\n"
        )
        f.write(line)
print(f"Saved input avgkin_loweps_dat file to: {avgratio_loweps_dat}")


# input file with ratios for higheps
avgratio_higheps_dat = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_aver_pi_%s_%s.dat" % (REPLAYPATH, physet_dir_name, PHY_SETTING, q_name_val, higheps_val)
# Read the CSV file into a DataFrame
with open(avgratio_higheps_dat, "w") as f:
    for idx, row in final_ratio_higheps_df.iterrows():
        line = (
            f"{row['avg_ratio']:<10.6f}"
            + " " * 3 +
            f"{row['avg_ratio_error']:<10.6f}"
            + " " * 3 +
            f"{int(row['tbin_number'])}"
            + " " * 3 +
            f"{int(row['phibin_number'])}\n"
        )
        f.write(line)
print(f"Saved input avgkin_higheps_dat file to: {avgratio_higheps_dat}")

print("-"*40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create input file for tbin interval
tbin_interval_dat = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_t_bin_interval" % (REPLAYPATH, physet_dir_name, PHY_SETTING)
# Simple extraction and conversion of Q value
q_str = PHY_SETTING.split('_')[0]  # 'Q3p85'
q_val = float(q_str[1:].replace('p', '.'))  # Remove 'Q', replace 'p' with '', then convert to int
#print(q_val)  # This will print 3.85 for Q3p85
# Extract total_tbins and total_phibins from the first row of the DataFrame
first_row_nbins = final_avgkin_df.iloc[0]
total_tbins = int(first_row_nbins["total_tbins"])
total_phibins = int(first_row_nbins["total_phibins"])
with open(tbin_interval_dat, "w") as f:
    # First line: Q2, total_tbins, total_phibins
    f.write(f"{q_val:<10.3f}" + " " * 2 + f"{total_tbins}" + " " * 3 + f"{total_phibins}\n")
    # Second line: all tbin edges (tmin and tmax, sorted, unique, 3 decimal places)
    t_edges = np.unique(np.concatenate([final_avgkin_df["t_min"].values, final_avgkin_df["t_max"].values]))
    t_edges = np.sort(t_edges)
    for t_edge in t_edges:
        f.write(f"{t_edge:.3f}" + " " * 3)
    f.write("\n")  # Newline at the end of the t edges
print(f"Saved input t_bin_interval.dat file to: {tbin_interval_dat}")

print("-"*40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create input file for Ebeam file
Eb_pion_dat = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_Eb_pion.dat" % (REPLAYPATH, physet_dir_name, PHY_SETTING)
# Extract beam energy values
loweps_beam_energy = final_avgkin_df.iloc[0]["Beam_Energy_loweps"] * 1000
higheps_beam_energy = final_avgkin_df.iloc[0]["Beam_Energy_higheps"] * 1000
with open(Eb_pion_dat, "w") as f:
    # First row: loweps beam energy, q_val, loweps_val
    f.write(f"{loweps_beam_energy:<10.4f}" + " " * 3 + f"{q_val:<10.3f}" + " " * 2 + f"{loweps}\n")
    f.write(f"{higheps_beam_energy:<10.4f}" + " " * 3 + f"{q_val:<10.3f}" + " " * 2 + f"{higheps}\n")
print(f"Saved input Eb_pion22_dat file to: {Eb_pion_dat}")


# Creating list of settings for LTSep analysis
list_settings_pion = "%s/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/%s/%s_list_settings_pion" % (REPLAYPATH, physet_dir_name, PHY_SETTING)
with open(list_settings_pion, "w") as f:
    tmn = final_avgkin_df["t_min"].min()
    tmx = final_avgkin_df["t_max"].max()
    nbin = 0
    # List of (eps value, description) for looping
    epsilons = [loweps, higheps]
    for idx, eps in enumerate(epsilons, 1):
        line = (
            f"{int(pol):+3d}" + " " * 3 +
            f"{float(q_val):<10.3f}" + " " * 3 +
            f"{float(eps):<10.3f}" + " " * 2 +
            f"{float(theta_c):+10.3f}" + " " * 3 +
            f"{float(tmn):<10.3f}" + " " * 3 +
            f"{float(tmx):<10.3f}" + " " * 3 +
            f"{int(nbin)}" + " " * 3 +
            f"{idx}"
        )
        f.write(line.rstrip() + "\n")  # .rstrip() removes trailing spaces on each line
print(f"Saved input list_settings_pion22 file to: {list_settings_pion}")
print("-"*40)

print ("Processing Complete")