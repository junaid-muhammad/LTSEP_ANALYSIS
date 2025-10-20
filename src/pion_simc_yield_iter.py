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

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=4:
    print("!!!!! ERROR !!!!!\n Expected 7 arguments\n Usage is with - ROOTfileSuffixs Beam Energy MaxEvents RunList CVSFile\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
SIMC_Suffix = sys.argv[3]
Iteration = sys.argv[4]

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
VOLATILEPATH=lt.VOLATILEPATH

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_iter%s" % (setting_name, Iteration)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/OUTPUT/Analysis/SIMC/%s/" % (USER, physet_dir_name)
MMCUT_CSV   = "%s/LTSEP_ANALYSIS/LTSep_CSVs/mm_offset_cut_csv/%s_std" % (REPLAYPATH, setting_name)
DCUT_CSV    = "%s/LTSEP_ANALYSIS/LTSep_CSVs/diamond_cut_csv/%s_std" % (REPLAYPATH, setting_name)
TBINCSVPATH = "%s/LTSEP_ANALYSIS/LTSep_CSVs/t_binning_csv/%s_std" % (REPLAYPATH, setting_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
#Pion_Analysis_Distributions = "%s/%s_%s_ProdCoin_Pion_Analysis_SIMCYield_Distributions.pdf" % (OUTPATH, PHY_SETTING, MaxEvent)

# Input file location and variables taking
rootFile_SIMC = "%s/%s.root" % (SIMCPATH, SIMC_Suffix)
mmcut_csv_file = "%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, setting_name)
dcut_csv_file = "%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, setting_name)
tbin_csv_file  = "%s/%s_tbinning_yields_pions.csv" % (TBINCSVPATH, setting_name)

###################################################################################################################################################

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

# SIMC Cuts for Pions Selection
HMS_Acceptance = lambda event: (event.hsdelta >= -8.0) & (event.hsdelta <= 8.0) & (event.hsxpfp >= -0.08) & (event.hsxpfp <= 0.08) & (event.hsypfp >= -0.045) & (event.hsypfp <= 0.045)
SHMS_Acceptance = lambda event: (event.ssdelta >= -10.0) & (event.ssdelta <= 20.0) & (event.ssxpfp >= -0.06) & (event.ssxpfp <= 0.06) & (event.ssypfp >= -0.04) & (event.ssypfp <= 0.04)

# Cuts for Pions Selection - Change according to your data
# Read the vertices from the CSV file
vertices = {}
try:
    with open(dcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            vertex_name = row["Vertex"]
            x_value = round(float(row["X"]), 3)  # Round to 3 decimal places
            y_value = round(float(row["Y"]), 3)  # Round to 3 decimal places
            vertices[vertex_name] = [x_value, y_value]
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error reading diamond cut vertices from {dcut_csv_file}: {e}")
    sys.exit(1)

# Assign the vertices
vertex1 = vertices["vertex1"]  # bottom-left
vertex2 = vertices["vertex2"]  # top-left
vertex3 = vertices["vertex3"]  # top-right
vertex4 = vertices["vertex4"]  # bottom-right

# Print the vertex values rounded to 3 decimals
#print("Diamond Cut Vertices (rounded to 3 decimals):")
#print(f"Vertex 1 (bottom-left): [{vertex1[0]:.3f}, {vertex1[1]:.3f}]")
#print(f"Vertex 2 (top-left): [{vertex2[0]:.3f}, {vertex2[1]:.3f}]")
#print(f"Vertex 3 (top-right): [{vertex3[0]:.3f}, {vertex3[1]:.3f}]")
#print(f"Vertex 4 (bottom-right): [{vertex4[0]:.3f}, {vertex4[1]:.3f}]")

# Define the diamond cut
cutg_diamond = ROOT.TCutG("cutg_diamond", 5)
cutg_diamond.SetVarX("Q2")
cutg_diamond.SetVarY("W")
cutg_diamond.SetPoint(0, vertex1[0], vertex1[1])  # bottom-left
cutg_diamond.SetPoint(1, vertex2[0], vertex2[1])  # top-left
cutg_diamond.SetPoint(2, vertex3[0], vertex3[1])  # top-right
cutg_diamond.SetPoint(3, vertex4[0], vertex4[1])  # bottom-right
cutg_diamond.SetPoint(4, vertex1[0], vertex1[1])  # bottom-left again to close the loop
Diamond_Cut = lambda event: (cutg_diamond.IsInside(event.Q2, event.W))

#------------------------------------------------------------------------------------------------------------------------------------

# Read the MMpi cut values from the CSV file
try:
    with open(mmcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        first_row = next(csv_reader)  # Get the first row
        if first_row:
            MM_Cut_lowvalue = float(first_row["MM_Cut_low"])  # Assign MM_Cut_low
            MM_Cut_highvalue = float(first_row["MM_Cut_high"])  # Assign MM_Cut_high
        else:
            raise ValueError(f"The CSV file {mmcut_csv_file} is empty or missing required columns.")
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
#print(f"MMpi_Offset = {MMpi_Offset:.6f}")
#print(f"MMpi_Cut_lowvalue = {MMpi_Cut_lowvalue}")
#print(f"MMpi_Cut_highvalue = {MMpi_Cut_highvalue}")
SIMC_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= event.missmass <= MM_Cut_highvalue)

#-----------------------------------------------------------------------------------------------------------------------------------------

# t-binning cut
# Read the CSV file
tbin_df = pd.read_csv(tbin_csv_file)

# Extract the t_min and t_max columns as arrays
t_min_values = tbin_df['t_min'].values
t_max_values = tbin_df['t_max'].values

# Define the cuts using row 1 and row 3 values
tbin_Cut1 = lambda event: (t_min_values[0] <= -event.t <= t_max_values[0])  # Row 1
tbin_Cut2 = lambda event: (t_min_values[1] <= -event.t <= t_max_values[1])  # Row 2
tbin_Cut3 = lambda event: (t_min_values[2] <= -event.t <= t_max_values[2])  # Row 3
tbin_Cut4 = lambda event: (t_min_values[3] <= -event.t <= t_max_values[3])  # Row 4
tbin_Cut5 = lambda event: (t_min_values[4] <= -event.t <= t_max_values[4])  # Row 5

# Bundle them into a list
tbin_cuts = [tbin_Cut1, tbin_Cut2, tbin_Cut3, tbin_Cut4, tbin_Cut5]

# Print the t-binning cuts with 3 decimal places
#print("\n t-binning cuts:")
#print(f"tbin_Cut1: t_min = {t_min_values[0]:.3f}, t_max = {t_max_values[0]:.3f}")
#print(f"tbin_Cut2: t_min = {t_min_values[1]:.3f}, t_max = {t_max_values[1]:.3f}")
#print(f"tbin_Cut3: t_min = {t_min_values[2]:.3f}, t_max = {t_max_values[2]:.3f}")
#print(f"tbin_Cut4: t_min = {t_min_values[3]:.3f}, t_max = {t_max_values[3]:.3f}")
#print(f"tbin_Cut5: t_min = {t_min_values[4]:.3f}, t_max = {t_max_values[4]:.3f}")
#print("-"*40)

# Define phi bins (15 bins from 0 to 360 degrees, each 24 degrees wide)
phi_bins = [i for i in range(0, 361, 24)]  # 0, 24, 48, ..., 360

# Calculate total number of t-bins and phi bins
total_tbins = len(tbin_cuts)
total_phibins = len(phi_bins) - 1
#print(f"Total t-bins: {total_tbins}, Total phi bins: {total_phibins}")

#########################################################################################################################################

# Grabs simc number of events and normalizaton factor
simc_hist = "%s/%s.hist" % (SIMCPATH, SIMC_Suffix)

f_simc = open(simc_hist)
for line in f_simc:
#    print(line)
    if "Ncontribute" in line:
        val = line.split("=")
        simc_nevents = int(val[1])
    if "normfac" in line:
        val = line.split("=")
        simc_normfactor = float(val[1])
if 'simc_nevents' and 'simc_normfactor' in locals():
    print('\nsimc_nevents = ',simc_nevents,'\nsimc_normfactor = ',simc_normfactor,'\n')
#    print('\n\ndata_charge = {:.4f} +/- {:.4f}'.format(data_charge, data_charge*eff_errProp_data),'\ndummy_charge = {:.4f} +/- {:.4f}'.format(dummy_charge, dummy_charge*eff_errProp_dummy),'\n\n')
else:
    print("ERROR: Invalid simc hist file %s" % simc_hist)
    sys.exit(1)
f_simc.close()

# Normalization factor Calculation for SIMC
normfac_simc = (simc_normfactor)/(simc_nevents)
print("-"*40)
print ("normfac_simc: ", normfac_simc)
print("-"*40)

# Define the output CSV file name
normfact_pions_path =  "%s/LTSEP_ANALYSIS/LTSep_CSVs/simc_yields_csv/%s/%s_Physics_SIMC_Normfact_Calculation.csv" % (REPLAYPATH, physet_dir_name, setting_name)

# Define CSV header
header = ["Physics_Setting", "simc_normfactor", "simc_nevents", "normfac_simc"]

# New row to be written
new_row = {
    "Physics_Setting": PHY_SETTING,
    "simc_normfactor": simc_normfactor,
    "simc_nevents": simc_nevents,
    "normfac_simc": normfac_simc
}

# Read existing CSV file
rows = []
try:
    with open(normfact_pions_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = list(csv_reader)
except FileNotFoundError:
    # If file doesn't exist yet
    rows = []

# Update row if Physics_Setting already exists
row_found = False
for row in rows:
    if row["Physics_Setting"] == PHY_SETTING:
        row.update(new_row)
        row_found = True
        break

# Append new row if not found
if not row_found:
    rows.append(new_row)

# Write to file
with open(normfact_pions_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(rows)

print(f"SIMC normalization info written to {normfact_pions_path}")

###################################################################################################################################################

nbins = 200

hist_config = {
    "MMpi":    ("Missing Mass Distribution (Cut_All); MM_{\pi}; Counts", nbins, 0.8, 1.2),
    "Q2":      ("Q^{2} Distribution (Cut_All); Q^{2}; Counts", nbins, 2.0, 5.0),
    "W":       ("W Distribution (Cut_All); W; Counts", nbins, 2.0, 3.0),
    "epsilon": ("Epsilon Distribution (Cut_All); #epsilon; Counts", nbins, 0.0, 1.0),
    "theta":   ("Theta Distribution (Cut_All); #theta; Counts", nbins, 0.0, 0.5),
}

def initialize_histograms(t_min_values, t_max_values, phi_bins, hist_config, category_suffixes):
    histograms_by_tphi_cut = []
    added_t_bin_vars = set()  # Track which (tmin, tmax, suffix, var) already had W/Q2 added

    for i in range(len(t_min_values)):
        tmin = t_min_values[i]
        tmax = t_max_values[i]
        t_label = f"t({tmin:.2f}-{tmax:.2f})".replace('.', 'p')
        for j in range(len(phi_bins) - 1):
            phimin = phi_bins[j]
            phimax = phi_bins[j + 1]
            phi_label = f"phi({phimin}-{phimax})"
            hist_set = {suffix: {} for suffix in category_suffixes}

            for var, (base_title, nbins, xmin, xmax) in hist_config.items():
                for suffix in category_suffixes:
                    if var == "MMpi":
                        # MMpi → t and φ bins
                        hname = f"{var}_{t_label}_{phi_label}_{suffix}_data_cut_all"
                        title = f"{base_title} (t ∈ [{tmin:.2f}, {tmax:.2f}], φ ∈ [{phimin}, {phimax}], Category: {suffix})"
                        hist_set[suffix][var] = ROOT.TH1D(hname, title, nbins, xmin, xmax)
                    elif var in ("Q2", "W", "epsilon", "theta") and (tmin, tmax, suffix, var) not in added_t_bin_vars:
                        # Q2, W → only once per t-bin across full phi
                        hname = f"{var}_{t_label}_{suffix}_data_cut_all"
                        title = f"{base_title} (t ∈ [{tmin:.2f}, {tmax:.2f}], φ ∈ [0, 360], Category: {suffix})"
                        hist_set[suffix][var] = ROOT.TH1D(hname, title, nbins, xmin, xmax)
                        added_t_bin_vars.add((tmin, tmax, suffix, var))

            histograms_by_tphi_cut.append((tmin, tmax, phimin, phimax, hist_set))
    return histograms_by_tphi_cut

categories_simc = ["simc", "norm_simc", "simc_raw"]

# Initialize histograms for data
data_histograms_by_tphi_cut = initialize_histograms(t_min_values, t_max_values, phi_bins, hist_config, categories_simc)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_SIMC = ROOT.TFile.Open(rootFile_SIMC, "READ")
Uncut_Pion_Events_SIMC_tree = infile_SIMC.Get("h10")
#nEntries_TBRANCH_SIMC  = Uncut_Proton_Events_SIMC_tree.GetEntries()

###################################################################################################################################################

# Fill the histograms with the data
for event in Uncut_Pion_Events_SIMC_tree:
    if HMS_Acceptance(event) and SHMS_Acceptance(event) and SIMC_MMpi_Cut(event) and Diamond_Cut(event):
        for tbin_cut, (tmin, tmax) in zip(tbin_cuts, zip(t_min_values, t_max_values)):
            if tbin_cut(event):
                # Fill Q² and W (once per t-bin)
                for tmin_hist, tmax_hist, _, _, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist:
                        if "Q2" in hist_set["simc"]:
                            hist_set["simc"]["Q2"].Fill(event.Q2, event.Weight)
                        if "W" in hist_set["simc"]:
                            hist_set["simc"]["W"].Fill(event.W, event.Weight)
                        if "epsilon" in hist_set["simc"]:
                            hist_set["simc"]["epsilon"].Fill(event.epsilon, event.Weight)
                        if "theta" in hist_set["simc"]:
                            hist_set["simc"]["theta"].Fill(event.thetapq, event.Weight)
                        break  # Only fill Q2/W once for this t-bin
                # Fill MMπ per t–φ bin
                phi_deg = event.phipq * (180 / math.pi) + 180  # Convert φ to degrees in [0, 360]
                for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist and phimin_hist <= phi_deg <= phimax_hist:
                        hist_set["simc"]["MMpi"].Fill(event.missmass, event.Weight)
                        hist_set["simc_raw"]["MMpi"].Fill(event.missmass)
                        break  # Only fill in the matching t–φ bin

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

############################################################################################################################################

# Scale randsub_data and randsub_dummy histograms with normalization factors
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        # Scale the simc histogram
        if "simc" in hist_set and var in hist_set["simc"]:
            simc_hist = hist_set["simc"][var]
            if simc_hist:
                norm_simc_hist = simc_hist.Clone()
                norm_simc_hist.Scale(normfac_simc)
                hist_set["norm_simc"][var] = norm_simc_hist

print("####################################")
print("###### SIMC Normalization done ######")
print("####################################\n")

#############################################################################################################################################

print(f"Average Q2 , W, phi, eps and theta values for each t-bin (with uncertainties)")

avg_kin_tbin = {}
for tmin, tmax in zip(t_min_values, t_max_values):
    for tmin_hist, tmax_hist, phimin_hist, _, hist_set in data_histograms_by_tphi_cut:
        if tmin == tmin_hist and tmax == tmax_hist and phimin_hist == 0:
            hist_Q = hist_set["norm_simc"]["Q2"]
            hist_W = hist_set["norm_simc"]["W"]
            hist_eps = hist_set["norm_simc"]["epsilon"]
            hist_theta = hist_set["norm_simc"]["theta"]

            avg_Q2 = hist_Q.GetMean()
            err_Q2 = hist_Q.GetMeanError()
            avg_W = hist_W.GetMean()
            err_W = hist_W.GetMeanError()
            avg_eps = hist_eps.GetMean()
            err_eps = hist_eps.GetMeanError()
            avg_theta = hist_theta.GetMean()
            err_theta = hist_theta.GetMeanError()

            tbin_key = f"t({tmin:.2f}-{tmax:.2f})"
            avg_kin_tbin[tbin_key] = {
                "avg_Q2": avg_Q2, "err_Q2": err_Q2,
                "avg_W": avg_W, "err_W": err_W,
                "avg_eps": avg_eps, "err_eps": err_eps,
                "avg_theta": avg_theta, "err_theta": err_theta
            }

            print(f"{tbin_key}: <Q²> = {avg_Q2:.4f} ± {err_Q2:.4f}, <W> = {avg_W:.4f} ± {err_W:.4f}, <eps> = {avg_eps:.4f} ± {err_eps:.4f}, <θ> = {avg_theta:.4f} ± {err_theta:.4f}")
            break

# Define output CSV file path for average kinematics
avg_kinematics_path = "%s/LTSEP_ANALYSIS/LTSep_CSVs/simc_yields_csv/%s/%s_Physics_Avg_SIMC_Kinematics.csv" % (REPLAYPATH, physet_dir_name, setting_name)

# Define the header
avg_kin_header = [
    "Physics_Setting", "total_tbins", "tbin_number", "t_min", "t_max", "t_central",
    "total_phibins", "phi_central", "avg_Q2", "avg_Q2_err", "avg_W", "avg_W_err",
    "avg_eps", "avg_eps_err", "avg_theta", "avg_theta_err"
]
# Read the existing CSV file
rows = []
try:
    with open(avg_kinematics_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = list(csv_reader)
except FileNotFoundError:
    rows = []

# Prepare new rows for this Physics_Setting
new_setting_rows = []
for i, (tmin, tmax) in enumerate(zip(t_min_values, t_max_values), start=1):
    tbin_key = f"t({tmin:.2f}-{tmax:.2f})"
    if tbin_key in avg_kin_tbin:
        avg_kin = avg_kin_tbin[tbin_key]
        new_setting_rows.append({
            "Physics_Setting": PHY_SETTING,
            "total_tbins": total_tbins,
            "tbin_number": i,
            "t_min": f"{tmin:.4f}",
            "t_max": f"{tmax:.4f}",
            "t_central": f"{(tmin + tmax) / 2:.4f}",
            "total_phibins": total_phibins,  # This should match the header
            "phi_central": f"{(0 + 360) / 2:.4f}",
            "avg_Q2": f"{avg_kin['avg_Q2']:.6f}",
            "avg_Q2_err": f"{avg_kin['err_Q2']:.6f}",
            "avg_W": f"{avg_kin['avg_W']:.6f}",
            "avg_W_err": f"{avg_kin['err_W']:.6f}",
            "avg_eps": f"{avg_kin['avg_eps']:.6f}",
            "avg_eps_err": f"{avg_kin['err_eps']:.6f}",
            "avg_theta": f"{avg_kin['avg_theta']:.6f}",
            "avg_theta_err": f"{avg_kin['err_theta']:.6f}"
        })

# Check if Physics_Setting already exists
setting_exists = any(
    row["Physics_Setting"] == PHY_SETTING and
    float(row["t_min"]) == float(tmin) and
    float(row["t_max"]) == float(tmax)
    for row in rows
)
# Replace or append
if setting_exists:
    rows = [row for row in rows if row["Physics_Setting"] != PHY_SETTING]
    rows.extend(new_setting_rows)
else:
    rows.extend(new_setting_rows)

# Write everything back to CSV
with open(avg_kinematics_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=avg_kin_header)
    writer.writeheader()
    writer.writerows(rows)

print(f"Average Q² and W values updated for '{PHY_SETTING}' in: {avg_kinematics_path}")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SIMC Yield Calculations
N_simc_raw_MMpi = {}

# Loop over all t-bins and phi-bins
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:
    # Create a unique key for the current t-bin and phi-bin
    bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"

    # Initialize the error array
    dN_simc_raw = array.array('d', [0.0])

    # Get the dummysub histogram for the current bin
    hist = hist_set["simc_raw"]["MMpi"]  # Replace "MMpi" with the variable of interest

    # Calculate the integral and error
    N_simc_raw = hist.IntegralAndError(1, nbins, dN_simc_raw, "")

    # Store the results in the dictionary
    N_simc_raw_MMpi[bin_key] = {
        "N_simc_raw": N_simc_raw,
        "error_raw": dN_simc_raw[0]
    }

# Dictionary to store integrals and errors for each t-bin and phi-bin
dN_simc_MMpi = {}

# Loop over all t-bins and phi-bins
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:
    # Create a unique key for the current t-bin and phi-bin
    bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"

    # Initialize the error array
    dN_simc = array.array('d', [0.0])

    # Get the dummysub histogram for the current bin
    hist = hist_set["norm_simc"]["MMpi"]  # Replace "MMpi" with the variable of interest

    # Calculate the integral and error
    N_simc = hist.IntegralAndError(1, nbins, dN_simc, "")

    # Store the results in the dictionary
    dN_simc_MMpi[bin_key] = {
        "integral": N_simc,
        "error": dN_simc[0]
    }

    # Print the result for debugging
#    print(f"SIMC Yield: {bin_key} = Yield: {N_simc:.8f}, Error: {dN_simc[0]:.8f}")
#print("-"*40)

# Dictionary to store integrals and errors for each t-bin and phi-bin
dN_error = {}

# Loop over all t-bins and phi-bins
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:
    # Create a unique key for the current t-bin and phi-bin
    bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"

    # Initialize error arrays
    dN_simc_error = array.array('d', [0.0])

    # Get the randsub_data and randsub_dummy histograms for the current bin
    simc_hist = hist_set["simc"]["MMpi"]  # Replace "MMpi" with the variable of interest

    # Calculate the integrals and errors
    N_simc = simc_hist.IntegralAndError(1, nbins, dN_simc_error, "")

    # Normalize the integrals
    N_simc_norm = N_simc * normfac_simc

    # Calculate the normalized errors
    if N_simc != 0:
        dN_simc_norm = N_simc_norm * (dN_simc_error[0]/N_simc)
    else:
        dN_simc_norm = N_simc_norm * (0)

    # Ensure dN_data_norm and dN_dummy_norm are not zero
#    if dN_simc_norm == 0:
#        dN_simc_norm = 1000.0

    # Store the results in the dictionary
    dN_error[bin_key] = {
        "N_simc": N_simc,
        "N_simc_norm": N_simc_norm,
        "dN_simc_norm": dN_simc_norm
    }

    # Print the result for debugging
#    print(f"SIMC Yield: {bin_key} = Yield: {N_simc_norm:.8f}, Error: {dN_simc_norm:.8f}")

# Print the results from both loops
print("=" * 40)
print("SIMC Yield Results:")
for bin_key in dN_simc_MMpi:
    # Use the correct key names from dN_data_MMP
    N_simc_norm = dN_simc_MMpi[bin_key]["integral"]
    dN_simc_norm = dN_error[bin_key]["dN_simc_norm"]
    print(f"SIMC Yield: {bin_key} = Yield: {N_simc_norm:.8f}, Error: {dN_simc_norm:.8f}")
print("=" * 40)

# Define the output CSV file name
yields_pions_path =  "%s/LTSEP_ANALYSIS/LTSep_CSVs/simc_yields_csv/%s/%s_Physics_SIMC_Yield.csv" % (REPLAYPATH, physet_dir_name, setting_name)

# Updated header with tbin_number and phibin_number included
header = ["Physics_Setting", "tbin_number", "t_min", "t_max", "phibin_number", "phi_min", "phi_max", "Counts", "simc_yield", "simc_yield_error", "%error/yield"]

# Read the existing CSV file into memory
rows = []
try:
    with open(yields_pions_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = list(csv_reader)
except FileNotFoundError:
    rows = []

# Gather all unique t-bin and phi-bin edges
tbin_edges = sorted(list({(tmin, tmax) for tmin, tmax, _, _, _ in data_histograms_by_tphi_cut}))

# Update or add rows
updated_rows = []
for tbin_number, (tmin, tmax) in enumerate(tbin_edges, start=1):
    # Filter phi bins for this t range, and sort
    phibins = sorted([(phimin, phimax) for tmin2, tmax2, phimin, phimax, _ in data_histograms_by_tphi_cut if tmin2 == tmin and tmax2 == tmax])
    for phibin_number, (phimin, phimax) in enumerate(phibins, start=1):
        bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"
        N_simc_raw = N_simc_raw_MMpi[bin_key]["N_simc_raw"]
        N_simc_norm = dN_simc_MMpi[bin_key]["integral"]
        dN_simc_norm = dN_error[bin_key]["dN_simc_norm"]

        # See if the row already exists and update it, else append new
        row_found = False
        for row in rows:
            if (
                row["Physics_Setting"] == PHY_SETTING
                and float(row["t_min"]) == tmin
                and float(row["t_max"]) == tmax
                and float(row["phi_min"]) == phimin
                and float(row["phi_max"]) == phimax
            ):
                row.update({
                    "Physics_Setting": PHY_SETTING,
                    "tbin_number": tbin_number,
                    "t_min": tmin,
                    "t_max": tmax,
                    "phibin_number": phibin_number,
                    "phi_min": phimin,
                    "phi_max": phimax,
                    "Counts": N_simc_raw,
                    "simc_yield": N_simc_norm,
                    "simc_yield_error": dN_simc_norm,
                    "%error/yield": (dN_simc_norm/N_simc_norm)*100 if N_simc_norm != 0 else dN_simc_norm
                })
                row_found = True
                break

        if not row_found:
            updated_rows.append({
                "Physics_Setting": PHY_SETTING,
                "tbin_number": tbin_number,
                "t_min": tmin,
                "t_max": tmax,
                "phibin_number": phibin_number,
                "phi_min": phimin,
                "phi_max": phimax,
                "Counts": N_simc_raw,
                "simc_yield": N_simc_norm,
                "simc_yield_error": dN_simc_norm,
                "%error/yield": (dN_simc_norm/N_simc_norm)*100 if N_simc_norm != 0 else dN_simc_norm
            })

# Combine updated rows with existing rows
rows.extend(updated_rows)

# Write everything back to the CSV
with open(yields_pions_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(rows)

print(f"simc yield results written to {yields_pions_path}")

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/OUTPUT/Analysis/PionLT/simc_yield/%s/%s_%s_Yield_SIMC.root" % (VOLATILEPATH, physet_dir_name, PHY_SETTING, MaxEvent) , "RECREATE")

# Create top-level directories for Q2 and W
q2_dir = outHistFile.mkdir("Cut_Q2_Pion_Events_Norm_SIMC")
w_dir = outHistFile.mkdir("Cut_W_Pion_Events_Norm_SIMC")
epsilon_dir = outHistFile.mkdir("Cut_Epsilon_Pion_Events_Norm_SIMC")
theta_dir = outHistFile.mkdir("Cut_Theta_Pion_Events_Norm_SIMC")

# Write histograms to t-bin-specific directories
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:

    # --- Write MMpi (binned in t and phi) ---
    t_dirname_MM = f"Cut_MM_Pion_Events_Norm_t{tmin:.2f}-t{tmax:.2f}_SIMC".replace(".", "p")
    tdir_MM = outHistFile.GetDirectory(t_dirname_MM)
    if not tdir_MM:
        tdir_MM = outHistFile.mkdir(t_dirname_MM)
    tdir_MM.cd()
    var_MM = "MMpi"
    for cat in ["norm_simc"]:
        hist = hist_set[cat].get(var_MM)
        if hist:
            hist.SetTitle(f"{var_MM}_t{tmin:.2f}-t{tmax:.2f}_phi{phimin}-phi{phimax}_{cat}".replace(".", "p"))
            hist.SetName(f"{var_MM}_t{tmin:.2f}-t{tmax:.2f}_phi{phimin}-phi{phimax}_{cat}".replace(".", "p"))
            hist.GetYaxis().SetTitle("Counts")
            hist.Write()

    # --- Write Q2 and W once per t-bin (no phi binning) ---
    for var, target_dir in zip(["Q2", "W", "epsilon", "theta"], [q2_dir, w_dir, epsilon_dir, theta_dir]):
        for cat in ["norm_simc"]:
            hist = hist_set[cat].get(var)
            if hist:
                target_dir.cd()
                hist.SetTitle(f"{var}_t{tmin:.2f}-t{tmax:.2f}_{cat}".replace(".", "p"))
                hist.SetName(f"{var}_t{tmin:.2f}-t{tmax:.2f}_{cat}".replace(".", "p"))
                hist.GetYaxis().SetTitle("Counts")
                hist.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_SIMC.Close()
outHistFile.Close()

print ("Processing Complete")