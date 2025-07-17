#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-07-15 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
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
from functools import reduce
import math as ma
import uncertainties as u

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - CSV file suffix \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
ITERATION = sys.argv[2]
SIMC_Suffix = sys.argv[3]

# Ensure ITERATION is greater than 01
if int(ITERATION) <= 1:
    print("!!!!! ERROR !!!!!\n ITERATION must be greater than 01\n!!!!! ERROR !!!!!")
    sys.exit(1)

# Creating argument for previous iteration
#ITERATION_PREV = f"{int(ITERATION)-1:02d}"
ITERATION_PREV = f"{ITERATION}"

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
XSECT_PARAMPATH = "%s/scripts/ltsep_analysis/src/fit_params" % (UTILPATH)
XSECT_OUTPATH = "%s/scripts/ltsep_analysis/src/output" % (UTILPATH)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/OUTPUT/Analysis/SIMC/" % (USER)

##################################################################################################################################################

# Read stuff from the main event tree
rootFile_SIMC = "%s/%s_iter%s/%s.root" % (SIMCPATH, PHY_SETTING, ITERATION_PREV, SIMC_Suffix)
infile_SIMC = ROOT.TFile.Open(rootFile_SIMC, "READ")
TBRANCH_SIMC = infile_SIMC.Get("h10")

# Loop over events and recalculate weights using model equations
total_events = TBRANCH_SIMC.GetEntries()

# Dynamically grab all branches from input tree
branch_arrays = {}
branch_types = {}
branches = TBRANCH_SIMC.GetListOfBranches()
for branch in branches:
    bname = branch.GetName()
    # Rename branches for grabbing
    if bname == 'sigcm':
        grab_bname = f'sigcm_prev_raw'
    elif bname == 'Weight':
        grab_bname = f'Weight_prev_raw'
    else:
        grab_bname = bname
    branch_arrays[grab_bname] = array.array('f', [0])
    branch_types[grab_bname] = 'F'
    TBRANCH_SIMC.SetBranchAddress(bname, branch_arrays[grab_bname])

# Create output file for new weights
outFile_SIMC = ROOT.TFile.Open("%s/%s.root" % (SIMCPATH, SIMC_Suffix.replace(".root","")), "RECREATE")
new_TBRANCH_SIMC = ROOT.TTree("h10", "Iteration_%s" % (ITERATION))
# Create all output branches same as input file
output_arrays = {}

# Create all output branches same as input file
for bname in branch_arrays:
    output_arrays[bname] = array.array('d', [0])
    new_TBRANCH_SIMC.Branch(bname, output_arrays[bname], f"{bname}/D")

# Add four more branches: Weight and sigcm
output_arrays[f'sigcm_prev_iter{ITERATION_PREV}'] = array.array('d', [0])
output_arrays[f'Weight_prev_iter{ITERATION_PREV}'] = array.array('d', [0])
output_arrays['Weight'] = array.array('d', [0])
output_arrays['sigcm'] = array.array('d', [0])
new_TBRANCH_SIMC.Branch(f'sigcm_prev_iter{ITERATION_PREV}', output_arrays[f'sigcm_prev_iter{ITERATION_PREV}'], f"sigcm_prev_iter{ITERATION_PREV}/D")
new_TBRANCH_SIMC.Branch(f'Weight_prev_iter{ITERATION_PREV}', output_arrays[f'Weight_prev_iter{ITERATION_PREV}'], f"Weight_prev_iter{ITERATION_PREV}/D")
new_TBRANCH_SIMC.Branch('Weight', output_arrays['Weight'], 'Weight/D')
new_TBRANCH_SIMC.Branch('sigcm', output_arrays['sigcm'], 'sigcm/D')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Fixed variables 
mp = 0.93827231  # Proton mass in GeV
m_pi = 0.13956995  # Charged pion mass in GeV
# Extract Q2 value from PHY_SETTING argument
Q2par = 0
Q2str = PHY_SETTING.split('_')[0]
if Q2str.startswith('Q'):
    Q2par = Q2str[1:].replace('p', '') 

if ITERATION_PREV >= "01":
    # Read fit parameters from file and assign to variables (first column only)
    param_file = "%s/iter%s/par.pl_%s" % (XSECT_PARAMPATH, ITERATION_PREV, Q2par)
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
    print(f"p1={p1}, p2={p2}, p3={p3}, p4={p4}, p5={p5}, p6={p6}, p7={p7}, p8={p8}, p9={p9}, p10={p10}, p11={p11}, p12={p12}, p13={p13}, p14={p14}, p15={p15}, p16={p16}")
    print ("Fit parameters loaded from previous iteration file for cross-section calculation:")

    # Store per-event sigcm_prev_iter values in a list for later use
    sigcm_prev_iter_list = []
    weight_prev_iter_list = []
    # Preparing variables for previous iteration cross-section calculation
    for i in range(total_events):
        TBRANCH_SIMC.GetEntry(i)
        Q2_evt = branch_arrays['Q2i'][0]
        W_evt = branch_arrays['Wi'][0]
        thetacm_rad_evt = branch_arrays['thetapq'][0] * (np.pi / 180.0)  # Convert to radians
        eps_evt = branch_arrays['epsilon'][0]
        phi_rad_evt = branch_arrays['phipqi'][0] * (np.pi / 180.0)  # Convert to radians
        abs_t_evt = branch_arrays['ti'][0]
        wt_prev_raw = branch_arrays['Weight_prev_raw'][0]
        sigcm_prev_evt = branch_arrays['sigcm_prev_raw'][0]

        # Funtional form for separated cross-section calculation
        sigma_T = np.float128((p1 / Q2_evt) + (p2 / (Q2_evt**2)))
        sigma_L = np.float128((p5 * Q2_evt * np.exp((p6 - p7 * np.log(Q2_evt)) * np.abs(abs_t_evt))) / (1 + p8 * Q2_evt + p9 * (Q2_evt**2))**2)
        sigma_LT = np.float128((np.exp(p10 + (p11 * np.abs(abs_t_evt) / np.sqrt(Q2_evt))) + p12 + (p13 / (Q2_evt**2))) * np.sin(thetacm_rad_evt))
        sigma_TT = np.float128((((p14 / (Q2_evt**2)) * (np.abs(abs_t_evt) / (np.abs(abs_t_evt) + m_pi**2)**2) * np.sin(thetacm_rad_evt)**2)))

        wfactor = np.float128(1.0/(W_evt**2 - mp**2)**2)

        diff_xsect_prev_iter = np.float128((1 / (2 * np.pi)) * (eps_evt * sigma_L + sigma_T + np.sqrt(2 * eps_evt * (eps_evt + 1)) * sigma_LT * np.cos(phi_rad_evt) + eps_evt * sigma_TT * np.cos(2 * phi_rad_evt)))
        diff_xsect_prev_iter *= wfactor

        mod_diff_xsect_prev_iter = np.float128(diff_xsect_prev_iter / 1e6)  # Convert to microbarns/MeV**2/rad

        # Print all input variables and cross section for first three events
        if i < 3:
            print(f"Event {i}: Q2={Q2_evt}, W={W_evt}, theta_cm(rad)={thetacm_rad_evt}, epsilon={eps_evt}, phi_cm={phi_rad_evt}, t={abs_t_evt}, cross_section={mod_diff_xsect_prev_iter}")

        # Calculate previous weight
        wt_prev_iter = np.float128(wt_prev_raw * (mod_diff_xsect_prev_iter / sigcm_prev_evt))

        # Save recalculated cross-section to the new output branch for each event
        output_arrays[f'sigcm_prev_iter{ITERATION_PREV}'][0] = np.float128(mod_diff_xsect_prev_iter)
        output_arrays[f'Weight_prev_iter{ITERATION_PREV}'][0] = np.float128(wt_prev_iter)
        sigcm_prev_iter_list.append(np.float128(mod_diff_xsect_prev_iter))
        weight_prev_iter_list.append(np.float128(wt_prev_iter))

        if i < 3:
            print(f"sigcm_prev_iter{ITERATION_PREV} for event {i}: {output_arrays[f'sigcm_prev_iter{ITERATION_PREV}'][0]}, raw: {branch_arrays['sigcm_prev_raw'][0]}")

else:
    print(f"cross-section calculation for previous iteration is not required")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if ITERATION >= "01":
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
    print(f"p1={p1}, p2={p2}, p3={p3}, p4={p4}, p5={p5}, p6={p6}, p7={p7}, p8={p8}, p9={p9}, p10={p10}, p11={p11}, p12={p12}, p13={p13}, p14={p14}, p15={p15}, p16={p16}")
    print ("Fit parameters loaded from new iteration file:")

    # Store per-event sigcm and Weight values in lists for later use
    sigcm_list = []
    weight_list = []
    # Preparing variables for previous iteration cross-section calculation
    for i in range(total_events):
        TBRANCH_SIMC.GetEntry(i)
        Q2_evt = branch_arrays['Q2i'][0]
        W_evt = branch_arrays['Wi'][0]
        thetacm_rad_evt = branch_arrays['thetapq'][0] * (np.pi / 180.0)  # Convert to radians
        eps_evt = branch_arrays['epsilon'][0]
        phi_rad_evt = branch_arrays['phipqi'][0] * (np.pi / 180.0)  # Convert to radians
        abs_t_evt = branch_arrays['ti'][0]
        # Use per-event value from previous iteration lists
        wt_prev_evt = weight_prev_iter_list[i]
        sigcm_prev_evt = sigcm_prev_iter_list[i]

        # Funtional form for separated cross-section calculation
        sigma_T = np.float128((p1 / Q2_evt) + (p2 / (Q2_evt**2)))
        sigma_L = np.float128((p5 * Q2_evt * np.exp((p6 - p7 * np.log(Q2_evt)) * np.abs(abs_t_evt))) / (1 + p8 * Q2_evt + p9 * (Q2_evt**2))**2)
        sigma_LT = np.float128((np.exp(p10 + (p11 * np.abs(abs_t_evt) / np.sqrt(Q2_evt))) + p12 + (p13 / (Q2_evt**2))) * np.sin(thetacm_rad_evt))
        sigma_TT = np.float128((((p14 / (Q2_evt**2)) * (np.abs(abs_t_evt) / (np.abs(abs_t_evt) + m_pi**2)**2) * np.sin(thetacm_rad_evt)**2)))

        wfactor = np.float128(1.0/(W_evt**2 - mp**2)**2)

        diff_xsect_new_iter = np.float128((1 / (2 * np.pi)) * (eps_evt * sigma_L + sigma_T + np.sqrt(2 * eps_evt * (eps_evt + 1)) * sigma_LT * np.cos(phi_rad_evt) + eps_evt * sigma_TT * np.cos(2 * phi_rad_evt)))

        diff_xsect_new_iter *= wfactor

        mod_diff_xsect_new_iter = np.float128(diff_xsect_new_iter / 1e6)  # Convert to microbarns/MeV**2/rad

        if i < 3:
            print(f"Event {i}: Q2={Q2_evt}, W={W_evt}, theta_cm(rad)={thetacm_rad_evt}, epsilon={eps_evt}, phi_cm={phi_rad_evt}, t={abs_t_evt}, cross_section_new={mod_diff_xsect_new_iter}")

        # Calculate new weight
        wt_new_iter = np.float128(wt_prev_evt * (mod_diff_xsect_new_iter / sigcm_prev_evt))

        # Store calculated cross-section and weight in output arrays for writing to output tree and append to lists
        output_arrays['sigcm'][0] = np.float128(mod_diff_xsect_new_iter)
        output_arrays['Weight'][0] = np.float128(wt_new_iter)
        sigcm_list.append(np.float128(mod_diff_xsect_new_iter))
        weight_list.append(np.float128(wt_new_iter))

        # Print comparison of previous and new sigcm for each event
        if i < 3:
            print(f"Event {i}: sigcm_prev_iter = {sigcm_prev_evt}, sigcm (new) = {output_arrays['sigcm'][0]}")
            print(f"Event {i}: Weight_prev_iter = {wt_prev_evt}, Weight (new) = {wt_new_iter}")

else:
    print(f"Iteration number is not correct")

################################################################################################################################################

# Track bad events and progress bar
bad_events = []
def progressBar(current, total, bar_length=25):
    percent = float(current) / total
    arrow = '-' * int(round(percent * bar_length)-1) + '>' if percent > 0 else ''
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f"\rProgress: [{arrow}{spaces}] {int(percent*100)}% ({current}/{total})")
    sys.stdout.flush()
    if current == total-1:
        print()

# Loop over all events and fill the new tree
for i in range(total_events):
    progressBar(i, total_events, bar_length=25)
    TBRANCH_SIMC.GetEntry(i)
    first_branch = list(branch_arrays.keys())[0]
    if branch_arrays[first_branch][0] == 0.0:
        bad_events.append(i)
        continue
    # Copy all branch values to output_arrays for this event
    for bname in branch_arrays:
        output_arrays[bname][0] = branch_arrays[bname][0]
    # Ensure sigcm_prev_iter{ITERATION}, sigcm, and Weight are set from the per-event lists before filling the tree
    if 'sigcm_prev_iter_list' in globals() or 'sigcm_prev_iter_list' in locals():
        output_arrays[f'sigcm_prev_iter{ITERATION_PREV}'][0] = sigcm_prev_iter_list[i]
    if 'weight_prev_iter_list' in globals() or 'weight_prev_iter_list' in locals():
        output_arrays[f'Weight_prev_iter{ITERATION_PREV}'][0] = weight_prev_iter_list[i]
    if 'sigcm_list' in globals() or 'sigcm_list' in locals():
        output_arrays['sigcm'][0] = sigcm_list[i]
    if 'weight_list' in globals() or 'weight_list' in locals():
        output_arrays['Weight'][0] = weight_list[i]
    new_TBRANCH_SIMC.Fill()
print(f"Total {total_events} events in the ROOT file.")
if bad_events:
    print(f": Number of bad events: {len(bad_events)}/{total_events}")

new_TBRANCH_SIMC.Write("h10",ROOT.TObject.kOverwrite)
outFile_SIMC.Close()
infile_SIMC.Close()
print(f"Reweighted SIMC file written to {outFile_SIMC.GetName()}")

print(f"\nThere were {len(bad_events)}/{total_events} bad events skipped...")

print("Processing Complete")