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
if len(sys.argv)-1!=9:
    print("!!!!! ERROR !!!!!\n Expected 9 arguments\n Usage is with - ROOTfileSuffixs Beam Energy MaxEvents RunList CVSFile\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
DATA_Suffix = sys.argv[3]
DUMMY_Suffix = sys.argv[4]
SIMC_Suffix = sys.argv[5]
DATA_RUN_LIST = sys.argv[6]
DUMMY_RUN_LIST = sys.argv[7]
CSV_FILE = sys.argv[8]
Iteration = sys.argv[9]

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
RUNLISTPATH = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_BATCH/InputRunLists/PionLT_2021_2022" % (USER)
EFF_CSV     = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/efficiencies" % (USER)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_iter%s" % (setting_name, Iteration)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/OUTPUT/Analysis/SIMC/%s/" % (USER, physet_dir_name)
MMCUT_CSV   = "%s/LTSEP_ANALYSIS/LTSep_CSVs/mm_offset_cut_csv/%s_std" % (REPLAYPATH, setting_name)
DCUT_CSV    = "%s/LTSEP_ANALYSIS/LTSep_CSVs/diamond_cut_csv/%s_std" % (REPLAYPATH, setting_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/LTSEP_ANALYSIS/src/plots/%s_%s_ProdCoin_Pion_Analysis_Ratio_Comparison_Distributions_iter%s.pdf" % (REPLAYPATH, PHY_SETTING, MaxEvent, Iteration)

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DATA_Suffix)
rootFile_DUMMY = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DUMMY_Suffix)
rootFile_SIMC = "%s/%s.root" % (SIMCPATH, SIMC_Suffix)
data_run_list = "%s/%s" % (RUNLISTPATH, DATA_RUN_LIST)
dummy_run_list = "%s/%s" % (RUNLISTPATH, DUMMY_RUN_LIST)
eff_csv_file = "%s/%s.csv" % (EFF_CSV, CSV_FILE)
mmcut_csv_file = "%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, setting_name)
dcut_csv_file = "%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, setting_name)

###############################################################################################################################################

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

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
    with open(mmcut_csv_file, mode='r', newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row[csv_reader.fieldnames[0]].strip() == PHY_SETTING:  # Match first column (Physics_Setting)
                MM_Offset = float(row["MM_Offset"].strip())
                MM_Cut_lowvalue = float(row["MM_Cut_low"].strip())
                MM_Cut_highvalue = float(row["MM_Cut_high"].strip())
                break
        else:
            raise ValueError(f"No matching Physics_Setting '{PHY_SETTING}' found in {mmcut_csv_file}.")

except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
print(f"MMpi_Offset = {MM_Offset:.6f}")
#print(f"MMpi_Cut_lowvalue = {MMpi_Cut_lowvalue}")
#print(f"MMpi_Cut_highvalue = {MMpi_Cut_highvalue}")
DATA_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset) <= MM_Cut_highvalue)
SIMC_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= event.missmass <= MM_Cut_highvalue)

# SIMC Cuts for Pions Selection
HMS_Acceptance = lambda event: (event.hsdelta >= -8.0) & (event.hsdelta <= 8.0) & (event.hsxpfp >= -0.08) & (event.hsxpfp <= 0.08) & (event.hsypfp >= -0.045) & (event.hsypfp <= 0.045)
SHMS_Acceptance = lambda event: (event.ssdelta >= -10.0) & (event.ssdelta <= 20.0) & (event.ssxpfp >= -0.06) & (event.ssxpfp <= 0.06) & (event.ssypfp >= -0.04) & (event.ssypfp <= 0.04)
#SHMS_Aero_Cut = lambda event: (event.paero_x_det > -55.0) & (event.paero_x_det < 55.0) & (event.paero_y_det > -50) & (event.paero_y_det < 50) # Aerogel tray n = 1.030
SHMS_Aero_Cut = lambda event: (event.paero_x_det > -45.0) & (event.paero_x_det < 45.0) & (event.paero_y_det > -30) & (event.paero_y_det < 30) # Aerogel tray n = 1.011

#########################################################################################################################################

# Section for grabing Prompt/Random selection parameters from PARAM filePARAMPATH = UTILPATH+"/DB/PARAM"
# Function to read run numbers from a given run list file
def read_run_numbers(run_list_file):
    try:
        with open(run_list_file, 'r') as f:
            return [int(line.strip()) for line in f if line.strip().isdigit()]
    except FileNotFoundError:
        print(f"!!!!! ERROR !!!!!\nRun list file not found: {run_list_file}\n!!!!! ERROR !!!!!")
        sys.exit(2)

# Read run numbers from the specified files
data_run_numbers = read_run_numbers(data_run_list)
dummy_run_numbers = read_run_numbers(dummy_run_list)

# Section for grabbing Prompt/Random selection parameters from PARAM file
PARAMPATH = UTILPATH+"/DB/PARAM/PionLT"
TimingCutFile = "%s/PionLT_Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!

TimingCutf = open(TimingCutFile)
try:
    TimingCutFile
except NameError:
    print("!!!!! ERRROR !!!!!\n One (or more) of the cut files not found!\n!!!!! ERRORR !!!!!")
    sys.exit(2)
print("\n Reading timing cuts from %s to select random windows \n" % TimingCutFile)

PromptWindow = [0, 0]
RandomWindows = [0, 0, 0, 0]
linenum = 0 # Count line number we're on
TempPar = -1 # To check later
for line in TimingCutf: # Read all lines in the cut file
    linenum += 1 # Add one to line number at start of loop
    if(linenum > 1): # Skip first line
        line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
        line = line.rstrip()
        array_list = line.split(",") # Convert line into an array_list, anything after a comma is a new entry

        # Check if any of the run numbers are within the specified range
        if any(run_num in range(int(array_list[0]), int(array_list[1]) + 1) for run_num in data_run_numbers + dummy_run_numbers):
            TempPar += 2  # If run number is in range, set to non -1 value
            BunchSpacing = float(array_list[2])
            CoinOffset = float(array_list[3])  # Coin offset value
            nSkip = float(array_list[4])  # Number of random windows skipped
            nWindows = float(array_list[5])  # Total number of random windows
            PromptPeak = float(array_list[6])  # Pion CT prompt peak position

if(TempPar == -1): # If value is still -1, run number provided din't match any ranges specified so exit
    print("!!!!! ERROR !!!!!\n Run number specified does not fall within a set of runs for which cuts are defined in %s\n!!!!! ERROR !!!!!" % TimingCutFile)
    sys.exit(3)
elif(TempPar > 1):
    print("!!! WARNING!!! Run number was found within the range of two (or more) line entries of %s !!! WARNING !!!" % TimingCutFile)
    print("The last matching entry will be treated as the input, you should ensure this is what you want")

# From our values from the file, reconstruct our windows
PromptWindow[0] = PromptPeak - (BunchSpacing/2) - CoinOffset
PromptWindow[1] = PromptPeak + (BunchSpacing/2) + CoinOffset
RandomWindows[0] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing) - ((nWindows/2)*BunchSpacing)
RandomWindows[1] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing)
RandomWindows[2] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) 
RandomWindows[3] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) + ((nWindows/2)*BunchSpacing)

###############################################################################################################################################

print("-"*40)
print ('Calculating effective charge for the data and dummy run list')
print("-"*40)

# Read CSV File and calculate total charge in mC for normalization
# Input runlists and csv files
data_run_list_file = (data_run_list)
dummy_run_list_file = (dummy_run_list)
eff_csv_file_name = (eff_csv_file)

# Read run numbers from the run list file
with open(data_run_list_file, 'r') as run_list_file_1:
    data_runs = [line.strip() for line in run_list_file_1 if line.strip()]
with open(dummy_run_list_file, 'r') as run_list_file_2:
    dummy_runs = [line.strip() for line in run_list_file_2 if line.strip()]

# Read CSV file using pandas
df = pd.read_csv(eff_csv_file_name)

# Filter DataFrame to include only rows with run numbers in the run list
filtered_data_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(data_runs)]
filtered_dummy_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(dummy_runs)]


# Intializing Variables to calculate the product of efficiency and charge
total_data_effective_charge = 0.0
total_dummy_effective_charge = 0.0
total_data_effective_charge_sum = 0.0
total_dummy_effective_charge_sum = 0.0
total_data_effective_charge_error_sum = 0.0
total_dummy_effective_charge_error_sum = 0.0
total_dummy_effective_charge_cal = 0.0
total_dummy_effective_charge_cal_error = 0.0

# Detector (HMS Cer + HMS Cal) Efficiencies
hms_Cer_detector_efficiency = 0.9981
hms_Cal_detector_efficiency = 0.9981
hms_Cer_detector_efficiency_error = 0.0001
hms_Cal_detector_efficiency_error = 0.0001

if setting_name == "Q3p85_W2p62_t0p21":
    #RF cut Efficiency
    RF_efficiency = 0.998
    RF_efficiency_error = 0.00007
    # Pion Absorption Correction
    pion_absorption_correction = 0.9654
    pion_absorption_correction_error = 0.00005
else:
    print("!!!!!\n Need to provide RF Eff and Absortion Corr for %s\n!!!!!" % (setting_name))
    RF_efficiency = 1.0
    RF_efficiency_error = 0.0
    pion_absorption_correction = 1.0
    pion_absorption_correction_error = 0.0
    print("Assuming RF Eff = 1.0 +/- 0.0 and Pion Absorption Corr = 1.0 +/- 0.0 \n")
    
row_data = []

# Print charge values for each run
for index, row in filtered_data_df.iterrows():
    data_charge = row['Corrected_Charge'] 
    data_charge_error = row['Corrected_Charge_ERROR']
    data_Boiling_factor = row['Target_BoilingCorr']
    data_Boiling_factor_error = row['Target_BoilingCorr_ERROR']
    data_edtm_livetime_Corr = row['Non_Scaler_EDTM_Live_Time_Corr']
    data_edtm_livetime_Corr_error = row['Non_Scaler_EDTM_Live_Time_Corr_ERROR']
    data_coinblocking_factor = row['CoinBlocking']
    data_coinblocking_factor_error = row['CoinBlocking_ERROR']

    data_hms_tracking_efficiency = row['HMS_Elec_SING_TRACK_EFF'] 
    data_shms_tracking_efficiency = row['SHMS_Pion_SING_TRACK_EFF'] 
    data_hms_tracking_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    data_shms_tracking_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']

    data_hms_hodo_3_of_4_efficiency = row['HMS_Hodo_3_of_4_EFF']
    data_shms_hodo_3_of_4_efficiency = row['SHMS_Hodo_3_of_4_EFF']
    data_hms_hodo_3_of_4_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    data_shms_hodo_3_of_4_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']
#    data_hms_hodo_3_of_4_efficiency_error = row['HMS_Hodo_3_of_4_EFF_ERROR']
#    data_shms_hodo_3_of_4_efficiency_error = row['SHMS_Hodo_3_of_4_EFF_ERROR']

    data_shms_aero_detector_efficiency = row['SHMS_Aero_COIN_Pion_Eff']
    data_shms_aero_detector_efficiency_error = row['SHMS_Aero_COIN_Pion_Eff_ERROR']

    data_product = (data_charge * pion_absorption_correction * data_hms_tracking_efficiency * data_shms_tracking_efficiency * RF_efficiency * hms_Cer_detector_efficiency * hms_Cal_detector_efficiency * data_hms_hodo_3_of_4_efficiency * data_shms_hodo_3_of_4_efficiency * data_shms_aero_detector_efficiency * data_edtm_livetime_Corr * data_Boiling_factor * data_coinblocking_factor)
    data_product_error = data_product * (math.sqrt((data_charge_error/data_charge)** 2 +  (pion_absorption_correction_error/pion_absorption_correction)**2 + (data_hms_tracking_efficiency_error/data_hms_tracking_efficiency)** 2 + (data_shms_tracking_efficiency_error/data_shms_tracking_efficiency)** 2 + (RF_efficiency_error/RF_efficiency)**2 + (data_edtm_livetime_Corr_error/data_edtm_livetime_Corr)** 2 + (hms_Cer_detector_efficiency_error/hms_Cer_detector_efficiency)** 2 + (hms_Cal_detector_efficiency_error/hms_Cal_detector_efficiency)** 2 + (data_hms_hodo_3_of_4_efficiency_error/data_hms_hodo_3_of_4_efficiency)** 2 + (data_shms_hodo_3_of_4_efficiency_error/data_shms_hodo_3_of_4_efficiency)**2 + (data_shms_aero_detector_efficiency_error/data_shms_aero_detector_efficiency)**2 + (data_Boiling_factor_error/data_Boiling_factor)** 2 + (data_coinblocking_factor_error/data_coinblocking_factor)**2)) 

    total_data_effective_charge_sum += data_product
    total_data_effective_charge_error_sum += (data_product_error)** 2

    print('Charge for data run: {:<10} Data Charge: {:.3f} Pion_Absorpt_Corr: {:.3f} HMS Tracking Eff: {:.3f} SHMS Tracking Eff: {:.3f} RF Eff: {:.3f} HMS Cer Detector Eff: {:.3f} HMS Cal Detector_Eff: {:.3f} HMS Hodo 3/4 Eff: {:.3f} SHMS Hodo 3/4 Eff: {:.3f} SHMS Aero Detector Eff: {:.3f} EDTM Live Time: {:.3f} Boiling Correction: {:.3f} Coin Blocking: {:.3f} Product: {:.3f}'.format(row["Run_Number"], data_charge, pion_absorption_correction, data_hms_tracking_efficiency, data_shms_tracking_efficiency, RF_efficiency, hms_Cer_detector_efficiency, hms_Cal_detector_efficiency, data_hms_hodo_3_of_4_efficiency, data_shms_hodo_3_of_4_efficiency, data_shms_aero_detector_efficiency, data_edtm_livetime_Corr, data_Boiling_factor, data_coinblocking_factor, data_product))

    row_data.append({'Run_Number':row["Run_Number"], 'charge':data_charge, 'charge_error':data_charge_error, 'pion_absorption_correction':pion_absorption_correction, 'HMS_Tracking_Eff':data_hms_tracking_efficiency, 'HMS_Tracking_Eff_error':data_hms_tracking_efficiency_error, 'SHMS_Tracking_Eff':data_shms_tracking_efficiency, 'SHMS_Tracking_Eff_error':data_shms_tracking_efficiency_error, 'RF_Eff':RF_efficiency, 'RF_Eff_error':RF_efficiency_error, 'HMS_Cer_Detector_Eff':hms_Cer_detector_efficiency, 'HMS_Cer_Detector_Eff_error':hms_Cer_detector_efficiency_error, 'HMS_Cal_Detector_Eff':hms_Cal_detector_efficiency, 'HMS_Cal_Detector_Eff_error':hms_Cal_detector_efficiency_error, 'HMS_Hodo_3_4_Eff':data_hms_hodo_3_of_4_efficiency, 'HMS_Hodo_3_4_Eff_error':data_hms_hodo_3_of_4_efficiency_error, 'SHMS_Hodo_3_4_Eff':data_shms_hodo_3_of_4_efficiency, 'SHMS_Hodo_3_4_Eff_error':data_shms_hodo_3_of_4_efficiency_error, 'SHMS Aerogel Eff':data_shms_aero_detector_efficiency, 'SHMS Aerogel Eff_error':data_shms_aero_detector_efficiency_error, 'EDTM_Live_Time':data_edtm_livetime_Corr, 'EDTM_Live_Time_error':data_edtm_livetime_Corr_error, 'Boiling_factor':data_Boiling_factor, 'Boiling_factor_error':data_Boiling_factor_error, 'Coin_Blocking':data_coinblocking_factor, 'Coin_Blocking_error':data_coinblocking_factor_error, 'effective_charge':data_product, 'effective_charge_error':data_product_error})

print("-"*40)

total_data_effective_charge = total_data_effective_charge_sum
total_data_effective_charge_error = math.sqrt(total_data_effective_charge_error_sum)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

row_dummy = []

for index, row in filtered_dummy_df.iterrows():
    dummy_charge = row['Corrected_Charge'] 
    dummy_charge_error = row['Corrected_Charge_ERROR']
    dummy_edtm_livetime_Corr = row['Non_Scaler_EDTM_Live_Time_Corr']
    dummy_edtm_livetime_Corr_error = row['Non_Scaler_EDTM_Live_Time_Corr_ERROR']
    dummy_coinblocking_factor = row['CoinBlocking']
    dummy_coinblocking_factor_error = row['CoinBlocking_ERROR']

    dummy_hms_tracking_efficiency = row['HMS_Elec_SING_TRACK_EFF'] 
    dummy_shms_tracking_efficiency = row['SHMS_Pion_SING_TRACK_EFF'] 
    dummy_hms_tracking_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    dummy_shms_tracking_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']

    dummy_hms_hodo_3_of_4_efficiency = row['HMS_Hodo_3_of_4_EFF']
    dummy_shms_hodo_3_of_4_efficiency = row['SHMS_Hodo_3_of_4_EFF']
    dummy_hms_hodo_3_of_4_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    dummy_shms_hodo_3_of_4_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']
#    dummy_hms_hodo_3_of_4_efficiency_error = row['HMS_Hodo_3_of_4_EFF_ERROR']
#    dummy_shms_hodo_3_of_4_efficiency_error = row['SHMS_Hodo_3_of_4_EFF_ERROR']

    dummy_shms_aero_detector_efficiency = row['SHMS_Aero_COIN_Pion_Eff']
    dummy_shms_aero_detector_efficiency_error = row['SHMS_Aero_COIN_Pion_Eff_ERROR']

    dummy_product = (dummy_charge * dummy_hms_tracking_efficiency * dummy_shms_tracking_efficiency * RF_efficiency * hms_Cer_detector_efficiency * hms_Cal_detector_efficiency * dummy_hms_hodo_3_of_4_efficiency * dummy_shms_hodo_3_of_4_efficiency * dummy_edtm_livetime_Corr * dummy_shms_aero_detector_efficiency * dummy_coinblocking_factor)
    dummy_product_error = dummy_product * (math.sqrt(((dummy_charge_error/dummy_charge) ** 2 + dummy_hms_tracking_efficiency_error/dummy_hms_tracking_efficiency) ** 2 + (dummy_shms_tracking_efficiency_error/dummy_shms_tracking_efficiency) ** 2 + (RF_efficiency_error/RF_efficiency) **2 + (dummy_edtm_livetime_Corr_error/dummy_edtm_livetime_Corr)** 2 + (hms_Cer_detector_efficiency_error/hms_Cer_detector_efficiency) ** 2 + (hms_Cal_detector_efficiency_error/hms_Cal_detector_efficiency) ** 2 + (dummy_hms_hodo_3_of_4_efficiency_error/dummy_hms_hodo_3_of_4_efficiency) ** 2 + (dummy_shms_hodo_3_of_4_efficiency_error/dummy_shms_hodo_3_of_4_efficiency) ** 2 + (dummy_shms_aero_detector_efficiency_error/dummy_shms_aero_detector_efficiency) ** 2 + (dummy_coinblocking_factor_error/dummy_coinblocking_factor)**2))

    total_dummy_effective_charge_sum += dummy_product
    total_dummy_effective_charge_error_sum += (dummy_product_error)** 2

    print('Charge for dummy run: {:<10} Dummy Charge: {:.3f} HMS Tracking Eff: {:.3f} SHMS Tracking Eff: {:.3f} RF Eff: {:.3f} HMS Cer Detector Eff: {:.3f} HMS Cal Detector_Eff: {:.3f} HMS Hodo 3/4 Eff: {:.3f} SHMS Hodo 3/4 Eff: {:.3f} EDTM Live Time: {:.3f} SHMS Aerogel Eff: {:.3f} Coin Blocking: {:.3f} Product: {:.3f}'.format(row["Run_Number"], dummy_charge, dummy_hms_tracking_efficiency, dummy_shms_tracking_efficiency, RF_efficiency, hms_Cer_detector_efficiency, hms_Cal_detector_efficiency, dummy_hms_hodo_3_of_4_efficiency, dummy_shms_hodo_3_of_4_efficiency, dummy_edtm_livetime_Corr, dummy_shms_aero_detector_efficiency, dummy_coinblocking_factor, dummy_product))

    row_dummy.append({'Run_Number':row["Run_Number"], 'charge':dummy_charge, 'charge_error':dummy_charge_error, 'HMS_Tracking_Eff':dummy_hms_tracking_efficiency, 'HMS_Tracking_Eff_error':dummy_hms_tracking_efficiency_error, 'SHMS_Tracking_Eff':dummy_shms_tracking_efficiency, 'SHMS_Tracking_Eff_error':dummy_shms_tracking_efficiency_error, 'RF_Eff':RF_efficiency, 'RF_Eff_error':RF_efficiency_error, 'HMS_Cer_Detector_Eff':hms_Cer_detector_efficiency, 'HMS_Cer_Detector_Eff_error':hms_Cer_detector_efficiency_error, 'HMS_Cal_Detector_Eff':hms_Cal_detector_efficiency, 'HMS_Cal_Detector_Eff_error':hms_Cal_detector_efficiency_error, 'HMS_Hodo_3_4_Eff':dummy_hms_hodo_3_of_4_efficiency, 'HMS_Hodo_3_4_Eff_error':dummy_hms_hodo_3_of_4_efficiency_error, 'SHMS_Hodo_3_4_Eff':dummy_shms_hodo_3_of_4_efficiency, 'SHMS_Hodo_3_4_Eff_error':dummy_shms_hodo_3_of_4_efficiency_error, 'SHMS Aerogel Eff':dummy_shms_aero_detector_efficiency, 'SHMS Aerogel Eff_error':dummy_shms_aero_detector_efficiency_error, 'EDTM_Live_Time':dummy_edtm_livetime_Corr, 'EDTM_Live_Time_error':dummy_edtm_livetime_Corr_error, 'Boiling_factor':'NA', 'Boiling_factor_error':'NA', 'Coin_Blocking':dummy_coinblocking_factor, 'Coin_Blocking_error':dummy_coinblocking_factor_error, 'effective_charge':dummy_product, 'effective_charge_error':dummy_product_error})

print("-"*40)

# Dummy Target Thickness Correction
#dummy_target_corr = 4.8579 # KaonLT
dummy_target_corr = 3.527 # PionLT
dummy_target_corr_error = 0.227 # PionLT

total_dummy_effective_charge_cal = total_dummy_effective_charge_sum
total_dummy_effective_charge_cal_error = math.sqrt(total_dummy_effective_charge_error_sum)

total_dummy_effective_charge = (total_dummy_effective_charge_cal * dummy_target_corr)
total_dummy_effective_charge_error = total_dummy_effective_charge * math.sqrt((total_dummy_effective_charge_cal_error/total_dummy_effective_charge_cal)**2 + (dummy_target_corr_error/dummy_target_corr)** 2)

print("\nTotal effective charge for the data run list: {:.5f} ± {:.5f}".format(total_data_effective_charge, total_data_effective_charge_error))
print("\nTotal effective charge for the dummy run list: {:.5f} ± {:.5f}".format(total_dummy_effective_charge, total_dummy_effective_charge_error))

# Normalization factor Calculation for Data and Dummy
normfac_data = 1.0/(total_data_effective_charge)
normfac_dummy = 1.0/(total_dummy_effective_charge)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
print ("normfac_data :", normfac_data)
print ("normfac_dummy: ", normfac_dummy)
print ("normfac_simc: ", normfac_simc)
print("-"*40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define the output CSV file name
eff_cal_output_csv_path = "%s/LTSEP_ANALYSIS/LTSep_CSVs/datasimc_ratios_csv/%s/%s_Physics_DataSIMC_Norm_Eff_Charge_Calculation.csv" % (REPLAYPATH, physet_dir_name, PHY_SETTING)

# Prepare the header for the CSV file
header = [
    "Physics Setting", "Run Type", "Run Number", "Charge", "Charge Error",
    "HMS Tracking Eff", "HMS Tracking Eff Error", "SHMS Tracking Eff", "SHMS Tracking Eff Error",
    "HMS Cer Detector Eff", "HMS Cer Detector Eff Error", "HMS Cal Detector Eff", "HMS Cal Detector Eff Error",
    "HMS Hodo 3/4 Eff", "HMS Hodo 3/4 Eff Error", "SHMS Hodo 3/4 Eff", "SHMS Hodo 3/4 Eff Error",
    "EDTM Live Time", "EDTM Live Time Error", "Boiling Factor", "Boiling Factor Error", "Dummy Target Correction",
    "Dummy Target Correction Error", "Effective Charge", "Effective Charge Error"
]

# Combine data and dummy rows
combined_rows = []

# Add data rows
for row in row_data:
    combined_rows.append({
        "Physics Setting": PHY_SETTING,
        "Run Type": "data",
        "Run Number": row["Run_Number"],
        "Charge": row["charge"],
        "Charge Error": row["charge_error"],
        "HMS Tracking Eff": row["HMS_Tracking_Eff"],
        "HMS Tracking Eff Error": row["HMS_Tracking_Eff_error"],
        "SHMS Tracking Eff": row["SHMS_Tracking_Eff"],
        "SHMS Tracking Eff Error": row["SHMS_Tracking_Eff_error"],
        "HMS Cer Detector Eff": row["HMS_Cer_Detector_Eff"],
        "HMS Cer Detector Eff Error": row["HMS_Cer_Detector_Eff_error"],
        "HMS Cal Detector Eff": row["HMS_Cal_Detector_Eff"],
        "HMS Cal Detector Eff Error": row["HMS_Cal_Detector_Eff_error"],
        "HMS Hodo 3/4 Eff": row["HMS_Hodo_3_4_Eff"],
        "HMS Hodo 3/4 Eff Error": row["HMS_Hodo_3_4_Eff_error"],
        "SHMS Hodo 3/4 Eff": row["SHMS_Hodo_3_4_Eff"],
        "SHMS Hodo 3/4 Eff Error": row["SHMS_Hodo_3_4_Eff_error"],
        "EDTM Live Time": row["EDTM_Live_Time"],
        "EDTM Live Time Error": row["EDTM_Live_Time_error"],
        "Boiling Factor": row["Boiling_factor"],
        "Boiling Factor Error": row["Boiling_factor_error"],
        "Dummy Target Correction": "NA",
        "Dummy Target Correction Error": "NA",
        "Effective Charge": row["effective_charge"],
        "Effective Charge Error": row["effective_charge_error"],
    })

# Add dummy rows
for row in row_dummy:
    combined_rows.append({
        "Physics Setting": PHY_SETTING,
        "Run Type": "dummy",
        "Run Number": row["Run_Number"],
        "Charge": row["charge"],
        "Charge Error": row["charge_error"],
        "HMS Tracking Eff": row["HMS_Tracking_Eff"],
        "HMS Tracking Eff Error": row["HMS_Tracking_Eff_error"],
        "SHMS Tracking Eff": row["SHMS_Tracking_Eff"],
        "SHMS Tracking Eff Error": row["SHMS_Tracking_Eff_error"],
        "HMS Cer Detector Eff": row["HMS_Cer_Detector_Eff"],
        "HMS Cer Detector Eff Error": row["HMS_Cer_Detector_Eff_error"],
        "HMS Cal Detector Eff": row["HMS_Cal_Detector_Eff"],
        "HMS Cal Detector Eff Error": row["HMS_Cal_Detector_Eff_error"],
        "HMS Hodo 3/4 Eff": row["HMS_Hodo_3_4_Eff"],
        "HMS Hodo 3/4 Eff Error": row["HMS_Hodo_3_4_Eff_error"],
        "SHMS Hodo 3/4 Eff": row["SHMS_Hodo_3_4_Eff"],
        "SHMS Hodo 3/4 Eff Error": row["SHMS_Hodo_3_4_Eff_error"],
        "EDTM Live Time": row["EDTM_Live_Time"],
        "EDTM Live Time Error": row["EDTM_Live_Time_error"],
        "Boiling Factor": row["Boiling_factor"],
        "Boiling Factor Error": row["Boiling_factor_error"],
        "Dummy Target Correction": dummy_target_corr,
        "Dummy Target Correction Error": dummy_target_corr_error,
        "Effective Charge": row["effective_charge"],
        "Effective Charge Error": row["effective_charge_error"],
    })

# Write the combined rows directly to the CSV file
with open(eff_cal_output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(combined_rows)

    # Skip two lines
    csv_file.write('\n\n')

    # Write the summary in Variable, Value format
    summary_writer = csv.writer(csv_file)
    summary_writer.writerow(["Variable", "Value"])
    summary_writer.writerow(["Total Data Effective Charge", f"{total_data_effective_charge}"])
    summary_writer.writerow(["Total Data Effective Charge Error", f"{total_data_effective_charge_error}"])
    summary_writer.writerow(["Total Dummy Effective Charge", f"{total_dummy_effective_charge:.5f}"])
    summary_writer.writerow(["Total Dummy Effective Charge Error", f"{total_dummy_effective_charge_error}"])
    summary_writer.writerow(["Normalization Data Factor (Data)", f"{normfac_data}"])
    summary_writer.writerow(["Normalization Dummy Factor (Dummy)", f"{normfac_dummy}"])
    summary_writer.writerow(["Normalization SIMC Factor", f"{normfac_simc}"])

print(f"Effective charge calculation factors written to {eff_cal_output_csv_path}")

###################################################################################################################################################
nbins = 200
nbins_p = 1000


# Defining Histograms for Pions
# Histograms having Cuts (Acceptance)
H_gtr_beta_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_accpt_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_xp_pions_data_accpt_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_yp_pions_data_accpt_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_dp_pions_data_accpt_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_p_pions_data_accpt_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_data_accpt_cut_all = ROOT.TH1D("H_dc_x_fp_pions_data_accpt_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_data_accpt_cut_all = ROOT.TH1D("H_dc_y_fp_pions_data_accpt_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_data_accpt_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_data_accpt_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_data_accpt_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_data_accpt_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_data_accpt_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_data_accpt_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_data_accpt_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_accpt_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_accpt_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_accpt_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_accpt_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_accpt_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_xp_pions_data_accpt_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_yp_pions_data_accpt_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_dp_pions_data_accpt_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_p_pions_data_accpt_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_data_accpt_cut_all = ROOT.TH1D("P_dc_x_fp_pions_data_accpt_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_data_accpt_cut_all = ROOT.TH1D("P_dc_y_fp_pions_data_accpt_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_data_accpt_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_data_accpt_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_data_accpt_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_data_accpt_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_data_accpt_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_data_accpt_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_data_accpt_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_accpt_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_accpt_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_accpt_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_accpt_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_accpt_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_accpt_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_accpt_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_accpt_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_accpt_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_data_accpt_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_accpt_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_accpt_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_data_accpt_cut_all = ROOT.TH1D("MMpi_pions_data_accpt_cut_all", "MIssing Mass data (dummysub_accpt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_data_accpt_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_accpt_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_data_accpt_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_accpt_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_data_accpt_cut_all = ROOT.TH1D("pmiss_pions_data_accpt_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_data_accpt_cut_all = ROOT.TH1D("pmiss_x_pions_data_accpt_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_data_accpt_cut_all = ROOT.TH1D("pmiss_y_pions_data_accpt_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_data_accpt_cut_all = ROOT.TH1D("pmiss_z_pions_data_accpt_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_data_accpt_cut_all = ROOT.TH1D("emiss_pions_data_accpt_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_data_accpt_cut_all = ROOT.TH1D("W_pions_data_accpt_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_data_accpt_cut_all = ROOT.TH1D("Q2_pions_data_accpt_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_data_accpt_cut_all = ROOT.TH1D("Epsilon_pions_data_accpt_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_data_accpt_cut_all = ROOT.TH1D("t_pions_data_accpt_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_data_accpt_cut_all = ROOT.TH1D("ph_q_pions_data_accpt_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)
'''
# Histograms having Cuts (Acceptance + PID + RF)
H_gtr_beta_pions_data_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_data_cut_all = ROOT.TH1D("H_gtr_xp_pions_data_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_data_cut_all = ROOT.TH1D("H_gtr_yp_pions_data_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_data_cut_all = ROOT.TH1D("H_gtr_dp_pions_data_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_data_cut_all = ROOT.TH1D("H_gtr_p_pions_data_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_data_cut_all = ROOT.TH1D("H_dc_x_fp_pions_data_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_data_cut_all = ROOT.TH1D("H_dc_y_fp_pions_data_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_data_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_data_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_data_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_data_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_data_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_data_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_data_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_data_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_data_cut_all = ROOT.TH1D("P_gtr_xp_pions_data_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_data_cut_all = ROOT.TH1D("P_gtr_yp_pions_data_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_data_cut_all = ROOT.TH1D("P_gtr_dp_pions_data_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_data_cut_all = ROOT.TH1D("P_gtr_p_pions_data_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_data_cut_all = ROOT.TH1D("P_dc_x_fp_pions_data_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_data_cut_all = ROOT.TH1D("P_dc_y_fp_pions_data_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_data_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_data_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_data_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_data_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_data_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_data_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_data_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_data_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_data_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_data_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_data_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_data_cut_all = ROOT.TH1D("MMpi_pions_data_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_data_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_data_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_data_cut_all = ROOT.TH1D("pmiss_pions_data_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_data_cut_all = ROOT.TH1D("pmiss_x_pions_data_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_data_cut_all = ROOT.TH1D("pmiss_y_pions_data_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_data_cut_all = ROOT.TH1D("pmiss_z_pions_data_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_data_cut_all = ROOT.TH1D("emiss_pions_data_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_data_cut_all = ROOT.TH1D("W_pions_data_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_data_cut_all = ROOT.TH1D("Q2_pions_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_data_cut_all = ROOT.TH1D("Epsilon_pions_data_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_data_cut_all = ROOT.TH1D("t_pions_data_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_data_cut_all = ROOT.TH1D("ph_q_pions_data_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)
'''
# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Data
H_gtr_beta_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_prompt_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_xp_pions_data_prompt_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_yp_pions_data_prompt_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_dp_pions_data_prompt_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_p_pions_data_prompt_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_data_prompt_cut_all = ROOT.TH1D("H_dc_x_fp_pions_data_prompt_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_data_prompt_cut_all = ROOT.TH1D("H_dc_y_fp_pions_data_prompt_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_data_prompt_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_data_prompt_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_data_prompt_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_data_prompt_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_data_prompt_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_data_prompt_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_data_prompt_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_prompt_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_prompt_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_prompt_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_prompt_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_prompt_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_xp_pions_data_prompt_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_yp_pions_data_prompt_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_dp_pions_data_prompt_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_p_pions_data_prompt_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_data_prompt_cut_all = ROOT.TH1D("P_dc_x_fp_pions_data_prompt_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_data_prompt_cut_all = ROOT.TH1D("P_dc_y_fp_pions_data_prompt_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_data_prompt_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_data_prompt_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_data_prompt_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_data_prompt_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_data_prompt_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_data_prompt_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_data_prompt_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_prompt_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_prompt_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_prompt_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_prompt_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_prompt_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_prompt_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_prompt_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_prompt_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_data_prompt_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_prompt_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_prompt_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_data_prompt_cut_all = ROOT.TH1D("MMpi_pions_data_prompt_cut_all", "MIssing Mass data (prompt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_data_prompt_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_prompt_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_data_prompt_cut_all = ROOT.TH1D("pmiss_pions_data_prompt_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_data_prompt_cut_all = ROOT.TH1D("pmiss_x_pions_data_prompt_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_data_prompt_cut_all = ROOT.TH1D("pmiss_y_pions_data_prompt_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_data_prompt_cut_all = ROOT.TH1D("pmiss_z_pions_data_prompt_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_data_prompt_cut_all = ROOT.TH1D("emiss_pions_data_prompt_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_data_prompt_cut_all = ROOT.TH1D("W_pions_data_prompt_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_data_prompt_cut_all = ROOT.TH1D("Q2_pions_data_prompt_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_data_prompt_cut_all = ROOT.TH1D("Epsilon_pions_data_prompt_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_data_prompt_cut_all = ROOT.TH1D("t_pions_data_prompt_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_data_prompt_cut_all = ROOT.TH1D("ph_q_pions_data_prompt_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Data
H_gtr_beta_pions_data_random_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_random_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_data_random_cut_all = ROOT.TH1D("H_gtr_xp_pions_data_random_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_data_random_cut_all = ROOT.TH1D("H_gtr_yp_pions_data_random_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_data_random_cut_all = ROOT.TH1D("H_gtr_dp_pions_data_random_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_data_random_cut_all = ROOT.TH1D("H_gtr_p_pions_data_random_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_data_random_cut_all = ROOT.TH1D("H_dc_x_fp_pions_data_random_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_data_random_cut_all = ROOT.TH1D("H_dc_y_fp_pions_data_random_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_data_random_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_data_random_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_data_random_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_data_random_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_data_random_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_data_random_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_data_random_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_random_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_data_random_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_random_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_random_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_random_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_random_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_random_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_data_random_cut_all = ROOT.TH1D("P_gtr_xp_pions_data_random_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_data_random_cut_all = ROOT.TH1D("P_gtr_yp_pions_data_random_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_data_random_cut_all = ROOT.TH1D("P_gtr_dp_pions_data_random_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_data_random_cut_all = ROOT.TH1D("P_gtr_p_pions_data_random_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_data_random_cut_all = ROOT.TH1D("P_dc_x_fp_pions_data_random_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_data_random_cut_all = ROOT.TH1D("P_dc_y_fp_pions_data_random_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_data_random_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_data_random_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_data_random_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_data_random_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_data_random_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_data_random_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_data_random_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_random_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_random_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_random_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_random_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_random_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_random_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_random_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_random_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_data_random_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_random_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_random_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_data_random_cut_all = ROOT.TH1D("MMpi_pions_data_random_cut_all", "MIssing Mass data (random_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_data_random_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_random_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_data_random_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_random_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_data_random_cut_all = ROOT.TH1D("pmiss_pions_data_random_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_data_random_cut_all = ROOT.TH1D("pmiss_x_pions_data_random_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_data_random_cut_all = ROOT.TH1D("pmiss_y_pions_data_random_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_data_random_cut_all = ROOT.TH1D("pmiss_z_pions_data_random_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_data_random_cut_all = ROOT.TH1D("emiss_pions_data_random_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_data_random_cut_all = ROOT.TH1D("W_pions_data_random_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_data_random_cut_all = ROOT.TH1D("Q2_pions_data_random_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_data_random_cut_all = ROOT.TH1D("Epsilon_pions_data_random_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_data_random_cut_all = ROOT.TH1D("t_pions_data_random_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_data_random_cut_all = ROOT.TH1D("ph_q_pions_data_random_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Dummy
H_gtr_beta_pions_dummy_prompt_cut_all = ROOT.TH1D("H_gtr_beta_pions_dummy_prompt_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_gtr_xp_pions_dummy_prompt_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_gtr_yp_pions_dummy_prompt_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_gtr_dp_pions_dummy_prompt_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_dummy_prompt_cut_all = ROOT.TH1D("H_gtr_p_pions_dummy_prompt_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_dc_x_fp_pions_dummy_prompt_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_dc_y_fp_pions_dummy_prompt_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_dummy_prompt_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_dummy_prompt_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_dummy_prompt_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_dummy_prompt_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_dummy_prompt_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_dummy_prompt_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_dummy_prompt_cut_all = ROOT.TH1D("H_cer_npeSum_pions_dummy_prompt_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_dummy_prompt_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_dummy_prompt_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_dummy_prompt_cut_all = ROOT.TH1D("P_gtr_beta_pions_dummy_prompt_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_gtr_xp_pions_dummy_prompt_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_gtr_yp_pions_dummy_prompt_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_gtr_dp_pions_dummy_prompt_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_dummy_prompt_cut_all = ROOT.TH1D("P_gtr_p_pions_dummy_prompt_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_dc_x_fp_pions_dummy_prompt_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_dc_y_fp_pions_dummy_prompt_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_dummy_prompt_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_dummy_prompt_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_dummy_prompt_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_dummy_prompt_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_dummy_prompt_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_dummy_prompt_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_dummy_prompt_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_dummy_prompt_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_dummy_prompt_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_dummy_prompt_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_dummy_prompt_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_dummy_prompt_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_dummy_prompt_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_dummy_prompt_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_dummy_prompt_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_dummy_prompt_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_dummy_prompt_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_dummy_prompt_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_dummy_prompt_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_dummy_prompt_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_prompt_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_dummy_prompt_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_dummy_prompt_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_dummy_prompt_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_dummy_prompt_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_dummy_prompt_cut_all = ROOT.TH1D("MMpi_pions_dummy_prompt_cut_all", "MIssing Mass dummy (prompt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_dummy_prompt_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_dummy_prompt_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_dummy_prompt_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_prompt_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_dummy_prompt_cut_all = ROOT.TH1D("pmiss_pions_dummy_prompt_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_dummy_prompt_cut_all = ROOT.TH1D("pmiss_x_pions_dummy_prompt_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_dummy_prompt_cut_all = ROOT.TH1D("pmiss_y_pions_dummy_prompt_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_dummy_prompt_cut_all = ROOT.TH1D("pmiss_z_pions_dummy_prompt_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_dummy_prompt_cut_all = ROOT.TH1D("emiss_pions_dummy_prompt_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_dummy_prompt_cut_all = ROOT.TH1D("W_pions_dummy_prompt_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_dummy_prompt_cut_all = ROOT.TH1D("Q2_pions_dummy_prompt_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_dummy_prompt_cut_all = ROOT.TH1D("Epsilon_pions_dummy_prompt_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_dummy_prompt_cut_all = ROOT.TH1D("t_pions_dummy_prompt_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_dummy_prompt_cut_all = ROOT.TH1D("ph_q_pions_dummy_prompt_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Dummy
H_gtr_beta_pions_dummy_random_cut_all = ROOT.TH1D("H_gtr_beta_pions_dummy_random_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_dummy_random_cut_all = ROOT.TH1D("H_gtr_xp_pions_dummy_random_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_dummy_random_cut_all = ROOT.TH1D("H_gtr_yp_pions_dummy_random_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_dummy_random_cut_all = ROOT.TH1D("H_gtr_dp_pions_dummy_random_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_dummy_random_cut_all = ROOT.TH1D("H_gtr_p_pions_dummy_random_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_dummy_random_cut_all = ROOT.TH1D("H_dc_x_fp_pions_dummy_random_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_dummy_random_cut_all = ROOT.TH1D("H_dc_y_fp_pions_dummy_random_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_dummy_random_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_dummy_random_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_dummy_random_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_dummy_random_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_dummy_random_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_dummy_random_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_dummy_random_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_dummy_random_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_dummy_random_cut_all = ROOT.TH1D("H_cer_npeSum_pions_dummy_random_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_dummy_random_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_dummy_random_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_dummy_random_cut_all = ROOT.TH1D("P_gtr_beta_pions_dummy_random_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_dummy_random_cut_all = ROOT.TH1D("P_gtr_xp_pions_dummy_random_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_dummy_random_cut_all = ROOT.TH1D("P_gtr_yp_pions_dummy_random_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_dummy_random_cut_all = ROOT.TH1D("P_gtr_dp_pions_dummy_random_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_dummy_random_cut_all = ROOT.TH1D("P_gtr_p_pions_dummy_random_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_dummy_random_cut_all = ROOT.TH1D("P_dc_x_fp_pions_dummy_random_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_dummy_random_cut_all = ROOT.TH1D("P_dc_y_fp_pions_dummy_random_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_dummy_random_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_dummy_random_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_dummy_random_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_dummy_random_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_dummy_random_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_dummy_random_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_dummy_random_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_dummy_random_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_dummy_random_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_dummy_random_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_dummy_random_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_dummy_random_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_dummy_random_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_dummy_random_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_dummy_random_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_dummy_random_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_dummy_random_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_dummy_random_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_dummy_random_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_dummy_random_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_dummy_random_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_random_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_dummy_random_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_dummy_random_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_dummy_random_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_dummy_random_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_dummy_random_cut_all = ROOT.TH1D("MMpi_pions_dummy_random_cut_all", "MIssing Mass dummy (random_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_dummy_random_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_dummy_random_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_dummy_random_cut_all = ROOT.TH1D("pmiss_pions_dummy_random_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_dummy_random_cut_all = ROOT.TH1D("pmiss_x_pions_dummy_random_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_dummy_random_cut_all = ROOT.TH1D("pmiss_y_pions_dummy_random_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_dummy_random_cut_all = ROOT.TH1D("pmiss_z_pions_dummy_random_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_dummy_random_cut_all = ROOT.TH1D("emiss_pions_dummy_random_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_dummy_random_cut_all = ROOT.TH1D("W_pions_dummy_random_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_dummy_random_cut_all = ROOT.TH1D("Q2_pions_dummy_random_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_dummy_random_cut_all = ROOT.TH1D("Epsilon_pions_dummy_random_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_dummy_random_cut_all = ROOT.TH1D("t_pions_dummy_random_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_dummy_random_cut_all = ROOT.TH1D("ph_q_pions_dummy_random_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Data
H_gtr_beta_pions_randsub_data_cut_all = ROOT.TH1D("H_gtr_beta_pions_randsub_data_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_randsub_data_cut_all = ROOT.TH1D("H_gtr_xp_pions_randsub_data_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_randsub_data_cut_all = ROOT.TH1D("H_gtr_yp_pions_randsub_data_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_randsub_data_cut_all = ROOT.TH1D("H_gtr_dp_pions_randsub_data_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_randsub_data_cut_all = ROOT.TH1D("H_gtr_p_pions_randsub_data_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_randsub_data_cut_all = ROOT.TH1D("H_dc_x_fp_pions_randsub_data_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_randsub_data_cut_all = ROOT.TH1D("H_dc_y_fp_pions_randsub_data_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_randsub_data_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_randsub_data_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_randsub_data_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_randsub_data_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_randsub_data_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_randsub_data_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_randsub_data_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_randsub_data_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_randsub_data_cut_all = ROOT.TH1D("H_cer_npeSum_pions_randsub_data_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_randsub_data_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_randsub_data_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_randsub_data_cut_all = ROOT.TH1D("P_gtr_beta_pions_randsub_data_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_randsub_data_cut_all = ROOT.TH1D("P_gtr_xp_pions_randsub_data_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_randsub_data_cut_all = ROOT.TH1D("P_gtr_yp_pions_randsub_data_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_randsub_data_cut_all = ROOT.TH1D("P_gtr_dp_pions_randsub_data_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_randsub_data_cut_all = ROOT.TH1D("P_gtr_p_pions_randsub_data_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_randsub_data_cut_all = ROOT.TH1D("P_dc_x_fp_pions_randsub_data_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_randsub_data_cut_all = ROOT.TH1D("P_dc_y_fp_pions_randsub_data_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_randsub_data_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_randsub_data_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_randsub_data_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_randsub_data_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_randsub_data_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_randsub_data_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_randsub_data_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_randsub_data_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_randsub_data_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_randsub_data_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_randsub_data_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_randsub_data_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_randsub_data_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_randsub_data_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_randsub_data_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_randsub_data_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_randsub_data_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_randsub_data_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_randsub_data_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_randsub_data_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_randsub_data_cut_all = ROOT.TH1D("P_aero_npeSum_pions_randsub_data_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_randsub_data_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_randsub_data_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_randsub_data_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_randsub_data_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_randsub_data_cut_all = ROOT.TH1D("MMpi_pions_randsub_data_cut_all", "MIssing Mass data (randsub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_randsub_data_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_randsub_data_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_randsub_data_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_randsub_data_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_randsub_data_cut_all = ROOT.TH1D("pmiss_pions_randsub_data_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_randsub_data_cut_all = ROOT.TH1D("pmiss_x_pions_randsub_data_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_randsub_data_cut_all = ROOT.TH1D("pmiss_y_pions_randsub_data_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_randsub_data_cut_all = ROOT.TH1D("pmiss_z_pions_randsub_data_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_randsub_data_cut_all = ROOT.TH1D("emiss_pions_randsub_data_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_randsub_data_cut_all = ROOT.TH1D("W_pions_randsub_data_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_randsub_data_cut_all = ROOT.TH1D("Q2_pions_randsub_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_randsub_data_cut_all = ROOT.TH1D("Epsilon_pions_randsub_data_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_randsub_data_cut_all = ROOT.TH1D("t_pions_randsub_data_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_randsub_data_cut_all = ROOT.TH1D("ph_q_pions_randsub_data_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Dummy
H_gtr_beta_pions_randsub_dummy_cut_all = ROOT.TH1D("H_gtr_beta_pions_randsub_dummy_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_gtr_xp_pions_randsub_dummy_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_gtr_yp_pions_randsub_dummy_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_gtr_dp_pions_randsub_dummy_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_randsub_dummy_cut_all = ROOT.TH1D("H_gtr_p_pions_randsub_dummy_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_dc_x_fp_pions_randsub_dummy_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_dc_y_fp_pions_randsub_dummy_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_randsub_dummy_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_randsub_dummy_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_randsub_dummy_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_randsub_dummy_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_randsub_dummy_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_randsub_dummy_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_randsub_dummy_cut_all = ROOT.TH1D("H_cer_npeSum_pions_randsub_dummy_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_randsub_dummy_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_randsub_dummy_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_randsub_dummy_cut_all = ROOT.TH1D("P_gtr_beta_pions_randsub_dummy_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_gtr_xp_pions_randsub_dummy_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_gtr_yp_pions_randsub_dummy_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_gtr_dp_pions_randsub_dummy_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_randsub_dummy_cut_all = ROOT.TH1D("P_gtr_p_pions_randsub_dummy_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_dc_x_fp_pions_randsub_dummy_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_dc_y_fp_pions_randsub_dummy_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_randsub_dummy_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_randsub_dummy_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_randsub_dummy_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_randsub_dummy_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_randsub_dummy_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_randsub_dummy_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_randsub_dummy_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_randsub_dummy_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_randsub_dummy_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_randsub_dummy_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_randsub_dummy_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_randsub_dummy_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_randsub_dummy_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_randsub_dummy_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_randsub_dummy_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_randsub_dummy_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_randsub_dummy_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_randsub_dummy_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_randsub_dummy_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_randsub_dummy_cut_all = ROOT.TH1D("P_aero_npeSum_pions_randsub_dummy_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_randsub_dummy_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_randsub_dummy_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_randsub_dummy_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_randsub_dummy_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_randsub_dummy_cut_all = ROOT.TH1D("MMpi_pions_randsub_dummy_cut_all", "MIssing Mass dummy (randsub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_randsub_dummy_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_randsub_dummy_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_randsub_dummy_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_randsub_dummy_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_randsub_dummy_cut_all = ROOT.TH1D("pmiss_pions_randsub_dummy_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_randsub_dummy_cut_all = ROOT.TH1D("pmiss_x_pions_randsub_dummy_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_randsub_dummy_cut_all = ROOT.TH1D("pmiss_y_pions_randsub_dummy_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_randsub_dummy_cut_all = ROOT.TH1D("pmiss_z_pions_randsub_dummy_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_randsub_dummy_cut_all = ROOT.TH1D("emiss_pions_randsub_dummy_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_randsub_dummy_cut_all = ROOT.TH1D("W_pions_randsub_dummy_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_randsub_dummy_cut_all = ROOT.TH1D("Q2_pions_randsub_dummy_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_randsub_dummy_cut_all = ROOT.TH1D("Epsilon_pions_randsub_dummy_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_randsub_dummy_cut_all = ROOT.TH1D("t_pions_randsub_dummy_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_randsub_dummy_cut_all = ROOT.TH1D("ph_q_pions_randsub_dummy_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + Norm) Data
H_gtr_beta_pions_normrandsub_data_cut_all = ROOT.TH1D("H_gtr_beta_pions_normrandsub_data_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_gtr_xp_pions_normrandsub_data_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_gtr_yp_pions_normrandsub_data_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_gtr_dp_pions_normrandsub_data_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_normrandsub_data_cut_all = ROOT.TH1D("H_gtr_p_pions_normrandsub_data_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_dc_x_fp_pions_normrandsub_data_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_dc_y_fp_pions_normrandsub_data_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_normrandsub_data_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_normrandsub_data_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_normrandsub_data_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_normrandsub_data_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_normrandsub_data_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_normrandsub_data_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_normrandsub_data_cut_all = ROOT.TH1D("H_cer_npeSum_pions_normrandsub_data_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_normrandsub_data_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_normrandsub_data_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_normrandsub_data_cut_all = ROOT.TH1D("P_gtr_beta_pions_normrandsub_data_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_gtr_xp_pions_normrandsub_data_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_gtr_yp_pions_normrandsub_data_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_gtr_dp_pions_normrandsub_data_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_normrandsub_data_cut_all = ROOT.TH1D("P_gtr_p_pions_normrandsub_data_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_dc_x_fp_pions_normrandsub_data_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_dc_y_fp_pions_normrandsub_data_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_normrandsub_data_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_normrandsub_data_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_normrandsub_data_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_normrandsub_data_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_normrandsub_data_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_normrandsub_data_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_normrandsub_data_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_normrandsub_data_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_normrandsub_data_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_normrandsub_data_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_normrandsub_data_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_normrandsub_data_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_normrandsub_data_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_normrandsub_data_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_normrandsub_data_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_normrandsub_data_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_normrandsub_data_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_normrandsub_data_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_normrandsub_data_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_normrandsub_data_cut_all = ROOT.TH1D("P_aero_npeSum_pions_normrandsub_data_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_normrandsub_data_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_normrandsub_data_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_normrandsub_data_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_normrandsub_data_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_normrandsub_data_cut_all = ROOT.TH1D("MMpi_pions_normrandsub_data_cut_all", "MIssing Mass data (datasub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_normrandsub_data_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_normrandsub_data_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_normrandsub_data_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_normrandsub_data_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_normrandsub_data_cut_all = ROOT.TH1D("pmiss_pions_normrandsub_data_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_normrandsub_data_cut_all = ROOT.TH1D("pmiss_x_pions_normrandsub_data_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_normrandsub_data_cut_all = ROOT.TH1D("pmiss_y_pions_normrandsub_data_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_normrandsub_data_cut_all = ROOT.TH1D("pmiss_z_pions_normrandsub_data_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_normrandsub_data_cut_all = ROOT.TH1D("emiss_pions_normrandsub_data_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_normrandsub_data_cut_all = ROOT.TH1D("W_pions_normrandsub_data_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_normrandsub_data_cut_all = ROOT.TH1D("Q2_pions_normrandsub_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_normrandsub_data_cut_all = ROOT.TH1D("Epsilon_pions_normrandsub_data_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_normrandsub_data_cut_all = ROOT.TH1D("t_pions_normrandsub_data_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_normrandsub_data_cut_all = ROOT.TH1D("ph_q_pions_normrandsub_data_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + Norm) Dummy
H_gtr_beta_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_gtr_beta_pions_normrandsub_dummy_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_gtr_xp_pions_normrandsub_dummy_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_gtr_yp_pions_normrandsub_dummy_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_gtr_dp_pions_normrandsub_dummy_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_gtr_p_pions_normrandsub_dummy_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_dc_x_fp_pions_normrandsub_dummy_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_dc_y_fp_pions_normrandsub_dummy_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_normrandsub_dummy_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_normrandsub_dummy_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_normrandsub_dummy_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_normrandsub_dummy_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_cer_npeSum_pions_normrandsub_dummy_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_normrandsub_dummy_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_normrandsub_dummy_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_gtr_beta_pions_normrandsub_dummy_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_gtr_xp_pions_normrandsub_dummy_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_gtr_yp_pions_normrandsub_dummy_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_gtr_dp_pions_normrandsub_dummy_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_gtr_p_pions_normrandsub_dummy_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_dc_x_fp_pions_normrandsub_dummy_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_dc_y_fp_pions_normrandsub_dummy_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_normrandsub_dummy_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_normrandsub_dummy_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_normrandsub_dummy_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_normrandsub_dummy_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_normrandsub_dummy_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_normrandsub_dummy_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_normrandsub_dummy_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_normrandsub_dummy_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_normrandsub_dummy_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_normrandsub_dummy_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_aero_npeSum_pions_normrandsub_dummy_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_normrandsub_dummy_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_normrandsub_dummy_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_normrandsub_dummy_cut_all = ROOT.TH1D("MMpi_pions_normrandsub_dummy_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_normrandsub_dummy_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_normrandsub_dummy_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_normrandsub_dummy_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_normrandsub_dummy_cut_all = ROOT.TH1D("pmiss_pions_normrandsub_dummy_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_normrandsub_dummy_cut_all = ROOT.TH1D("pmiss_x_pions_normrandsub_dummy_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_normrandsub_dummy_cut_all = ROOT.TH1D("pmiss_y_pions_normrandsub_dummy_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_normrandsub_dummy_cut_all = ROOT.TH1D("pmiss_z_pions_normrandsub_dummy_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_normrandsub_dummy_cut_all = ROOT.TH1D("emiss_pions_normrandsub_dummy_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_normrandsub_dummy_cut_all = ROOT.TH1D("W_pions_normrandsub_dummy_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_normrandsub_dummy_cut_all = ROOT.TH1D("Q2_pions_normrandsub_dummy_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_normrandsub_dummy_cut_all = ROOT.TH1D("Epsilon_pions_normrandsub_dummy_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_normrandsub_dummy_cut_all = ROOT.TH1D("t_pions_normrandsub_dummy_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_normrandsub_dummy_cut_all = ROOT.TH1D("ph_q_pions_normrandsub_dummy_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# Histograms having Cuts (Acceptance + PID + RF + Prompt Selection) Norm Dummy Subtraction Data
H_gtr_beta_pions_normdummysub_data_cut_all = ROOT.TH1D("H_gtr_beta_pions_normdummysub_data_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_gtr_xp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_gtr_xp_pions_normdummysub_data_cut_all", "HMS xptar; HMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
H_gtr_yp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_gtr_yp_pions_normdummysub_data_cut_all", "HMS yptar; HMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
H_gtr_dp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_gtr_dp_pions_normdummysub_data_cut_all", "HMS #delta; HMS_gtr_dp; Counts", nbins, -15, 15)
H_gtr_p_pions_normdummysub_data_cut_all = ROOT.TH1D("H_gtr_p_pions_normdummysub_data_cut_all", "HMS p; HMS_gtr_p; Counts", nbins, 4, 8)
H_dc_x_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_dc_x_fp_pions_normdummysub_data_cut_all", "HMS x_fp'; HMS_dc_x_fp; Counts", nbins, -100, 100)
H_dc_y_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_dc_y_fp_pions_normdummysub_data_cut_all", "HMS y_fp'; HMS_dc_y_fp; Counts", nbins, -100, 100)
H_dc_xp_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_dc_xp_fp_pions_normdummysub_data_cut_all", "HMS xp_fp'; HMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
H_dc_yp_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("H_dc_yp_fp_pions_normdummysub_data_cut_all", "HMS yp_fp'; HMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
H_cal_etotnorm_pions_normdummysub_data_cut_all = ROOT.TH1D("H_cal_etotnorm_pions_normdummysub_data_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", nbins, 0.2, 1.8)
H_cal_etottracknorm_pions_normdummysub_data_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_normdummysub_data_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", nbins, 0.2, 1.8)
H_cer_npeSum_pions_normdummysub_data_cut_all = ROOT.TH1D("H_cer_npeSum_pions_normdummysub_data_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_normdummysub_data_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_normdummysub_data_cut_all", "HMS RFTime; HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_normdummysub_data_cut_all = ROOT.TH1D("P_gtr_beta_pions_normdummysub_data_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", nbins, 0.5, 1.3)
P_gtr_xp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_gtr_xp_pions_normdummysub_data_cut_all", "SHMS xptar; SHMS_gtr_xptar; Counts", nbins, -0.2, 0.2)
P_gtr_yp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_gtr_yp_pions_normdummysub_data_cut_all", "SHMS yptar; SHMS_gtr_yptar; Counts", nbins, -0.2, 0.2)
P_gtr_dp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_gtr_dp_pions_normdummysub_data_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", nbins, -30, 30)
P_gtr_p_pions_normdummysub_data_cut_all = ROOT.TH1D("P_gtr_p_pions_normdummysub_data_cut_all", "SHMS p; SHMS_gtr_p; Counts", nbins, 1, 7)
P_dc_x_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_dc_x_fp_pions_normdummysub_data_cut_all", "SHMS x_fp'; SHMS_dc_x_fp; Counts", nbins, -100, 100)
P_dc_y_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_dc_y_fp_pions_normdummysub_data_cut_all", "SHMS y_fp'; SHMS_dc_y_fp; Counts", nbins, -100, 100)
P_dc_xp_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_dc_xp_fp_pions_normdummysub_data_cut_all", "SHMS xp_fp'; SHMS_dc_xp_fp; Counts", nbins, -0.2, 0.2)
P_dc_yp_fp_pions_normdummysub_data_cut_all = ROOT.TH1D("P_dc_yp_fp_pions_normdummysub_data_cut_all", "SHMS yp_fp'; SHMS_dc_yp_fp; Counts", nbins, -0.2, 0.2)
P_cal_etotnorm_pions_normdummysub_data_cut_all = ROOT.TH1D("P_cal_etotnorm_pions_normdummysub_data_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", nbins, 0, 1)
P_cal_etottracknorm_pions_normdummysub_data_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_normdummysub_data_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", nbins, 0, 1.6)
P_hgcer_npeSum_pions_normdummysub_data_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_normdummysub_data_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_normdummysub_data_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_normdummysub_data_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_normdummysub_data_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_normdummysub_data_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", nbins, -50, 50)
P_ngcer_npeSum_pions_normdummysub_data_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_normdummysub_data_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_normdummysub_data_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_normdummysub_data_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_normdummysub_data_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_normdummysub_data_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", nbins, -50, 50)
P_aero_npeSum_pions_normdummysub_data_cut_all = ROOT.TH1D("P_aero_npeSum_pions_normdummysub_data_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_normdummysub_data_cut_all = ROOT.TH1D("P_aero_xAtAero_pions_normdummysub_data_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_normdummysub_data_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_normdummysub_data_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", nbins, -50, 50)
MMpi_pions_normdummysub_data_cut_all = ROOT.TH1D("MMpi_pions_normdummysub_data_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)
P_RFTime_Dist_pions_normdummysub_data_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_normdummysub_data_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC2_pions_normdummysub_data_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_normdummysub_data_cut_all", "Electron-Pion CTime; e p Coin_Time; Counts", nbins, -50, 50)
pmiss_pions_normdummysub_data_cut_all = ROOT.TH1D("pmiss_pions_normdummysub_data_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_normdummysub_data_cut_all = ROOT.TH1D("pmiss_x_pions_normdummysub_data_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_normdummysub_data_cut_all = ROOT.TH1D("pmiss_y_pions_normdummysub_data_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_normdummysub_data_cut_all = ROOT.TH1D("pmiss_z_pions_normdummysub_data_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_normdummysub_data_cut_all = ROOT.TH1D("emiss_pions_normdummysub_data_cut_all", "Energy Distribution; emiss; Counts", nbins, -0.5, 2.5)
W_pions_normdummysub_data_cut_all = ROOT.TH1D("W_pions_normdummysub_data_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
Q2_pions_normdummysub_data_cut_all = ROOT.TH1D("Q2_pions_normdummysub_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Epsilon_pions_normdummysub_data_cut_all = ROOT.TH1D("Epsilon_pions_normdummysub_data_cut_all", "Epsilon Distribution; \epsilon; Counts", nbins, 0, 2.0)
t_pions_normdummysub_data_cut_all = ROOT.TH1D("t_pions_normdummysub_data_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
ph_q_pions_normdummysub_data_cut_all = ROOT.TH1D("ph_q_pions_normdummysub_data_cut_all", "phi Distribution; \phi; Counts", nbins, 0, 360)

# SIMC Histograms with Cuts
H_hsdelta_pions_simc_cut_all = ROOT.TH1D("H_hsdelta_pions_simc_cut_all", "HMS #delta; HMS_#delta; Counts", nbins, -15, 15)
H_hsxptar_pions_simc_cut_all = ROOT.TH1D("H_hsxptar_pions_simc_cut_all", "HMS xptar; HMS_xptar; Counts", nbins, -0.2, 0.2)
H_hsyptar_pions_simc_cut_all = ROOT.TH1D("H_hsyptar_pions_simc_cut_all", "HMS yptar; HMS_yptar; Counts", nbins, -0.2, 0.2)
H_hsytar_pions_simc_cut_all = ROOT.TH1D("H_hsytar_pions_simc_cut_all", "HMS ytar; HMS_ytar; Counts", nbins, -20, 20)
H_hsxfp_pions_simc_cut_all = ROOT.TH1D("H_hsxfp_pions_simc_cut_all", "HMS x_fp'; HMS_xfp; Counts", nbins, -100, 100)
H_hsyfp_pions_simc_cut_all = ROOT.TH1D("H_hsyfp_pions_simc_cut_all", "HMS y_fp'; HMS_yfp; Counts", nbins, -100, 100)
H_hsxpfp_pions_simc_cut_all = ROOT.TH1D("H_hsxpfp_pions_simc_cut_all", "HMS xp_fp'; HMS_xpfp; Counts", nbins, -0.2, 0.2)
H_hsypfp_pions_simc_cut_all = ROOT.TH1D("H_hsypfp_pions_simc_cut_all", "HMS yp_fp'; HMS_ypfp; Counts", nbins, -0.2, 0.2)
P_ssdelta_pions_simc_cut_all = ROOT.TH1D("P_ssdelta_pions_simc_cut_all", "SHMS #delta; SHMS_#delta; Counts", nbins, -30, 30)
P_ssxptar_pions_simc_cut_all = ROOT.TH1D("P_ssxptar_pions_simc_cut_all", "SHMS xptar; SHMS_xptar; Counts", nbins, -0.2, 0.2)
P_ssyptar_pions_simc_cut_all = ROOT.TH1D("P_ssyptar_pions_simc_cut_all", "SHMS yptar; SHMS_yptar; Counts", nbins, -0.2, 0.2)
P_ssytar_pions_simc_cut_all = ROOT.TH1D("P_ssytar_pions_simc_cut_all", "SHMS ytar; SHMS_ytar; Counts", nbins, -20, 20)
P_ssxfp_pions_simc_cut_all = ROOT.TH1D("P_ssxfp_pions_simc_cut_all", "SHMS x_fp'; SHMS_xfp; Counts", nbins, -100, 100)
P_ssyfp_pions_simc_cut_all = ROOT.TH1D("P_ssyfp_pions_simc_cut_all", "SHMS y_fp'; SHMS_yfp; Counts", nbins, -100, 100)
P_ssxpfp_pions_simc_cut_all = ROOT.TH1D("P_ssxpfp_pions_simc_cut_all", "SHMS xp_fp'; SHMS_xpfp; Counts", nbins, -0.2, 0.2)
P_ssypfp_pions_simc_cut_all = ROOT.TH1D("P_ssypfp_pions_simc_cut_all", "SHMS yp_fp'; SHMS_ypfp; Counts", nbins, -0.2, 0.2)
q_pions_simc_cut_all = ROOT.TH1D("q_pions_simc_cut_all", "q Distribution; q; Counts", nbins, 0, 10)
nu_pions_simc_cut_all = ROOT.TH1D("nu_pions_simc_cut_all", "nu Distribution; nu; Counts", nbins, 0, 10)
Q2_pions_simc_cut_all = ROOT.TH1D("Q2_pions_simc_cut_all", "#Q^2 Distribution; Q^{2}; Counts", nbins, 0, 10)
epsilon_pions_simc_cut_all = ROOT.TH1D("epsilon_pions_simc_cut_all", "epsilon Distribution; #epsilon; Counts", nbins, 0, 2)
thetapq_pions_simc_cut_all = ROOT.TH1D("thetapq_pions_simc_cut_all", "thetapq Distribution; thetapq; Counts", nbins, -1, 1)
phipq_pions_simc_cut_all = ROOT.TH1D("phipq_pions_simc_cut_all", "phipq Distribution; phipq; Counts", nbins, 0, 360)
pmiss_pions_simc_cut_all = ROOT.TH1D("pmiss_pions_simc_cut_all", "Momentum Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_x_pions_simc_cut_all = ROOT.TH1D("pmiss_x_pions_simc_cut_all", "Momentum_x Distribution; pmiss_x; Counts", nbins, -0.8, 0.8)
pmiss_y_pions_simc_cut_all = ROOT.TH1D("pmiss_y_pions_simc_cut_all", "Momentum_y Distribution; pmiss_y; Counts", nbins, -0.8, 0.8)
pmiss_z_pions_simc_cut_all = ROOT.TH1D("pmiss_z_pions_simc_cut_all", "Momentum_z Distribution; pmiss_z; Counts", nbins, -2.0, 2.0)
emiss_pions_simc_cut_all = ROOT.TH1D("emiss_pions_simc_cut_all", "Energy Distribution; emiss; Counts",nbins, -0.5, 2.5)
W_pions_simc_cut_all = ROOT.TH1D("W_pions_simc_cut_all", "W Distribution; W; Counts", nbins_p, 0, 10.0)
t_pions_simc_cut_all = ROOT.TH1D("t_pions_simc_cut_all", "t Distribution; t; Counts", nbins, -2.0, 2.0)
MMpi_pions_simc_cut_all = ROOT.TH1D("MMpi_pions_simc_cut_all", "MIssing Mass SIMC (cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms for Error Calculations
pmiss_pions_data_prompt_cut_all_error = ROOT.TH1D("pmiss_pions_data_prompt_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0) 
pmiss_pions_dummy_prompt_cut_all_error = ROOT.TH1D("pmiss_pions_dummy_prompt_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_pions_data_random_cut_all_error = ROOT.TH1D("pmiss_pions_data_random_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0) 
pmiss_pions_dummy_random_cut_all_error = ROOT.TH1D("pmiss_pions_dummy_random_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_pions_randsub_data_cut_all_error = ROOT.TH1D("pmiss_pions_randsub_data_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_pions_randsub_dummy_cut_all_error = ROOT.TH1D("pmiss_pions_randsub_dummy_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0)
pmiss_pions_simc_cut_all_error = ROOT.TH1D("pmiss_pions_simc_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins, 0.0, 2.0)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")
infile_SIMC = ROOT.TFile.Open(rootFile_SIMC, "READ")

#Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_All_Data_tree = infile_DATA.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DATA  = Cut_Pion_Events_Prompt_Data_tree.GetEntries()

#Uncut_Pion_Events_Dummy_tree = infile_DUMMY.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Accpt")
#Cut_Pion_Events_All_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DUMMY  = Cut_Pion_Events_Prompt_Dummy_tree.GetEntries()

Uncut_Pion_Events_SIMC_tree = infile_SIMC.Get("h10")
#nEntries_TBRANCH_SIMC  = Uncut_Proton_Events_SIMC_tree.GetEntries()

###################################################################################################################################################

#Fill histograms for Cut All Data

for event in Cut_Pion_Events_Accpt_Data_tree:
#    if event.MMpi<0.10:
        H_gtr_beta_pions_data_accpt_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_data_accpt_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_data_accpt_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_data_accpt_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_data_accpt_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_data_accpt_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_data_accpt_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_data_accpt_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_data_accpt_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_data_accpt_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_data_accpt_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_accpt_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_accpt_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_accpt_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_data_accpt_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_data_accpt_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_data_accpt_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_data_accpt_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_data_accpt_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_data_accpt_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_data_accpt_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_data_accpt_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_data_accpt_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_data_accpt_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_accpt_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_accpt_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_accpt_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_accpt_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_accpt_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_accpt_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_accpt_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_accpt_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_accpt_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_data_accpt_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_data_accpt_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_data_accpt_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_data_accpt_cut_all.Fill(event.pmiss)
        pmiss_x_pions_data_accpt_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_data_accpt_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_data_accpt_cut_all.Fill(event.pmiss_z)
        emiss_pions_data_accpt_cut_all.Fill(event.emiss)
        W_pions_data_accpt_cut_all.Fill(event.W)
        Q2_pions_data_accpt_cut_all.Fill(event.Q2)
        Epsilon_pions_data_accpt_cut_all.Fill(event.epsilon)
        t_pions_data_accpt_cut_all.Fill(-event.MandelT)
        ph_q_pions_data_accpt_cut_all.Fill(event.ph_q)
#    ibin += 1
'''
for event in Cut_Pion_Events_All_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        H_gtr_beta_pions_data_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_data_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_data_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_data_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_data_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_data_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_data_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_data_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_data_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_data_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_data_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_data_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_data_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_data_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_data_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_data_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_data_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_data_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_data_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_data_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_data_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_data_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_data_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_data_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_data_cut_all.Fill(event.pmiss)
        pmiss_x_pions_data_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_data_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_data_cut_all.Fill(event.pmiss_z)
        emiss_pions_data_cut_all.Fill(event.emiss)
        W_pions_data_cut_all.Fill(event.W)
        Q2_pions_data_cut_all.Fill(event.Q2)
        Epsilon_pions_data_cut_all.Fill(event.epsilon)
        t_pions_data_cut_all.Fill(-event.MandelT)
        ph_q_pions_data_cut_all.Fill(event.ph_q)
#    ibin += 1
'''
#Fill histograms for Prompt Data
for event in Cut_Pion_Events_Prompt_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        H_gtr_beta_pions_data_prompt_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_data_prompt_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_data_prompt_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_data_prompt_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_data_prompt_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_data_prompt_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_data_prompt_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_data_prompt_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_data_prompt_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_data_prompt_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_data_prompt_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_prompt_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_prompt_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_prompt_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_data_prompt_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_data_prompt_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_data_prompt_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_data_prompt_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_data_prompt_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_data_prompt_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_data_prompt_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_data_prompt_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_data_prompt_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_data_prompt_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_prompt_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_prompt_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_prompt_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_prompt_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_prompt_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_prompt_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_prompt_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_prompt_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_prompt_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_data_prompt_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_data_prompt_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_data_prompt_cut_all.Fill(event.pmiss)
        pmiss_x_pions_data_prompt_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_data_prompt_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_data_prompt_cut_all.Fill(event.pmiss_z)
        emiss_pions_data_prompt_cut_all.Fill(event.emiss)
        W_pions_data_prompt_cut_all.Fill(event.W)
        phi_deg = event.ph_q * (180 / math.pi)  # Convert φ to degrees
        if phi_deg < 0:
            phi_deg += 360  # Ensure φ is in [0, 360]
        ph_q_pions_data_prompt_cut_all.Fill(phi_deg)
        Q2_pions_data_prompt_cut_all.Fill(event.Q2)
        Epsilon_pions_data_prompt_cut_all.Fill(event.epsilon)
        t_pions_data_prompt_cut_all.Fill(-event.MandelT)
        pmiss_pions_data_prompt_cut_all_error.Fill(event.pmiss)
#    ibin += 1

#Fill histograms for Random Data
for event in Cut_Pion_Events_Random_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        H_gtr_beta_pions_data_random_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_data_random_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_data_random_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_data_random_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_data_random_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_data_random_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_data_random_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_data_random_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_data_random_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_data_random_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_data_random_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_random_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_random_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_random_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_data_random_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_data_random_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_data_random_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_data_random_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_data_random_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_data_random_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_data_random_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_data_random_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_data_random_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_data_random_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_random_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_random_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_random_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_random_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_random_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_random_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_random_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_random_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_random_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_data_random_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_data_random_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_data_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_data_random_cut_all.Fill(event.pmiss)
        pmiss_x_pions_data_random_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_data_random_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_data_random_cut_all.Fill(event.pmiss_z)
        emiss_pions_data_random_cut_all.Fill(event.emiss)
        W_pions_data_random_cut_all.Fill(event.W)
        phi_deg = event.ph_q * (180 / math.pi)  # Convert φ to degrees
        if phi_deg < 0:
            phi_deg += 360  # Ensure φ is in [0, 360]
        ph_q_pions_data_random_cut_all.Fill(phi_deg)
        Q2_pions_data_random_cut_all.Fill(event.Q2)
        Epsilon_pions_data_random_cut_all.Fill(event.epsilon)
        t_pions_data_random_cut_all.Fill(-event.MandelT)
        pmiss_pions_data_random_cut_all_error.Fill(event.pmiss)
#    ibin += 1

# Fill histograms from Prompt Dummy
#ibin = 1
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        H_gtr_beta_pions_dummy_prompt_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_dummy_prompt_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_dummy_prompt_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_dummy_prompt_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_dummy_prompt_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_dummy_prompt_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_dummy_prompt_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_dummy_prompt_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_dummy_prompt_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_dummy_prompt_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_dummy_prompt_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_dummy_prompt_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_dummy_prompt_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_dummy_prompt_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_dummy_prompt_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_dummy_prompt_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_dummy_prompt_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_dummy_prompt_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_dummy_prompt_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_dummy_prompt_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_dummy_prompt_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_dummy_prompt_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_dummy_prompt_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_dummy_prompt_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_dummy_prompt_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_dummy_prompt_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_dummy_prompt_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_dummy_prompt_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_dummy_prompt_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_dummy_prompt_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_dummy_prompt_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_dummy_prompt_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_dummy_prompt_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_dummy_prompt_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_dummy_prompt_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_dummy_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_dummy_prompt_cut_all.Fill(event.pmiss)
        pmiss_x_pions_dummy_prompt_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_dummy_prompt_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_dummy_prompt_cut_all.Fill(event.pmiss_z)
        emiss_pions_dummy_prompt_cut_all.Fill(event.emiss)
        W_pions_dummy_prompt_cut_all.Fill(event.W)
        phi_deg = event.ph_q * (180 / math.pi)  # Convert φ to degrees
        if phi_deg < 0:
            phi_deg += 360  # Ensure φ is in [0, 360]
        ph_q_pions_dummy_prompt_cut_all.Fill(phi_deg)
        Q2_pions_dummy_prompt_cut_all.Fill(event.Q2)
        Epsilon_pions_dummy_prompt_cut_all.Fill(event.epsilon)
        t_pions_dummy_prompt_cut_all.Fill(-event.MandelT)
        pmiss_pions_dummy_prompt_cut_all_error.Fill(event.pmiss)

#    ibin += 1

# Fill histograms from Random Dummy
#ibin = 1
for event in Cut_Pion_Events_Random_Dummy_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        H_gtr_beta_pions_dummy_random_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_pions_dummy_random_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_pions_dummy_random_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_pions_dummy_random_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_pions_dummy_random_cut_all.Fill(event.H_gtr_p)
        H_dc_x_fp_pions_dummy_random_cut_all.Fill(event.H_dc_x_fp)
        H_dc_y_fp_pions_dummy_random_cut_all.Fill(event.H_dc_y_fp)
        H_dc_xp_fp_pions_dummy_random_cut_all.Fill(event.H_dc_xp_fp)
        H_dc_yp_fp_pions_dummy_random_cut_all.Fill(event.H_dc_yp_fp)
        H_cal_etotnorm_pions_dummy_random_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_pions_dummy_random_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_dummy_random_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_dummy_random_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_dummy_random_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_pions_dummy_random_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_pions_dummy_random_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_pions_dummy_random_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_pions_dummy_random_cut_all.Fill(event.P_gtr_p)
        P_dc_x_fp_pions_dummy_random_cut_all.Fill(event.P_dc_x_fp)
        P_dc_y_fp_pions_dummy_random_cut_all.Fill(event.P_dc_y_fp)
        P_dc_xp_fp_pions_dummy_random_cut_all.Fill(event.P_dc_xp_fp)
        P_dc_yp_fp_pions_dummy_random_cut_all.Fill(event.P_dc_yp_fp)
        P_cal_etotnorm_pions_dummy_random_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_pions_dummy_random_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_dummy_random_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_dummy_random_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_dummy_random_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_dummy_random_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_dummy_random_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_dummy_random_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_dummy_random_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_dummy_random_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_dummy_random_cut_all.Fill(event.P_aero_yAtAero)
        MMpi_pions_dummy_random_cut_all.Fill(event.MMpi + MM_Offset)
        P_RFTime_Dist_pions_dummy_random_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        pmiss_pions_dummy_random_cut_all.Fill(event.pmiss)
        pmiss_x_pions_dummy_random_cut_all.Fill(event.pmiss_x)
        pmiss_y_pions_dummy_random_cut_all.Fill(event.pmiss_y)
        pmiss_z_pions_dummy_random_cut_all.Fill(event.pmiss_z)
        emiss_pions_dummy_random_cut_all.Fill(event.emiss)
        W_pions_dummy_random_cut_all.Fill(event.W)
        phi_deg = event.ph_q * (180 / math.pi)  # Convert φ to degrees
        if phi_deg < 0:
            phi_deg += 360  # Ensure φ is in [0, 360]
        ph_q_pions_dummy_random_cut_all.Fill(phi_deg)
        Q2_pions_dummy_random_cut_all.Fill(event.Q2)
        Epsilon_pions_dummy_random_cut_all.Fill(event.epsilon)
        t_pions_dummy_random_cut_all.Fill(-event.MandelT)
        pmiss_pions_dummy_random_cut_all_error.Fill(event.pmiss)
#    ibin += 1

# Fill histograms from SIMC ROOT File
for event in Uncut_Pion_Events_SIMC_tree:
    if HMS_Acceptance(event) and SHMS_Acceptance(event) and SHMS_Aero_Cut(event) and SIMC_MMpi_Cut(event) and Diamond_Cut(event):
            H_hsdelta_pions_simc_cut_all.Fill(event.hsdelta, event.Weight)
            H_hsxptar_pions_simc_cut_all.Fill(event.hsxptar, event.Weight)
            H_hsyptar_pions_simc_cut_all.Fill(event.hsyptar, event.Weight)
            H_hsytar_pions_simc_cut_all.Fill(event.hsytar, event.Weight)
            H_hsxfp_pions_simc_cut_all.Fill(event.hsxfp, event.Weight)
            H_hsyfp_pions_simc_cut_all.Fill(event.hsyfp, event.Weight)
            H_hsxpfp_pions_simc_cut_all.Fill(event.hsxpfp, event.Weight)
            H_hsypfp_pions_simc_cut_all.Fill(event.hsypfp, event.Weight)
            P_ssdelta_pions_simc_cut_all.Fill(event.ssdelta, event.Weight)
            P_ssxptar_pions_simc_cut_all.Fill(event.ssxptar, event.Weight)
            P_ssyptar_pions_simc_cut_all.Fill(event.ssyptar, event.Weight)
            P_ssytar_pions_simc_cut_all.Fill(event.ssytar, event.Weight)
            P_ssxfp_pions_simc_cut_all.Fill(event.ssxfp, event.Weight)
            P_ssyfp_pions_simc_cut_all.Fill(event.ssyfp, event.Weight)
            P_ssxpfp_pions_simc_cut_all.Fill(event.ssxpfp, event.Weight)
            P_ssypfp_pions_simc_cut_all.Fill(event.ssypfp, event.Weight)
            q_pions_simc_cut_all.Fill(event.q, event.Weight)
            nu_pions_simc_cut_all.Fill(event.nu, event.Weight)
            Q2_pions_simc_cut_all.Fill(event.Q2, event.Weight)
            epsilon_pions_simc_cut_all.Fill(event.epsilon, event.Weight)
            thetapq_pions_simc_cut_all.Fill(event.thetapq, event.Weight)
            phi_deg = event.phipq * (180 / math.pi)  # Convert φ to degrees
            if phi_deg < 0:
                phi_deg += 360  # Ensure φ is in [0, 360]
            phipq_pions_simc_cut_all.Fill(phi_deg, event.Weight)
            pmiss_pions_simc_cut_all.Fill(event.Pm, event.Weight)
            pmiss_x_pions_simc_cut_all.Fill(event.pmper, event.Weight)
            pmiss_y_pions_simc_cut_all.Fill(event.pmoop, event.Weight)
            pmiss_z_pions_simc_cut_all.Fill(event.pmpar, event.Weight)
            emiss_pions_simc_cut_all.Fill(event.Em, event.Weight)
            W_pions_simc_cut_all.Fill(event.W, event.Weight)
            t_pions_simc_cut_all.Fill(-event.t, event.Weight)
            MMpi_pions_simc_cut_all.Fill(event.missmass, event.Weight)
            pmiss_pions_simc_cut_all_error.Fill(event.Pm, event.Weight)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
H_gtr_beta_pions_data_random_cut_all.Scale(1.0/nWindows)
H_gtr_xp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_gtr_yp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_gtr_dp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_gtr_p_pions_data_random_cut_all.Scale(1.0/nWindows)
H_dc_x_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_dc_y_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_dc_xp_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_dc_yp_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cal_etotnorm_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cal_etottracknorm_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
H_RFTime_Dist_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_beta_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_xp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_yp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_dp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_p_pions_data_random_cut_all.Scale(1.0/nWindows)
P_dc_x_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_dc_y_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_dc_xp_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_dc_yp_fp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_cal_etotnorm_pions_data_random_cut_all.Scale(1.0/nWindows)
P_cal_etottracknorm_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_xAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_yAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_xAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_yAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_xAtAero_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_yAtAero_pions_data_random_cut_all.Scale(1.0/nWindows)
MMpi_pions_data_random_cut_all.Scale(1.0/nWindows)
P_RFTime_Dist_pions_data_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_data_random_cut_all.Scale(1.0/nWindows)
pmiss_pions_data_random_cut_all.Scale(1.0/nWindows)
pmiss_x_pions_data_random_cut_all.Scale(1.0/nWindows)
pmiss_y_pions_data_random_cut_all.Scale(1.0/nWindows)
pmiss_z_pions_data_random_cut_all.Scale(1.0/nWindows)
emiss_pions_data_random_cut_all.Scale(1.0/nWindows)
W_pions_data_random_cut_all.Scale(1.0/nWindows)
Q2_pions_data_random_cut_all.Scale(1.0/nWindows)
Epsilon_pions_data_random_cut_all.Scale(1.0/nWindows)
t_pions_data_random_cut_all.Scale(1.0/nWindows)
ph_q_pions_data_random_cut_all.Scale(1.0/nWindows)
pmiss_pions_data_random_cut_all_error.Scale(1.0/nWindows)

H_gtr_beta_pions_randsub_data_cut_all.Add(H_gtr_beta_pions_data_prompt_cut_all, H_gtr_beta_pions_data_random_cut_all, 1, -1)
H_gtr_xp_pions_randsub_data_cut_all.Add(H_gtr_xp_pions_data_prompt_cut_all, H_gtr_xp_pions_data_random_cut_all, 1, -1)
H_gtr_yp_pions_randsub_data_cut_all.Add(H_gtr_yp_pions_data_prompt_cut_all, H_gtr_yp_pions_data_random_cut_all, 1, -1)
H_gtr_dp_pions_randsub_data_cut_all.Add(H_gtr_dp_pions_data_prompt_cut_all, H_gtr_dp_pions_data_random_cut_all, 1, -1)
H_gtr_p_pions_randsub_data_cut_all.Add(H_gtr_p_pions_data_prompt_cut_all, H_gtr_p_pions_data_random_cut_all, 1, -1)
H_dc_x_fp_pions_randsub_data_cut_all.Add(H_dc_x_fp_pions_data_prompt_cut_all, H_dc_x_fp_pions_data_random_cut_all, 1, -1)
H_dc_y_fp_pions_randsub_data_cut_all.Add(H_dc_y_fp_pions_data_prompt_cut_all, H_dc_y_fp_pions_data_random_cut_all, 1, -1)
H_dc_xp_fp_pions_randsub_data_cut_all.Add(H_dc_xp_fp_pions_data_prompt_cut_all, H_dc_xp_fp_pions_data_random_cut_all, 1, -1)
H_dc_yp_fp_pions_randsub_data_cut_all.Add(H_dc_yp_fp_pions_data_prompt_cut_all, H_dc_yp_fp_pions_data_random_cut_all, 1, -1)
H_cal_etotnorm_pions_randsub_data_cut_all.Add(H_cal_etotnorm_pions_data_prompt_cut_all, H_cal_etotnorm_pions_data_random_cut_all, 1, -1)
H_cal_etottracknorm_pions_randsub_data_cut_all.Add(H_cal_etottracknorm_pions_data_prompt_cut_all, H_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
H_cer_npeSum_pions_randsub_data_cut_all.Add(H_cer_npeSum_pions_data_prompt_cut_all, H_cer_npeSum_pions_data_random_cut_all, 1, -1)
H_RFTime_Dist_pions_randsub_data_cut_all.Add(H_RFTime_Dist_pions_data_prompt_cut_all, H_RFTime_Dist_pions_data_random_cut_all, 1, -1)
P_gtr_beta_pions_randsub_data_cut_all.Add(P_gtr_beta_pions_data_prompt_cut_all, P_gtr_beta_pions_data_random_cut_all, 1, -1)
P_gtr_xp_pions_randsub_data_cut_all.Add(P_gtr_xp_pions_data_prompt_cut_all, P_gtr_xp_pions_data_random_cut_all, 1, -1)
P_gtr_yp_pions_randsub_data_cut_all.Add(P_gtr_yp_pions_data_prompt_cut_all, P_gtr_yp_pions_data_random_cut_all, 1, -1)
P_gtr_dp_pions_randsub_data_cut_all.Add(P_gtr_dp_pions_data_prompt_cut_all, P_gtr_dp_pions_data_random_cut_all, 1, -1)
P_gtr_p_pions_randsub_data_cut_all.Add(P_gtr_p_pions_data_prompt_cut_all, P_gtr_p_pions_data_random_cut_all, 1, -1)
P_dc_x_fp_pions_randsub_data_cut_all.Add(P_dc_x_fp_pions_data_prompt_cut_all, P_dc_x_fp_pions_data_random_cut_all, 1, -1)
P_dc_y_fp_pions_randsub_data_cut_all.Add(P_dc_y_fp_pions_data_prompt_cut_all, P_dc_y_fp_pions_data_random_cut_all, 1, -1)
P_dc_xp_fp_pions_randsub_data_cut_all.Add(P_dc_xp_fp_pions_data_prompt_cut_all, P_dc_xp_fp_pions_data_random_cut_all, 1, -1)
P_dc_yp_fp_pions_randsub_data_cut_all.Add(P_dc_yp_fp_pions_data_prompt_cut_all, P_dc_yp_fp_pions_data_random_cut_all, 1, -1)
P_cal_etotnorm_pions_randsub_data_cut_all.Add(P_cal_etotnorm_pions_data_prompt_cut_all, P_cal_etotnorm_pions_data_random_cut_all, 1, -1)
P_cal_etottracknorm_pions_randsub_data_cut_all.Add(P_cal_etottracknorm_pions_data_prompt_cut_all, P_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
P_hgcer_npeSum_pions_randsub_data_cut_all.Add(P_hgcer_npeSum_pions_data_prompt_cut_all, P_hgcer_npeSum_pions_data_random_cut_all, 1, -1)
P_hgcer_xAtCer_pions_randsub_data_cut_all.Add(P_hgcer_xAtCer_pions_data_prompt_cut_all, P_hgcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_hgcer_yAtCer_pions_randsub_data_cut_all.Add(P_hgcer_yAtCer_pions_data_prompt_cut_all, P_hgcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_npeSum_pions_randsub_data_cut_all.Add(P_ngcer_npeSum_pions_data_prompt_cut_all, P_ngcer_npeSum_pions_data_random_cut_all, 1, -1)
P_ngcer_xAtCer_pions_randsub_data_cut_all.Add(P_ngcer_xAtCer_pions_data_prompt_cut_all, P_ngcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_yAtCer_pions_randsub_data_cut_all.Add(P_ngcer_yAtCer_pions_data_prompt_cut_all, P_ngcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_aero_npeSum_pions_randsub_data_cut_all.Add(P_aero_npeSum_pions_data_prompt_cut_all, P_aero_npeSum_pions_data_random_cut_all, 1, -1)
P_aero_xAtAero_pions_randsub_data_cut_all.Add(P_aero_xAtAero_pions_data_prompt_cut_all, P_aero_xAtAero_pions_data_random_cut_all, 1, -1)
P_aero_yAtAero_pions_randsub_data_cut_all.Add(P_aero_yAtAero_pions_data_prompt_cut_all, P_aero_yAtAero_pions_data_random_cut_all, 1, -1)
MMpi_pions_randsub_data_cut_all.Add(MMpi_pions_data_prompt_cut_all, MMpi_pions_data_random_cut_all, 1, -1)
P_RFTime_Dist_pions_randsub_data_cut_all.Add(P_RFTime_Dist_pions_data_prompt_cut_all, P_RFTime_Dist_pions_data_random_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_randsub_data_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_cut_all, 1, -1)
pmiss_pions_randsub_data_cut_all.Add(pmiss_pions_data_prompt_cut_all, pmiss_pions_data_random_cut_all, 1, -1)
pmiss_x_pions_randsub_data_cut_all.Add(pmiss_x_pions_data_prompt_cut_all, pmiss_x_pions_data_random_cut_all, 1, -1)
pmiss_y_pions_randsub_data_cut_all.Add(pmiss_y_pions_data_prompt_cut_all, pmiss_y_pions_data_random_cut_all, 1, -1)
pmiss_z_pions_randsub_data_cut_all.Add(pmiss_z_pions_data_prompt_cut_all, pmiss_z_pions_data_random_cut_all, 1, -1)
emiss_pions_randsub_data_cut_all.Add(emiss_pions_data_prompt_cut_all, emiss_pions_data_random_cut_all, 1, -1)
W_pions_randsub_data_cut_all.Add(W_pions_data_prompt_cut_all, W_pions_data_random_cut_all, 1, -1)
Q2_pions_randsub_data_cut_all.Add(Q2_pions_data_prompt_cut_all, Q2_pions_data_random_cut_all, 1, -1)
Epsilon_pions_randsub_data_cut_all.Add(Epsilon_pions_data_prompt_cut_all, Epsilon_pions_data_random_cut_all, 1, -1)
t_pions_randsub_data_cut_all.Add(t_pions_data_prompt_cut_all, t_pions_data_random_cut_all, 1, -1)
ph_q_pions_randsub_data_cut_all.Add(ph_q_pions_data_prompt_cut_all, ph_q_pions_data_random_cut_all, 1, -1)
pmiss_pions_randsub_data_cut_all_error.Add(pmiss_pions_data_prompt_cut_all_error, pmiss_pions_data_random_cut_all_error, 1, -1)

H_gtr_beta_pions_normrandsub_data_cut_all.Add(H_gtr_beta_pions_data_prompt_cut_all, H_gtr_beta_pions_data_random_cut_all, 1, -1)
H_gtr_xp_pions_normrandsub_data_cut_all.Add(H_gtr_xp_pions_data_prompt_cut_all, H_gtr_xp_pions_data_random_cut_all, 1, -1)
H_gtr_yp_pions_normrandsub_data_cut_all.Add(H_gtr_yp_pions_data_prompt_cut_all, H_gtr_yp_pions_data_random_cut_all, 1, -1)
H_gtr_dp_pions_normrandsub_data_cut_all.Add(H_gtr_dp_pions_data_prompt_cut_all, H_gtr_dp_pions_data_random_cut_all, 1, -1)
H_gtr_p_pions_normrandsub_data_cut_all.Add(H_gtr_p_pions_data_prompt_cut_all, H_gtr_p_pions_data_random_cut_all, 1, -1)
H_dc_x_fp_pions_normrandsub_data_cut_all.Add(H_dc_x_fp_pions_data_prompt_cut_all, H_dc_x_fp_pions_data_random_cut_all, 1, -1)
H_dc_y_fp_pions_normrandsub_data_cut_all.Add(H_dc_y_fp_pions_data_prompt_cut_all, H_dc_y_fp_pions_data_random_cut_all, 1, -1)
H_dc_xp_fp_pions_normrandsub_data_cut_all.Add(H_dc_xp_fp_pions_data_prompt_cut_all, H_dc_xp_fp_pions_data_random_cut_all, 1, -1)
H_dc_yp_fp_pions_normrandsub_data_cut_all.Add(H_dc_yp_fp_pions_data_prompt_cut_all, H_dc_yp_fp_pions_data_random_cut_all, 1, -1)
H_cal_etotnorm_pions_normrandsub_data_cut_all.Add(H_cal_etotnorm_pions_data_prompt_cut_all, H_cal_etotnorm_pions_data_random_cut_all, 1, -1)
H_cal_etottracknorm_pions_normrandsub_data_cut_all.Add(H_cal_etottracknorm_pions_data_prompt_cut_all, H_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
H_cer_npeSum_pions_normrandsub_data_cut_all.Add(H_cer_npeSum_pions_data_prompt_cut_all, H_cer_npeSum_pions_data_random_cut_all, 1, -1)
H_RFTime_Dist_pions_normrandsub_data_cut_all.Add(H_RFTime_Dist_pions_data_prompt_cut_all, H_RFTime_Dist_pions_data_random_cut_all, 1, -1)
P_gtr_beta_pions_normrandsub_data_cut_all.Add(P_gtr_beta_pions_data_prompt_cut_all, P_gtr_beta_pions_data_random_cut_all, 1, -1)
P_gtr_xp_pions_normrandsub_data_cut_all.Add(P_gtr_xp_pions_data_prompt_cut_all, P_gtr_xp_pions_data_random_cut_all, 1, -1)
P_gtr_yp_pions_normrandsub_data_cut_all.Add(P_gtr_yp_pions_data_prompt_cut_all, P_gtr_yp_pions_data_random_cut_all, 1, -1)
P_gtr_dp_pions_normrandsub_data_cut_all.Add(P_gtr_dp_pions_data_prompt_cut_all, P_gtr_dp_pions_data_random_cut_all, 1, -1)
P_gtr_p_pions_normrandsub_data_cut_all.Add(P_gtr_p_pions_data_prompt_cut_all, P_gtr_p_pions_data_random_cut_all, 1, -1)
P_dc_x_fp_pions_normrandsub_data_cut_all.Add(P_dc_x_fp_pions_data_prompt_cut_all, P_dc_x_fp_pions_data_random_cut_all, 1, -1)
P_dc_y_fp_pions_normrandsub_data_cut_all.Add(P_dc_y_fp_pions_data_prompt_cut_all, P_dc_y_fp_pions_data_random_cut_all, 1, -1)
P_dc_xp_fp_pions_normrandsub_data_cut_all.Add(P_dc_xp_fp_pions_data_prompt_cut_all, P_dc_xp_fp_pions_data_random_cut_all, 1, -1)
P_dc_yp_fp_pions_normrandsub_data_cut_all.Add(P_dc_yp_fp_pions_data_prompt_cut_all, P_dc_yp_fp_pions_data_random_cut_all, 1, -1)
P_cal_etotnorm_pions_normrandsub_data_cut_all.Add(P_cal_etotnorm_pions_data_prompt_cut_all, P_cal_etotnorm_pions_data_random_cut_all, 1, -1)
P_cal_etottracknorm_pions_normrandsub_data_cut_all.Add(P_cal_etottracknorm_pions_data_prompt_cut_all, P_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
P_hgcer_npeSum_pions_normrandsub_data_cut_all.Add(P_hgcer_npeSum_pions_data_prompt_cut_all, P_hgcer_npeSum_pions_data_random_cut_all, 1, -1)
P_hgcer_xAtCer_pions_normrandsub_data_cut_all.Add(P_hgcer_xAtCer_pions_data_prompt_cut_all, P_hgcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_hgcer_yAtCer_pions_normrandsub_data_cut_all.Add(P_hgcer_yAtCer_pions_data_prompt_cut_all, P_hgcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_npeSum_pions_normrandsub_data_cut_all.Add(P_ngcer_npeSum_pions_data_prompt_cut_all, P_ngcer_npeSum_pions_data_random_cut_all, 1, -1)
P_ngcer_xAtCer_pions_normrandsub_data_cut_all.Add(P_ngcer_xAtCer_pions_data_prompt_cut_all, P_ngcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_yAtCer_pions_normrandsub_data_cut_all.Add(P_ngcer_yAtCer_pions_data_prompt_cut_all, P_ngcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_aero_npeSum_pions_normrandsub_data_cut_all.Add(P_aero_npeSum_pions_data_prompt_cut_all, P_aero_npeSum_pions_data_random_cut_all, 1, -1)
P_aero_xAtAero_pions_normrandsub_data_cut_all.Add(P_aero_xAtAero_pions_data_prompt_cut_all, P_aero_xAtAero_pions_data_random_cut_all, 1, -1)
P_aero_yAtAero_pions_normrandsub_data_cut_all.Add(P_aero_yAtAero_pions_data_prompt_cut_all, P_aero_yAtAero_pions_data_random_cut_all, 1, -1)
MMpi_pions_normrandsub_data_cut_all.Add(MMpi_pions_data_prompt_cut_all, MMpi_pions_data_random_cut_all, 1, -1)
P_RFTime_Dist_pions_normrandsub_data_cut_all.Add(P_RFTime_Dist_pions_data_prompt_cut_all, P_RFTime_Dist_pions_data_random_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_normrandsub_data_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_cut_all, 1, -1)
pmiss_pions_normrandsub_data_cut_all.Add(pmiss_pions_data_prompt_cut_all, pmiss_pions_data_random_cut_all, 1, -1)
pmiss_x_pions_normrandsub_data_cut_all.Add(pmiss_x_pions_data_prompt_cut_all, pmiss_x_pions_data_random_cut_all, 1, -1)
pmiss_y_pions_normrandsub_data_cut_all.Add(pmiss_y_pions_data_prompt_cut_all, pmiss_y_pions_data_random_cut_all, 1, -1)
pmiss_z_pions_normrandsub_data_cut_all.Add(pmiss_z_pions_data_prompt_cut_all, pmiss_z_pions_data_random_cut_all, 1, -1)
emiss_pions_normrandsub_data_cut_all.Add(emiss_pions_data_prompt_cut_all, emiss_pions_data_random_cut_all, 1, -1)
W_pions_normrandsub_data_cut_all.Add(W_pions_data_prompt_cut_all, W_pions_data_random_cut_all, 1, -1)
Q2_pions_normrandsub_data_cut_all.Add(Q2_pions_data_prompt_cut_all, Q2_pions_data_random_cut_all, 1, -1)
Epsilon_pions_normrandsub_data_cut_all.Add(Epsilon_pions_data_prompt_cut_all, Epsilon_pions_data_random_cut_all, 1, -1)
t_pions_normrandsub_data_cut_all.Add(t_pions_data_prompt_cut_all, t_pions_data_random_cut_all, 1, -1)
ph_q_pions_normrandsub_data_cut_all.Add(ph_q_pions_data_prompt_cut_all, ph_q_pions_data_random_cut_all, 1, -1)

H_gtr_beta_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_gtr_xp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_gtr_yp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_gtr_dp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_gtr_p_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_dc_x_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_dc_y_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_dc_xp_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_dc_yp_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_cal_etotnorm_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_cal_etottracknorm_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_cer_npeSum_pions_dummy_random_cut_all.Scale(1.0/nWindows)
H_RFTime_Dist_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_gtr_beta_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_gtr_xp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_gtr_yp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_gtr_dp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_gtr_p_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_dc_x_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_dc_y_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_dc_xp_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_dc_yp_fp_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_cal_etotnorm_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_cal_etottracknorm_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_hgcer_npeSum_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_hgcer_xAtCer_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_hgcer_yAtCer_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_ngcer_npeSum_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_ngcer_xAtCer_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_ngcer_yAtCer_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_aero_xAtAero_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_aero_yAtAero_pions_dummy_random_cut_all.Scale(1.0/nWindows)
MMpi_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_RFTime_Dist_pions_dummy_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all.Scale(1.0/nWindows)
pmiss_pions_dummy_random_cut_all.Scale(1.0/nWindows)
pmiss_x_pions_dummy_random_cut_all.Scale(1.0/nWindows)
pmiss_y_pions_dummy_random_cut_all.Scale(1.0/nWindows)
pmiss_z_pions_dummy_random_cut_all.Scale(1.0/nWindows)
emiss_pions_dummy_random_cut_all.Scale(1.0/nWindows)
W_pions_dummy_random_cut_all.Scale(1.0/nWindows)
Q2_pions_dummy_random_cut_all.Scale(1.0/nWindows)
Epsilon_pions_dummy_random_cut_all.Scale(1.0/nWindows)
t_pions_dummy_random_cut_all.Scale(1.0/nWindows)
ph_q_pions_dummy_random_cut_all.Scale(1.0/nWindows)
pmiss_pions_dummy_random_cut_all_error.Scale(1.0/nWindows)

H_gtr_beta_pions_randsub_dummy_cut_all.Add(H_gtr_beta_pions_dummy_prompt_cut_all, H_gtr_beta_pions_dummy_random_cut_all, 1, -1)
H_gtr_xp_pions_randsub_dummy_cut_all.Add(H_gtr_xp_pions_dummy_prompt_cut_all, H_gtr_xp_pions_dummy_random_cut_all, 1, -1)
H_gtr_yp_pions_randsub_dummy_cut_all.Add(H_gtr_yp_pions_dummy_prompt_cut_all, H_gtr_yp_pions_dummy_random_cut_all, 1, -1)
H_gtr_dp_pions_randsub_dummy_cut_all.Add(H_gtr_dp_pions_dummy_prompt_cut_all, H_gtr_dp_pions_dummy_random_cut_all, 1, -1)
H_gtr_p_pions_randsub_dummy_cut_all.Add(H_gtr_p_pions_dummy_prompt_cut_all, H_gtr_p_pions_dummy_random_cut_all, 1, -1)
H_dc_x_fp_pions_randsub_dummy_cut_all.Add(H_dc_x_fp_pions_dummy_prompt_cut_all, H_dc_x_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_y_fp_pions_randsub_dummy_cut_all.Add(H_dc_y_fp_pions_dummy_prompt_cut_all, H_dc_y_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_xp_fp_pions_randsub_dummy_cut_all.Add(H_dc_xp_fp_pions_dummy_prompt_cut_all, H_dc_xp_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_yp_fp_pions_randsub_dummy_cut_all.Add(H_dc_yp_fp_pions_dummy_prompt_cut_all, H_dc_yp_fp_pions_dummy_random_cut_all, 1, -1)
H_cal_etotnorm_pions_randsub_dummy_cut_all.Add(H_cal_etotnorm_pions_dummy_prompt_cut_all, H_cal_etotnorm_pions_dummy_random_cut_all, 1, -1)
H_cal_etottracknorm_pions_randsub_dummy_cut_all.Add(H_cal_etottracknorm_pions_dummy_prompt_cut_all, H_cal_etottracknorm_pions_dummy_random_cut_all, 1, -1)
H_cer_npeSum_pions_randsub_dummy_cut_all.Add(H_cer_npeSum_pions_dummy_prompt_cut_all, H_cer_npeSum_pions_dummy_random_cut_all, 1, -1)
H_RFTime_Dist_pions_randsub_dummy_cut_all.Add(H_RFTime_Dist_pions_dummy_prompt_cut_all, H_RFTime_Dist_pions_dummy_random_cut_all, 1, -1)
P_gtr_beta_pions_randsub_dummy_cut_all.Add(P_gtr_beta_pions_dummy_prompt_cut_all, P_gtr_beta_pions_dummy_random_cut_all, 1, -1)
P_gtr_xp_pions_randsub_dummy_cut_all.Add(P_gtr_xp_pions_dummy_prompt_cut_all, P_gtr_xp_pions_dummy_random_cut_all, 1, -1)
P_gtr_yp_pions_randsub_dummy_cut_all.Add(P_gtr_yp_pions_dummy_prompt_cut_all, P_gtr_yp_pions_dummy_random_cut_all, 1, -1)
P_gtr_dp_pions_randsub_dummy_cut_all.Add(P_gtr_dp_pions_dummy_prompt_cut_all, P_gtr_dp_pions_dummy_random_cut_all, 1, -1)
P_gtr_p_pions_randsub_dummy_cut_all.Add(P_gtr_p_pions_dummy_prompt_cut_all, P_gtr_p_pions_dummy_random_cut_all, 1, -1)
P_dc_x_fp_pions_randsub_dummy_cut_all.Add(P_dc_x_fp_pions_dummy_prompt_cut_all, P_dc_x_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_y_fp_pions_randsub_dummy_cut_all.Add(P_dc_y_fp_pions_dummy_prompt_cut_all, P_dc_y_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_xp_fp_pions_randsub_dummy_cut_all.Add(P_dc_xp_fp_pions_dummy_prompt_cut_all, P_dc_xp_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_yp_fp_pions_randsub_dummy_cut_all.Add(P_dc_yp_fp_pions_dummy_prompt_cut_all, P_dc_yp_fp_pions_dummy_random_cut_all, 1, -1)
P_cal_etotnorm_pions_randsub_dummy_cut_all.Add(P_cal_etotnorm_pions_dummy_prompt_cut_all, P_cal_etotnorm_pions_dummy_random_cut_all, 1, -1)
P_cal_etottracknorm_pions_randsub_dummy_cut_all.Add(P_cal_etottracknorm_pions_dummy_prompt_cut_all, P_cal_etottracknorm_pions_dummy_random_cut_all, 1, -1)
P_hgcer_npeSum_pions_randsub_dummy_cut_all.Add(P_hgcer_npeSum_pions_dummy_prompt_cut_all, P_hgcer_npeSum_pions_dummy_random_cut_all, 1, -1)
P_hgcer_xAtCer_pions_randsub_dummy_cut_all.Add(P_hgcer_xAtCer_pions_dummy_prompt_cut_all, P_hgcer_xAtCer_pions_dummy_random_cut_all, 1, -1)
P_hgcer_yAtCer_pions_randsub_dummy_cut_all.Add(P_hgcer_yAtCer_pions_dummy_prompt_cut_all, P_hgcer_yAtCer_pions_dummy_random_cut_all, 1, -1)
P_ngcer_npeSum_pions_randsub_dummy_cut_all.Add(P_ngcer_npeSum_pions_dummy_prompt_cut_all, P_ngcer_npeSum_pions_dummy_random_cut_all, 1, -1)
P_ngcer_xAtCer_pions_randsub_dummy_cut_all.Add(P_ngcer_xAtCer_pions_dummy_prompt_cut_all, P_ngcer_xAtCer_pions_dummy_random_cut_all, 1, -1)
P_ngcer_yAtCer_pions_randsub_dummy_cut_all.Add(P_ngcer_yAtCer_pions_dummy_prompt_cut_all, P_ngcer_yAtCer_pions_dummy_random_cut_all, 1, -1)
P_aero_npeSum_pions_randsub_dummy_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_cut_all, P_aero_npeSum_pions_dummy_random_cut_all, 1, -1)
P_aero_xAtAero_pions_randsub_dummy_cut_all.Add(P_aero_xAtAero_pions_dummy_prompt_cut_all, P_aero_xAtAero_pions_dummy_random_cut_all, 1, -1)
P_aero_yAtAero_pions_randsub_dummy_cut_all.Add(P_aero_yAtAero_pions_dummy_prompt_cut_all, P_aero_yAtAero_pions_dummy_random_cut_all, 1, -1)
MMpi_pions_randsub_dummy_cut_all.Add(MMpi_pions_dummy_prompt_cut_all, MMpi_pions_dummy_random_cut_all, 1, -1)
P_RFTime_Dist_pions_randsub_dummy_cut_all.Add(P_RFTime_Dist_pions_dummy_prompt_cut_all, P_RFTime_Dist_pions_dummy_random_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_randsub_dummy_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all, 1, -1)
pmiss_pions_randsub_dummy_cut_all.Add(pmiss_pions_dummy_prompt_cut_all, pmiss_pions_dummy_random_cut_all, 1, -1)
pmiss_x_pions_randsub_dummy_cut_all.Add(pmiss_x_pions_dummy_prompt_cut_all, pmiss_x_pions_dummy_random_cut_all, 1, -1)
pmiss_y_pions_randsub_dummy_cut_all.Add(pmiss_y_pions_dummy_prompt_cut_all, pmiss_y_pions_dummy_random_cut_all, 1, -1)
pmiss_z_pions_randsub_dummy_cut_all.Add(pmiss_z_pions_dummy_prompt_cut_all, pmiss_z_pions_dummy_random_cut_all, 1, -1)
emiss_pions_randsub_dummy_cut_all.Add(emiss_pions_dummy_prompt_cut_all, emiss_pions_dummy_random_cut_all, 1, -1)
W_pions_randsub_dummy_cut_all.Add(W_pions_dummy_prompt_cut_all, W_pions_dummy_random_cut_all, 1, -1)
Q2_pions_randsub_dummy_cut_all.Add(Q2_pions_dummy_prompt_cut_all, Q2_pions_dummy_random_cut_all, 1, -1)
Epsilon_pions_randsub_dummy_cut_all.Add(Epsilon_pions_dummy_prompt_cut_all, Epsilon_pions_dummy_random_cut_all, 1, -1)
t_pions_randsub_dummy_cut_all.Add(t_pions_dummy_prompt_cut_all, t_pions_dummy_random_cut_all, 1, -1)
ph_q_pions_randsub_dummy_cut_all.Add(ph_q_pions_dummy_prompt_cut_all, ph_q_pions_dummy_random_cut_all, 1, -1)
pmiss_pions_randsub_dummy_cut_all_error.Add(pmiss_pions_dummy_prompt_cut_all_error, pmiss_pions_dummy_random_cut_all_error, 1, -1)

H_gtr_beta_pions_normrandsub_dummy_cut_all.Add(H_gtr_beta_pions_dummy_prompt_cut_all, H_gtr_beta_pions_dummy_random_cut_all, 1, -1)
H_gtr_xp_pions_normrandsub_dummy_cut_all.Add(H_gtr_xp_pions_dummy_prompt_cut_all, H_gtr_xp_pions_dummy_random_cut_all, 1, -1)
H_gtr_yp_pions_normrandsub_dummy_cut_all.Add(H_gtr_yp_pions_dummy_prompt_cut_all, H_gtr_yp_pions_dummy_random_cut_all, 1, -1)
H_gtr_dp_pions_normrandsub_dummy_cut_all.Add(H_gtr_dp_pions_dummy_prompt_cut_all, H_gtr_dp_pions_dummy_random_cut_all, 1, -1)
H_gtr_p_pions_normrandsub_dummy_cut_all.Add(H_gtr_p_pions_dummy_prompt_cut_all, H_gtr_p_pions_dummy_random_cut_all, 1, -1)
H_dc_x_fp_pions_normrandsub_dummy_cut_all.Add(H_dc_x_fp_pions_dummy_prompt_cut_all, H_dc_x_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_y_fp_pions_normrandsub_dummy_cut_all.Add(H_dc_y_fp_pions_dummy_prompt_cut_all, H_dc_y_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_xp_fp_pions_normrandsub_dummy_cut_all.Add(H_dc_xp_fp_pions_dummy_prompt_cut_all, H_dc_xp_fp_pions_dummy_random_cut_all, 1, -1)
H_dc_yp_fp_pions_normrandsub_dummy_cut_all.Add(H_dc_yp_fp_pions_dummy_prompt_cut_all, H_dc_yp_fp_pions_dummy_random_cut_all, 1, -1)
H_cal_etotnorm_pions_normrandsub_dummy_cut_all.Add(H_cal_etotnorm_pions_dummy_prompt_cut_all, H_cal_etotnorm_pions_dummy_random_cut_all, 1, -1)
H_cal_etottracknorm_pions_normrandsub_dummy_cut_all.Add(H_cal_etottracknorm_pions_dummy_prompt_cut_all, H_cal_etottracknorm_pions_dummy_random_cut_all, 1, -1)
H_cer_npeSum_pions_normrandsub_dummy_cut_all.Add(H_cer_npeSum_pions_dummy_prompt_cut_all, H_cer_npeSum_pions_dummy_random_cut_all, 1, -1)
H_RFTime_Dist_pions_normrandsub_dummy_cut_all.Add(H_RFTime_Dist_pions_dummy_prompt_cut_all, H_RFTime_Dist_pions_dummy_random_cut_all, 1, -1)
P_gtr_beta_pions_normrandsub_dummy_cut_all.Add(P_gtr_beta_pions_dummy_prompt_cut_all, P_gtr_beta_pions_dummy_random_cut_all, 1, -1)
P_gtr_xp_pions_normrandsub_dummy_cut_all.Add(P_gtr_xp_pions_dummy_prompt_cut_all, P_gtr_xp_pions_dummy_random_cut_all, 1, -1)
P_gtr_yp_pions_normrandsub_dummy_cut_all.Add(P_gtr_yp_pions_dummy_prompt_cut_all, P_gtr_yp_pions_dummy_random_cut_all, 1, -1)
P_gtr_dp_pions_normrandsub_dummy_cut_all.Add(P_gtr_dp_pions_dummy_prompt_cut_all, P_gtr_dp_pions_dummy_random_cut_all, 1, -1)
P_gtr_p_pions_normrandsub_dummy_cut_all.Add(P_gtr_p_pions_dummy_prompt_cut_all, P_gtr_p_pions_dummy_random_cut_all, 1, -1)
P_dc_x_fp_pions_normrandsub_dummy_cut_all.Add(P_dc_x_fp_pions_dummy_prompt_cut_all, P_dc_x_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_y_fp_pions_normrandsub_dummy_cut_all.Add(P_dc_y_fp_pions_dummy_prompt_cut_all, P_dc_y_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_xp_fp_pions_normrandsub_dummy_cut_all.Add(P_dc_xp_fp_pions_dummy_prompt_cut_all, P_dc_xp_fp_pions_dummy_random_cut_all, 1, -1)
P_dc_yp_fp_pions_normrandsub_dummy_cut_all.Add(P_dc_yp_fp_pions_dummy_prompt_cut_all, P_dc_yp_fp_pions_dummy_random_cut_all, 1, -1)
P_cal_etotnorm_pions_normrandsub_dummy_cut_all.Add(P_cal_etotnorm_pions_dummy_prompt_cut_all, P_cal_etotnorm_pions_dummy_random_cut_all, 1, -1)
P_cal_etottracknorm_pions_normrandsub_dummy_cut_all.Add(P_cal_etottracknorm_pions_dummy_prompt_cut_all, P_cal_etottracknorm_pions_dummy_random_cut_all, 1, -1)
P_hgcer_npeSum_pions_normrandsub_dummy_cut_all.Add(P_hgcer_npeSum_pions_dummy_prompt_cut_all, P_hgcer_npeSum_pions_dummy_random_cut_all, 1, -1)
P_hgcer_xAtCer_pions_normrandsub_dummy_cut_all.Add(P_hgcer_xAtCer_pions_dummy_prompt_cut_all, P_hgcer_xAtCer_pions_dummy_random_cut_all, 1, -1)
P_hgcer_yAtCer_pions_normrandsub_dummy_cut_all.Add(P_hgcer_yAtCer_pions_dummy_prompt_cut_all, P_hgcer_yAtCer_pions_dummy_random_cut_all, 1, -1)
P_ngcer_npeSum_pions_normrandsub_dummy_cut_all.Add(P_ngcer_npeSum_pions_dummy_prompt_cut_all, P_ngcer_npeSum_pions_dummy_random_cut_all, 1, -1)
P_ngcer_xAtCer_pions_normrandsub_dummy_cut_all.Add(P_ngcer_xAtCer_pions_dummy_prompt_cut_all, P_ngcer_xAtCer_pions_dummy_random_cut_all, 1, -1)
P_ngcer_yAtCer_pions_normrandsub_dummy_cut_all.Add(P_ngcer_yAtCer_pions_dummy_prompt_cut_all, P_ngcer_yAtCer_pions_dummy_random_cut_all, 1, -1)
P_aero_npeSum_pions_normrandsub_dummy_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_cut_all, P_aero_npeSum_pions_dummy_random_cut_all, 1, -1)
P_aero_xAtAero_pions_normrandsub_dummy_cut_all.Add(P_aero_xAtAero_pions_dummy_prompt_cut_all, P_aero_xAtAero_pions_dummy_random_cut_all, 1, -1)
P_aero_yAtAero_pions_normrandsub_dummy_cut_all.Add(P_aero_yAtAero_pions_dummy_prompt_cut_all, P_aero_yAtAero_pions_dummy_random_cut_all, 1, -1)
MMpi_pions_normrandsub_dummy_cut_all.Add(MMpi_pions_dummy_prompt_cut_all, MMpi_pions_dummy_random_cut_all, 1, -1)
P_RFTime_Dist_pions_normrandsub_dummy_cut_all.Add(P_RFTime_Dist_pions_dummy_prompt_cut_all, P_RFTime_Dist_pions_dummy_random_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_normrandsub_dummy_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_cut_all, 1, -1)
pmiss_pions_normrandsub_dummy_cut_all.Add(pmiss_pions_dummy_prompt_cut_all, pmiss_pions_dummy_random_cut_all, 1, -1)
pmiss_x_pions_normrandsub_dummy_cut_all.Add(pmiss_x_pions_dummy_prompt_cut_all, pmiss_x_pions_dummy_random_cut_all, 1, -1)
pmiss_y_pions_normrandsub_dummy_cut_all.Add(pmiss_y_pions_dummy_prompt_cut_all, pmiss_y_pions_dummy_random_cut_all, 1, -1)
pmiss_z_pions_normrandsub_dummy_cut_all.Add(pmiss_z_pions_dummy_prompt_cut_all, pmiss_z_pions_dummy_random_cut_all, 1, -1)
emiss_pions_normrandsub_dummy_cut_all.Add(emiss_pions_dummy_prompt_cut_all, emiss_pions_dummy_random_cut_all, 1, -1)
W_pions_normrandsub_dummy_cut_all.Add(W_pions_dummy_prompt_cut_all, W_pions_dummy_random_cut_all, 1, -1)
Q2_pions_normrandsub_dummy_cut_all.Add(Q2_pions_dummy_prompt_cut_all, Q2_pions_dummy_random_cut_all, 1, -1)
Epsilon_pions_normrandsub_dummy_cut_all.Add(Epsilon_pions_dummy_prompt_cut_all, Epsilon_pions_dummy_random_cut_all, 1, -1)
t_pions_normrandsub_dummy_cut_all.Add(t_pions_dummy_prompt_cut_all, t_pions_dummy_random_cut_all, 1, -1)
ph_q_pions_normrandsub_dummy_cut_all.Add(ph_q_pions_dummy_prompt_cut_all, ph_q_pions_dummy_random_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Data Normalization
H_gtr_beta_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_gtr_xp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_gtr_yp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_gtr_dp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_gtr_p_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_dc_x_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_dc_y_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_dc_xp_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_dc_yp_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_cal_etotnorm_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_cal_etottracknorm_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_cer_npeSum_pions_normrandsub_data_cut_all.Scale(normfac_data)
H_RFTime_Dist_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_gtr_beta_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_gtr_xp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_gtr_yp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_gtr_dp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_gtr_p_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_dc_x_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_dc_y_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_dc_xp_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_dc_yp_fp_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_cal_etotnorm_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_cal_etottracknorm_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_hgcer_npeSum_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_hgcer_xAtCer_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_hgcer_yAtCer_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_ngcer_npeSum_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_ngcer_xAtCer_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_ngcer_yAtCer_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_aero_npeSum_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_aero_xAtAero_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_aero_yAtAero_pions_normrandsub_data_cut_all.Scale(normfac_data)
MMpi_pions_normrandsub_data_cut_all.Scale(normfac_data)
P_RFTime_Dist_pions_normrandsub_data_cut_all.Scale(normfac_data)
CTime_ePiCoinTime_ROC2_pions_normrandsub_data_cut_all.Scale(normfac_data)
pmiss_pions_normrandsub_data_cut_all.Scale(normfac_data)
pmiss_x_pions_normrandsub_data_cut_all.Scale(normfac_data)
pmiss_y_pions_normrandsub_data_cut_all.Scale(normfac_data)
pmiss_z_pions_normrandsub_data_cut_all.Scale(normfac_data)
emiss_pions_normrandsub_data_cut_all.Scale(normfac_data)
W_pions_normrandsub_data_cut_all.Scale(normfac_data)
Q2_pions_normrandsub_data_cut_all.Scale(normfac_data)
Epsilon_pions_normrandsub_data_cut_all.Scale(normfac_data)
t_pions_normrandsub_data_cut_all.Scale(normfac_data)
ph_q_pions_normrandsub_data_cut_all.Scale(normfac_data)

# Dummy Normalization
H_gtr_beta_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_gtr_xp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_gtr_yp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_gtr_dp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_gtr_p_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_dc_x_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_dc_y_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_dc_xp_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_dc_yp_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_cal_etotnorm_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_cal_etottracknorm_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_cer_npeSum_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
H_RFTime_Dist_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_gtr_beta_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_gtr_xp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_gtr_yp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_gtr_dp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_gtr_p_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_dc_x_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_dc_y_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_dc_xp_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_dc_yp_fp_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_cal_etotnorm_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_cal_etottracknorm_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_hgcer_npeSum_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_hgcer_xAtCer_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_hgcer_yAtCer_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_ngcer_npeSum_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_ngcer_xAtCer_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_ngcer_yAtCer_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_aero_npeSum_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_aero_xAtAero_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_aero_yAtAero_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
MMpi_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
P_RFTime_Dist_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
CTime_ePiCoinTime_ROC2_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
pmiss_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
pmiss_x_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
pmiss_y_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
pmiss_z_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
emiss_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
W_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
Q2_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
Epsilon_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
t_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)
ph_q_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)

# SIMC Normalization
H_hsdelta_pions_simc_cut_all.Scale(normfac_simc)
H_hsxptar_pions_simc_cut_all.Scale(normfac_simc)
H_hsyptar_pions_simc_cut_all.Scale(normfac_simc)
H_hsytar_pions_simc_cut_all.Scale(normfac_simc)
H_hsxfp_pions_simc_cut_all.Scale(normfac_simc)
H_hsyfp_pions_simc_cut_all.Scale(normfac_simc)
H_hsxpfp_pions_simc_cut_all.Scale(normfac_simc)
H_hsypfp_pions_simc_cut_all.Scale(normfac_simc)
P_ssdelta_pions_simc_cut_all.Scale(normfac_simc)
P_ssxptar_pions_simc_cut_all.Scale(normfac_simc)
P_ssyptar_pions_simc_cut_all.Scale(normfac_simc)
P_ssytar_pions_simc_cut_all.Scale(normfac_simc)
P_ssxfp_pions_simc_cut_all.Scale(normfac_simc)
P_ssyfp_pions_simc_cut_all.Scale(normfac_simc)
P_ssxpfp_pions_simc_cut_all.Scale(normfac_simc)
P_ssypfp_pions_simc_cut_all.Scale(normfac_simc)
q_pions_simc_cut_all.Scale(normfac_simc)
nu_pions_simc_cut_all.Scale(normfac_simc)
Q2_pions_simc_cut_all.Scale(normfac_simc)
epsilon_pions_simc_cut_all.Scale(normfac_simc)
thetapq_pions_simc_cut_all.Scale(normfac_simc)
phipq_pions_simc_cut_all.Scale(normfac_simc)
pmiss_pions_simc_cut_all.Scale(normfac_simc)
pmiss_x_pions_simc_cut_all.Scale(normfac_simc)
pmiss_y_pions_simc_cut_all.Scale(normfac_simc)
pmiss_z_pions_simc_cut_all.Scale(normfac_simc)
emiss_pions_simc_cut_all.Scale(normfac_simc)
W_pions_simc_cut_all.Scale(normfac_simc)
t_pions_simc_cut_all.Scale(normfac_simc)
MMpi_pions_simc_cut_all.Scale(normfac_simc)

############################################################################################################################################

# Dummy Subtraction
H_gtr_beta_pions_normdummysub_data_cut_all.Add(H_gtr_beta_pions_normrandsub_data_cut_all, H_gtr_beta_pions_normrandsub_dummy_cut_all, 1, -1)
H_gtr_xp_pions_normdummysub_data_cut_all.Add(H_gtr_xp_pions_normrandsub_data_cut_all, H_gtr_xp_pions_normrandsub_dummy_cut_all, 1, -1)
H_gtr_yp_pions_normdummysub_data_cut_all.Add(H_gtr_yp_pions_normrandsub_data_cut_all, H_gtr_yp_pions_normrandsub_dummy_cut_all, 1, -1)
H_gtr_dp_pions_normdummysub_data_cut_all.Add(H_gtr_dp_pions_normrandsub_data_cut_all, H_gtr_dp_pions_normrandsub_dummy_cut_all, 1, -1)
H_gtr_p_pions_normdummysub_data_cut_all.Add(H_gtr_p_pions_normrandsub_data_cut_all, H_gtr_p_pions_normrandsub_dummy_cut_all, 1, -1)
H_dc_x_fp_pions_normdummysub_data_cut_all.Add(H_dc_x_fp_pions_normrandsub_data_cut_all, H_dc_x_fp_pions_normrandsub_dummy_cut_all, 1, -1)
H_dc_y_fp_pions_normdummysub_data_cut_all.Add(H_dc_y_fp_pions_normrandsub_data_cut_all, H_dc_y_fp_pions_normrandsub_dummy_cut_all, 1, -1)
H_dc_xp_fp_pions_normdummysub_data_cut_all.Add(H_dc_xp_fp_pions_normrandsub_data_cut_all, H_dc_xp_fp_pions_normrandsub_dummy_cut_all, 1, -1)
H_dc_yp_fp_pions_normdummysub_data_cut_all.Add(H_dc_yp_fp_pions_normrandsub_data_cut_all, H_dc_yp_fp_pions_normrandsub_dummy_cut_all, 1, -1)
H_cal_etotnorm_pions_normdummysub_data_cut_all.Add(H_cal_etotnorm_pions_normrandsub_data_cut_all, H_cal_etotnorm_pions_normrandsub_dummy_cut_all, 1, -1)
H_cal_etottracknorm_pions_normdummysub_data_cut_all.Add(H_cal_etottracknorm_pions_normrandsub_data_cut_all, H_cal_etottracknorm_pions_normrandsub_dummy_cut_all, 1, -1)
H_cer_npeSum_pions_normdummysub_data_cut_all.Add(H_cer_npeSum_pions_normrandsub_data_cut_all, H_cer_npeSum_pions_normrandsub_dummy_cut_all, 1, -1)
H_RFTime_Dist_pions_normdummysub_data_cut_all.Add(H_RFTime_Dist_pions_normrandsub_data_cut_all, H_RFTime_Dist_pions_normrandsub_dummy_cut_all, 1, -1)
P_gtr_beta_pions_normdummysub_data_cut_all.Add(P_gtr_beta_pions_normrandsub_data_cut_all, P_gtr_beta_pions_normrandsub_dummy_cut_all, 1, -1)
P_gtr_xp_pions_normdummysub_data_cut_all.Add(P_gtr_xp_pions_normrandsub_data_cut_all, P_gtr_xp_pions_normrandsub_dummy_cut_all, 1, -1)
P_gtr_yp_pions_normdummysub_data_cut_all.Add(P_gtr_yp_pions_normrandsub_data_cut_all, P_gtr_yp_pions_normrandsub_dummy_cut_all, 1, -1)
P_gtr_dp_pions_normdummysub_data_cut_all.Add(P_gtr_dp_pions_normrandsub_data_cut_all, P_gtr_dp_pions_normrandsub_dummy_cut_all, 1, -1)
P_gtr_p_pions_normdummysub_data_cut_all.Add(P_gtr_p_pions_normrandsub_data_cut_all, P_gtr_p_pions_normrandsub_dummy_cut_all, 1, -1)
P_dc_x_fp_pions_normdummysub_data_cut_all.Add(P_dc_x_fp_pions_normrandsub_data_cut_all, P_dc_x_fp_pions_normrandsub_dummy_cut_all, 1, -1)
P_dc_y_fp_pions_normdummysub_data_cut_all.Add(P_dc_y_fp_pions_normrandsub_data_cut_all, P_dc_y_fp_pions_normrandsub_dummy_cut_all, 1, -1)
P_dc_xp_fp_pions_normdummysub_data_cut_all.Add(P_dc_xp_fp_pions_normrandsub_data_cut_all, P_dc_xp_fp_pions_normrandsub_dummy_cut_all, 1, -1)
P_dc_yp_fp_pions_normdummysub_data_cut_all.Add(P_dc_yp_fp_pions_normrandsub_data_cut_all, P_dc_yp_fp_pions_normrandsub_dummy_cut_all, 1, -1)
P_cal_etotnorm_pions_normdummysub_data_cut_all.Add(P_cal_etotnorm_pions_normrandsub_data_cut_all, P_cal_etotnorm_pions_normrandsub_dummy_cut_all, 1, -1)
P_cal_etottracknorm_pions_normdummysub_data_cut_all.Add(P_cal_etottracknorm_pions_normrandsub_data_cut_all, P_cal_etottracknorm_pions_normrandsub_dummy_cut_all, 1, -1)
P_hgcer_npeSum_pions_normdummysub_data_cut_all.Add(P_hgcer_npeSum_pions_normrandsub_data_cut_all, P_hgcer_npeSum_pions_normrandsub_dummy_cut_all, 1, -1)
P_hgcer_xAtCer_pions_normdummysub_data_cut_all.Add(P_hgcer_xAtCer_pions_normrandsub_data_cut_all, P_hgcer_xAtCer_pions_normrandsub_dummy_cut_all, 1, -1)
P_hgcer_yAtCer_pions_normdummysub_data_cut_all.Add(P_hgcer_yAtCer_pions_normrandsub_data_cut_all, P_hgcer_yAtCer_pions_normrandsub_dummy_cut_all, 1, -1)
P_ngcer_npeSum_pions_normdummysub_data_cut_all.Add(P_ngcer_npeSum_pions_normrandsub_data_cut_all, P_ngcer_npeSum_pions_normrandsub_dummy_cut_all, 1, -1)
P_ngcer_xAtCer_pions_normdummysub_data_cut_all.Add(P_ngcer_xAtCer_pions_normrandsub_data_cut_all, P_ngcer_xAtCer_pions_normrandsub_dummy_cut_all, 1, -1)
P_ngcer_yAtCer_pions_normdummysub_data_cut_all.Add(P_ngcer_yAtCer_pions_normrandsub_data_cut_all, P_ngcer_yAtCer_pions_normrandsub_dummy_cut_all, 1, -1)
P_aero_npeSum_pions_normdummysub_data_cut_all.Add(P_aero_npeSum_pions_normrandsub_data_cut_all, P_aero_npeSum_pions_normrandsub_dummy_cut_all, 1, -1)
P_aero_xAtAero_pions_normdummysub_data_cut_all.Add(P_aero_xAtAero_pions_normrandsub_data_cut_all, P_aero_xAtAero_pions_normrandsub_dummy_cut_all, 1, -1)
P_aero_yAtAero_pions_normdummysub_data_cut_all.Add(P_aero_yAtAero_pions_normrandsub_data_cut_all, P_aero_yAtAero_pions_normrandsub_dummy_cut_all, 1, -1)
MMpi_pions_normdummysub_data_cut_all.Add(MMpi_pions_normrandsub_data_cut_all, MMpi_pions_normrandsub_dummy_cut_all, 1, -1)
P_RFTime_Dist_pions_normdummysub_data_cut_all.Add(P_RFTime_Dist_pions_normrandsub_data_cut_all, P_RFTime_Dist_pions_normrandsub_dummy_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_normdummysub_data_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_normrandsub_data_cut_all, CTime_ePiCoinTime_ROC2_pions_normrandsub_dummy_cut_all, 1, -1)
pmiss_pions_normdummysub_data_cut_all.Add(pmiss_pions_normrandsub_data_cut_all, pmiss_pions_normrandsub_dummy_cut_all, 1, -1)
pmiss_x_pions_normdummysub_data_cut_all.Add(pmiss_x_pions_normrandsub_data_cut_all, pmiss_x_pions_normrandsub_dummy_cut_all, 1, -1)
pmiss_y_pions_normdummysub_data_cut_all.Add(pmiss_y_pions_normrandsub_data_cut_all, pmiss_y_pions_normrandsub_dummy_cut_all, 1, -1)
pmiss_z_pions_normdummysub_data_cut_all.Add(pmiss_z_pions_normrandsub_data_cut_all, pmiss_z_pions_normrandsub_dummy_cut_all, 1, -1)
emiss_pions_normdummysub_data_cut_all.Add(emiss_pions_normrandsub_data_cut_all, emiss_pions_normrandsub_dummy_cut_all, 1, -1)
W_pions_normdummysub_data_cut_all.Add(W_pions_normrandsub_data_cut_all, W_pions_normrandsub_dummy_cut_all, 1, -1)
Q2_pions_normdummysub_data_cut_all.Add(Q2_pions_normrandsub_data_cut_all, Q2_pions_normrandsub_dummy_cut_all, 1, -1)
Epsilon_pions_normdummysub_data_cut_all.Add(Epsilon_pions_normrandsub_data_cut_all, Epsilon_pions_normrandsub_dummy_cut_all, 1, -1)
t_pions_normdummysub_data_cut_all.Add(t_pions_normrandsub_data_cut_all, t_pions_normrandsub_dummy_cut_all, 1, -1)
ph_q_pions_normdummysub_data_cut_all.Add(ph_q_pions_normrandsub_data_cut_all, ph_q_pions_normrandsub_dummy_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

###########################################################################################################################################

# Yield Calculations
dN_data_hms_dp = array.array('d', [0.0])
dN_data_shms_dp = array.array('d', [0.0])
dN_data_hms_xp = array.array('d', [0.0])
dN_data_hms_yp = array.array('d', [0.0])
dN_data_shms_xp = array.array('d', [0.0])
dN_data_shms_yp = array.array('d', [0.0])
dN_data_hms_xpfp = array.array('d', [0.0])
dN_data_hms_ypfp = array.array('d', [0.0])
dN_data_shms_xpfp = array.array('d', [0.0])
dN_data_shms_ypfp = array.array('d', [0.0])
dN_data_pmiss = array.array('d', [0.0])
dN_data_emiss = array.array('d', [0.0])
dN_data_W = array.array('d', [0.0])
dN_data_MMpi = array.array('d', [0.0])

dN_simc_hsdelta = array.array('d', [0.0])
dN_simc_ssdelta = array.array('d', [0.0])
dN_simc_hsxptar = array.array('d', [0.0])
dN_simc_hsyptar = array.array('d', [0.0])
dN_simc_ssxptar = array.array('d', [0.0])
dN_simc_ssyptar = array.array('d', [0.0])
dN_simc_hsxpfp = array.array('d', [0.0])
dN_simc_hsypfp = array.array('d', [0.0])
dN_simc_ssxpfp = array.array('d', [0.0])
dN_simc_ssypfp = array.array('d', [0.0])
dN_simc_pmiss = array.array('d', [0.0])
dN_simc_emiss = array.array('d', [0.0])
dN_simc_W = array.array('d', [0.0])
dN_simc_MMpi = array.array('d', [0.0])

N_data_hms_dp = H_gtr_dp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_hms_dp,"")
N_data_shms_dp = P_gtr_dp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_shms_dp,"")
N_data_hms_xp = H_gtr_xp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_hms_xp,"")
N_data_hms_yp = H_gtr_yp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_hms_yp,"")
N_data_shms_xp = P_gtr_xp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_shms_xp,"")
N_data_shms_yp = P_gtr_yp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_shms_yp,"")
N_data_hms_xpfp = H_dc_xp_fp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_hms_xpfp,"")
N_data_hms_ypfp = H_dc_yp_fp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_hms_ypfp,"")
N_data_shms_xpfp = P_dc_xp_fp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_shms_xpfp,"")
N_data_shms_ypfp = P_dc_yp_fp_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_shms_ypfp,"")
N_data_pmiss = pmiss_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_pmiss,"")
N_data_emiss = emiss_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins,dN_data_emiss,"")
N_data_W = W_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins_p,dN_data_W,"")
N_data_MMpi = MMpi_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins_p,dN_data_MMpi,"")

N_simc_hsdelta = H_hsdelta_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_hsdelta,"")
N_simc_ssdelta = P_ssdelta_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_ssdelta,"")
N_simc_hsxptar = H_hsxptar_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_hsxptar,"")
N_simc_hsyptar = H_hsyptar_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_hsyptar,"")
N_simc_ssxptar = P_ssxptar_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_ssxptar,"")
N_simc_ssyptar = P_ssyptar_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_ssyptar,"")
N_simc_hsxpfp = H_hsxpfp_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_hsxpfp,"")
N_simc_hsypfp = H_hsypfp_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_hsypfp,"")
N_simc_ssxpfp = P_ssxpfp_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_ssxpfp,"")
N_simc_ssypfp = P_ssypfp_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_ssypfp,"")
N_simc_pmiss = pmiss_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_pmiss,"")
N_simc_emiss = emiss_pions_simc_cut_all.IntegralAndError(1,nbins,dN_simc_emiss,"")
N_simc_W = W_pions_simc_cut_all.IntegralAndError(1,nbins_p,dN_simc_W,"")
N_simc_MMpi = MMpi_pions_simc_cut_all.IntegralAndError(1,nbins_p,dN_simc_MMpi,"")

# Physics Yield Calculations
physics_yield_hdelta = N_data_hms_dp
physics_yield_pdelta = N_data_shms_dp
physics_yield_hxptar = N_data_hms_xp
physics_yield_hyptar = N_data_hms_yp
physics_yield_pxptar = N_data_shms_xp
physics_yield_pyptar = N_data_shms_yp
physics_yield_hxpfp = N_data_hms_xpfp
physics_yield_hypfp = N_data_hms_ypfp
physics_yield_pxpfp = N_data_shms_xpfp
physics_yield_pypfp = N_data_shms_ypfp
physics_yield_pmiss = N_data_pmiss
physics_yield_emiss = N_data_emiss
physics_yield_W = N_data_W
physics_yield_MMpi = N_data_MMpi
physics_yield = N_data_MMpi


# Physics SIMC Calculations
simc_yield_hdelta = N_simc_hsdelta
simc_yield_pdelta = N_simc_ssdelta
simc_yield_hxptar = N_simc_hsxptar
simc_yield_hyptar = N_simc_hsyptar
simc_yield_pxptar = N_simc_ssxptar
simc_yield_pyptar = N_simc_ssyptar
simc_yield_hxpfp = N_simc_hsxpfp
simc_yield_hypfp = N_simc_hsypfp
simc_yield_pxpfp = N_simc_ssxpfp
simc_yield_pypfp = N_simc_ssypfp
simc_yield_pmiss = N_simc_pmiss
simc_yield_emiss = N_simc_emiss
simc_yield_W = N_simc_W
simc_yield_MMpi = N_simc_MMpi
simc_yield = N_simc_MMpi


# Data/SIMC Ratio Calculations
dataSimcRatio_hdelta = N_data_hms_dp/N_simc_hsdelta
dataSimcRatio_pdelta = N_data_shms_dp/N_simc_ssdelta
dataSimcRatio_hxptar = N_data_hms_xp/N_simc_hsxptar
dataSimcRatio_hyptar = N_data_hms_yp/N_simc_hsyptar
dataSimcRatio_pxptar = N_data_shms_xp/N_simc_ssxptar
dataSimcRatio_pyptar = N_data_shms_yp/N_simc_ssyptar
dataSimcRatio_hxpfp = N_data_hms_xpfp/N_simc_hsxpfp
dataSimcRatio_hypfp = N_data_hms_ypfp/N_simc_hsypfp
dataSimcRatio_pxpfp = N_data_shms_xpfp/N_simc_ssxpfp
dataSimcRatio_pypfp = N_data_shms_ypfp/N_simc_ssypfp
dataSimcRatio_pmiss = N_data_pmiss/N_simc_pmiss
dataSimcRatio_emiss = N_data_emiss/N_simc_emiss
dataSimcRatio_W = N_data_W/N_simc_W
dataSimcRatio_MMpi = N_data_MMpi/N_simc_MMpi
DATASIMC_Ratio = N_data_MMpi/N_simc_MMpi

# Error Propagation
dN_data_error = array.array('d', [0.0])
dN_dummy_error = array.array('d', [0.0])
dN_simc_error = array.array('d', [0.0])

N_data = pmiss_pions_randsub_data_cut_all_error.IntegralAndError(1,nbins,dN_data_error,"")
N_dummy = pmiss_pions_randsub_dummy_cut_all_error.IntegralAndError(1,nbins,dN_dummy_error,"")
N_simc = pmiss_pions_simc_cut_all_error.IntegralAndError(1,nbins,dN_simc_error,"")

# Data/SIMC Ratio Calculations with errors
N_data_norm = N_data / total_data_effective_charge
N_dummy_norm = N_dummy / total_dummy_effective_charge
N_simc_norm = N_simc * normfac_simc
N_data_dummy_sub_norm = N_data_norm - N_dummy_norm

physics_yield_cal = N_data_dummy_sub_norm
simc_yield_cal = N_simc_norm

dN_data_norm = N_data_norm * ma.sqrt((dN_data_error[0]/N_data)**2  + (total_data_effective_charge_error/total_data_effective_charge)**2)
dN_dummy_norm = N_dummy_norm * ma.sqrt((dN_dummy_error[0]/N_dummy)**2 + (total_dummy_effective_charge_error/total_dummy_effective_charge)**2)
dN_data_dummy_sub_norm = ma.sqrt((dN_data_norm)**2 + (dN_dummy_norm)**2)
dN_simc_norm = N_simc_norm * (dN_simc_error[0]/N_simc)

physics_yield_err = dN_data_dummy_sub_norm
simc_yield_err = dN_simc_norm

dataSimcRatio = (N_data_dummy_sub_norm) / N_simc_norm
dataSimcRatio_err = dataSimcRatio * ma.sqrt((dN_data_dummy_sub_norm/N_data_dummy_sub_norm)**2 + (dN_simc_norm/N_simc_norm)**2)

print("="*40)
print("Physics Yield hdelta = {:.3f} +/- {:.3f}".format(physics_yield_hdelta, physics_yield_err))
print("Physics Yield pdelta = {:.3f} +/- {:.3f}".format(physics_yield_pdelta, physics_yield_err))
print("Physics Yield hxptar = {:.3f} +/- {:.3f}".format(physics_yield_hxptar, physics_yield_err))
print("Physics Yield hyptar = {:.3f} +/- {:.3f}".format(physics_yield_hyptar, physics_yield_err))
print("Physics Yield pxptar = {:.3f} +/- {:.3f}".format(physics_yield_pxptar, physics_yield_err))
print("Physics Yield pyptar = {:.3f} +/- {:.3f}".format(physics_yield_pyptar, physics_yield_err))
print("Physics Yield hxpfp = {:.3f} +/- {:.3f}".format(physics_yield_hxpfp, physics_yield_err))
print("Physics Yield hypfp = {:.3f} +/- {:.3f}".format(physics_yield_hypfp, physics_yield_err))
print("Physics Yield pxpfp = {:.3f} +/- {:.3f}".format(physics_yield_pxpfp, physics_yield_err))
print("Physics Yield pypfp = {:.3f} +/- {:.3f}".format(physics_yield_pypfp, physics_yield_err))
print("Physics Yield pmiss = {:.3f} +/- {:.3f}".format(physics_yield_pmiss, physics_yield_err))
print("Physics Yield emiss = {:.3f} +/- {:.3f}".format(physics_yield_emiss, physics_yield_err))
print("Physics Yield W = {:.3f} +/- {:.3f}".format(physics_yield_W, physics_yield_err))
print("Physics Yield MMpi = {:.3f} +/- {:.3f}".format(physics_yield_MMpi, physics_yield_err))
print("="*40)
print("SIMC Yield hdelta = {:.3f} +/- {:.3f}".format(simc_yield_hdelta, simc_yield_err))
print("SIMC Yield pdelta = {:.3f} +/- {:.3f}".format(simc_yield_pdelta, simc_yield_err))
print("SIMC Yield hxptar = {:.3f} +/- {:.3f}".format(simc_yield_hxptar, simc_yield_err))
print("SIMC Yield hyptar = {:.3f} +/- {:.3f}".format(simc_yield_hyptar, simc_yield_err))
print("SIMC Yield pxptar = {:.3f} +/- {:.3f}".format(simc_yield_pxptar, simc_yield_err))
print("SIMC Yield pyptar = {:.3f} +/- {:.3f}".format(simc_yield_pyptar, simc_yield_err))
print("SIMC Yield hxpfp = {:.3f} +/- {:.3f}".format(simc_yield_hxpfp, simc_yield_err))
print("SIMC Yield hypfp = {:.3f} +/- {:.3f}".format(simc_yield_hypfp, simc_yield_err))
print("SIMC Yield pxpfp = {:.3f} +/- {:.3f}".format(simc_yield_pxpfp, simc_yield_err))
print("SIMC Yield pypfp = {:.3f} +/- {:.3f}".format(simc_yield_pypfp, simc_yield_err))
print("SIMC Yield pmiss = {:.3f} +/- {:.3f}".format(simc_yield_pmiss, simc_yield_err))
print("SIMC Yield emiss = {:.3f} +/- {:.3f}".format(simc_yield_emiss, simc_yield_err))
print("SIMC Yield W = {:.3f} +/- {:.3f}".format(simc_yield_W, simc_yield_err))
print("SIMC Yield MMpi = {:.3f} +/- {:.3f}".format(simc_yield, simc_yield_err))
print("="*40)
print("Data/SIMC ratio hdelta = {:.3f} +/- {:.3f}".format(dataSimcRatio_hdelta, dataSimcRatio_err))
print("Data/SIMC ratio pdelta = {:.3f} +/- {:.3f}".format(dataSimcRatio_pdelta, dataSimcRatio_err))
print("Data/SIMC ratio hxptar = {:.3f} +/- {:.3f}".format(dataSimcRatio_hxptar, dataSimcRatio_err))
print("Data/SIMC ratio hyptar = {:.3f} +/- {:.3f}".format(dataSimcRatio_hyptar, dataSimcRatio_err))
print("Data/SIMC ratio pxptar = {:.3f} +/- {:.3f}".format(dataSimcRatio_pxptar, dataSimcRatio_err))
print("Data/SIMC ratio pyptar = {:.3f} +/- {:.3f}".format(dataSimcRatio_pyptar, dataSimcRatio_err))
print("Data/SIMC ratio hxpfp = {:.3f} +/- {:.3f}".format(dataSimcRatio_hxpfp, dataSimcRatio_err))
print("Data/SIMC ratio hypfp = {:.3f} +/- {:.3f}".format(dataSimcRatio_hypfp, dataSimcRatio_err))
print("Data/SIMC ratio pxpfp = {:.3f} +/- {:.3f}".format(dataSimcRatio_pxpfp, dataSimcRatio_err))
print("Data/SIMC ratio pypfp = {:.3f} +/- {:.3f}".format(dataSimcRatio_pypfp, dataSimcRatio_err))
print("Data/SIMC ratio pmiss = {:.3f} +/- {:.3f}".format(dataSimcRatio_pmiss, dataSimcRatio_err))
print("Data/SIMC ratio emiss = {:.3f} +/- {:.3f}".format(dataSimcRatio_emiss, dataSimcRatio_err))
print("Data/SIMC ratio W = {:.3f} +/- {:.3f}".format(dataSimcRatio_W, dataSimcRatio_err))
print("Data/SIMC ratio MMpi = {:.3f} +/- {:.3f}".format(dataSimcRatio_MMpi, dataSimcRatio_err))
print("="*40)

# Define the output CSV file name
yield_cal_output_csv_path = "%s/LTSEP_ANALYSIS/LTSep_CSVs/datasimc_ratios_csv/%s/%s_Physics_DataSIMC_Ratio_Calculation.csv" % (REPLAYPATH, physet_dir_name, setting_name)

# Prepare the header for the CSV file
header = [
    "Physics Setting", "Physics Yield", "Physics Yield Error",
    "SIMC Yield", "SIMC Yield Error", "Data/SIMC Ratio", "Data/SIMC Ratio Error"
]

# Prepare data row as a dictionary
new_row = {
    "Physics Setting": PHY_SETTING,
    "Physics Yield": f"{physics_yield}",
    "Physics Yield Error": f"{physics_yield_err}",
    "SIMC Yield": f"{simc_yield}",
    "SIMC Yield Error": f"{simc_yield_err}",
    "Data/SIMC Ratio": f"{DATASIMC_Ratio}",
    "Data/SIMC Ratio Error": f"{dataSimcRatio_err}",
}

# Read existing rows if file exists
existing_rows = []
if os.path.exists(yield_cal_output_csv_path):
    with open(yield_cal_output_csv_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        existing_rows = list(csv_reader)

# Update existing row or append new
row_found = False
for row in existing_rows:
    if row["Physics Setting"] == PHY_SETTING:
        row.update(new_row)
        row_found = True
        break

if not row_found:
    existing_rows.append(new_row)

# Write back all rows to CSV
with open(yield_cal_output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(existing_rows)  # Write all rows including updated or new one

print(f"DATA/SIMC Ratios written to {yield_cal_output_csv_path}")

#############################################################################################################################################

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Format the dynamic strings first
phy_setting_text = "{}".format(PHY_SETTING)
Beam_Energy_S, HMS_p, HMS_theta, SHMS_p, SHMS_theta  = filtered_data_df['Beam_Energy'].iloc[0],filtered_data_df['HMS_P_Central'].iloc[0],filtered_data_df['HMS_Angle'].iloc[0],filtered_data_df['SHMS_P_Central'].iloc[0],filtered_data_df['SHMS_Angle'].iloc[0]
ratio_text = "DATA/SIMC Ratio = {:.3f} ± {:.3f}".format(DATASIMC_Ratio, dataSimcRatio_err)

# Scale factor for SIMC plots
#match_scale = 0

# Saving histograms in PDF
c1_kin = TCanvas("c1_kin", "Variables Distributions", 100, 0, 1200, 1200)
c1_kin.Divide(2,2)
c1_kin.cd(1)
c1_kin_text_lines = [
#    ROOT.TText(0.5, 0.9, "Pion ProdCoin Setting"),
    ROOT.TText(0.5, 0.9, phy_setting_text),
    ROOT.TText(0.5, 0.8, 'Beam Energy = {:.3f}'.format(Beam_Energy_S)),
    ROOT.TText(0.5, 0.7, 'HMS_p = {:.3f}'.format(HMS_p)),
    ROOT.TText(0.5, 0.6, 'HMS_theta = {:.3f}'.format(HMS_theta)),
    ROOT.TText(0.5, 0.5, 'SHMS_p = {:.3f}'.format(SHMS_p)),
    ROOT.TText(0.5, 0.4, 'SHMS_theta = {:.3f}'.format(SHMS_theta)),
    ROOT.TText(0.5, 0.3, "Red = SIMC"),
    ROOT.TText(0.5, 0.2, "Blue = DATA"),
    ROOT.TText(0.5, 0.1, ratio_text)

]
for c1_kin_text in c1_kin_text_lines:
    c1_kin_text.SetTextSize(0.053)
    c1_kin_text.SetTextAlign(22)
    c1_kin_text.SetTextColor(ROOT.kGreen + 4)
    if c1_kin_text.GetTitle() == "Red = SIMC":
       c1_kin_text.SetTextColor(ROOT.kRed)  # Setting text color to red
    if c1_kin_text.GetTitle() == "Blue = DATA":
       c1_kin_text.SetTextColor(ROOT.kBlue)  # Setting text color to red
    c1_kin_text.Draw()

c1_kin.cd(2)
MMpi_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.88, 1.07)
MMpi_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
MMpi_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
MMpi_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
MMpi_pions_normdummysub_data_cut_all.Draw("E1")
#MMpi_pions_simc_cut_all.Scale(match_scale)
MMpi_pions_simc_cut_all.GetXaxis().SetRangeUser(0.88, 1.07)
MMpi_pions_simc_cut_all.GetYaxis().CenterTitle()
MMpi_pions_simc_cut_all.SetLineColor(ROOT.kRed)
MMpi_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
MMpi_pions_simc_cut_all.Draw("E1same")
# Add legend if needed
legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend.AddEntry(MMpi_pions_normdummysub_data_cut_all, "DATA", "lep")
legend.AddEntry(MMpi_pions_simc_cut_all, "SIMC", "lep")
legend.Draw()
c1_kin.cd(3)
W_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(2.4, 2.8)
W_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
W_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
W_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
W_pions_normdummysub_data_cut_all.Draw("E1")
#W_pions_simc_cut_all.Scale(match_scale)
W_pions_simc_cut_all.GetXaxis().SetRangeUser(2.4, 2.8)
W_pions_simc_cut_all.GetYaxis().CenterTitle()
W_pions_simc_cut_all.SetLineColor(ROOT.kRed)
W_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
W_pions_simc_cut_all.Draw("E1same")
# Legend for W plot
legend_W = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_W.AddEntry(W_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_W.AddEntry(W_pions_simc_cut_all, "SIMC", "lep")
legend_W.Draw()
c1_kin.cd(4)
t_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.1, 0.7)
t_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
t_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
t_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
t_pions_normdummysub_data_cut_all.Draw("E1")
#t_pions_simc_cut_all.Scale(match_scale)
t_pions_simc_cut_all.GetXaxis().SetRangeUser(0.1, 0.7)
t_pions_simc_cut_all.GetYaxis().CenterTitle()
t_pions_simc_cut_all.SetLineColor(ROOT.kRed)
t_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
t_pions_simc_cut_all.Draw("E1same")
# Legend for t plot
legend_t = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_t.AddEntry(t_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_t.AddEntry(t_pions_simc_cut_all, "SIMC", "lep")
legend_t.Draw()
c1_kin.Print(Pion_Analysis_Distributions + '(')

c2_kin = TCanvas("c2_kin", "Variables Distributions", 100, 0, 1200, 1200)
c2_kin.Divide(2,2)
c2_kin.cd(1)
Q2_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(3.0, 4.6)
Q2_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
Q2_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
Q2_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
Q2_pions_normdummysub_data_cut_all.Draw("E1")
#Q2_pions_simc_cut_all.Scale(match_scale)
Q2_pions_simc_cut_all.GetXaxis().SetRangeUser(3.0, 4.6)
Q2_pions_simc_cut_all.GetYaxis().CenterTitle()
Q2_pions_simc_cut_all.SetLineColor(ROOT.kRed)
Q2_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
Q2_pions_simc_cut_all.Draw("E1same")
# Legend for Q² plot
legend_Q2 = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_Q2.AddEntry(Q2_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_Q2.AddEntry(Q2_pions_simc_cut_all, "SIMC", "lep")
legend_Q2.Draw()
c2_kin.cd(2)
ph_q_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.0, 360.0)
#ph_q_pions_normdummysub_data_cut_all.GetYaxis().SetRangeUser(-0.06, 0.1)
ph_q_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
ph_q_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
ph_q_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
ph_q_pions_normdummysub_data_cut_all.Draw("E1")
#phipq_pions_simc_cut_all.Scale(match_scale)
phipq_pions_simc_cut_all.GetXaxis().SetRangeUser(0.0, 360.0)
#phipq_pions_simc_cut_all.GetYaxis().SetRangeUser(-0.06, 0.1)
phipq_pions_simc_cut_all.GetYaxis().CenterTitle()
phipq_pions_simc_cut_all.SetLineColor(ROOT.kRed)
phipq_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
phipq_pions_simc_cut_all.Draw("E1same")
# Legend for phi Plot
legend_phipq = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_phipq.AddEntry(ph_q_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_phipq.AddEntry(phipq_pions_simc_cut_all, "SIMC", "lep")
legend_phipq.Draw()
c2_kin.cd(3)
H_gtr_dp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-9, 9)
#H_gtr_dp_pions_normdummysub_data_cut_all.GetYaxis().SetRangeUser(0, 0.08)
H_gtr_dp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_gtr_dp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_gtr_dp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_gtr_dp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsdelta_pions_simc_cut_all.Scale(match_scale)
H_hsdelta_pions_simc_cut_all.GetXaxis().SetRangeUser(-9, 9)
#H_hsdelta_pions_simc_cut_all.GetYaxis().SetRangeUser(0, 0.08)
H_hsdelta_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsdelta_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsdelta_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsdelta_pions_simc_cut_all.Draw("E1same")
# Legend for HMS Delta Plot
legend_hsdelta = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsdelta.AddEntry(H_gtr_dp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsdelta.AddEntry(H_hsdelta_pions_simc_cut_all, "SIMC", "lep")
legend_hsdelta.Draw()
c2_kin.cd(4)
P_gtr_dp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-10, 6.0)
P_gtr_dp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_gtr_dp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_gtr_dp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_gtr_dp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssdelta_pions_simc_cut_all.Scale(match_scale)
P_ssdelta_pions_simc_cut_all.GetXaxis().SetRangeUser(-10, 6.0)
P_ssdelta_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssdelta_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssdelta_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssdelta_pions_simc_cut_all.Draw("E1same")
# Legend for SS Delta Plot
legend_ssdelta = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssdelta.AddEntry(P_gtr_dp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssdelta.AddEntry(P_ssdelta_pions_simc_cut_all, "SIMC", "lep")
legend_ssdelta.Draw()
c2_kin.Print(Pion_Analysis_Distributions )

c3_kin = TCanvas("c3_kin", "Variables Distributions", 100, 0, 1200, 1200)
c3_kin.Divide(2,2)
c3_kin.cd(1)
H_gtr_xp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
H_gtr_xp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_gtr_xp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_gtr_xp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_gtr_xp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsxptar_pions_simc_cut_all.Scale(match_scale)
H_hsxptar_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
H_hsxptar_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsxptar_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsxptar_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsxptar_pions_simc_cut_all.Draw("E1same")
# Legend for HMS Xp Plot
legend_hsxptar = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsxptar.AddEntry(H_gtr_xp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsxptar.AddEntry(H_hsxptar_pions_simc_cut_all, "SIMC", "lep")
legend_hsxptar.Draw()
c3_kin.cd(2)
P_gtr_xp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_gtr_xp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_gtr_xp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_gtr_xp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_gtr_xp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssxptar_pions_simc_cut_all.Scale(match_scale)
P_ssxptar_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_ssxptar_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssxptar_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssxptar_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssxptar_pions_simc_cut_all.Draw("E1same")
# Legend for P Xp Plot
legend_ssxptar = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssxptar.AddEntry(P_gtr_xp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssxptar.AddEntry(P_ssxptar_pions_simc_cut_all, "SIMC", "lep")
legend_ssxptar.Draw()
c3_kin.cd(3)
H_gtr_yp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
H_gtr_yp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_gtr_yp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_gtr_yp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_gtr_yp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsyptar_pions_simc_cut_all.Scale(match_scale)
H_hsyptar_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
H_hsyptar_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsyptar_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsyptar_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsyptar_pions_simc_cut_all.Draw("E1same")
# Legend for HMS Yp Plot
legend_hsyptar = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsyptar.AddEntry(H_gtr_yp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsyptar.AddEntry(H_hsyptar_pions_simc_cut_all, "SIMC", "lep")
legend_hsyptar.Draw()
c3_kin.cd(4)
P_gtr_yp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
P_gtr_yp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_gtr_yp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_gtr_yp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_gtr_yp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssyptar_pions_simc_cut_all.Scale(match_scale)
P_ssyptar_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
P_ssyptar_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssyptar_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssyptar_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssyptar_pions_simc_cut_all.Draw("E1same")
# Legend for P Yp Plot
legend_ssyptar = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssyptar.AddEntry(P_gtr_yp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssyptar.AddEntry(P_ssyptar_pions_simc_cut_all, "SIMC", "lep")
legend_ssyptar.Draw()
c3_kin.Print(Pion_Analysis_Distributions)

c4_kin = TCanvas("c4_kin", "Variables Distributions", 100, 0, 1200, 1200)
c4_kin.Divide(2,2)
c4_kin.cd(1)
pmiss_y_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.4, 0.4)
pmiss_y_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
pmiss_y_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
pmiss_y_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
pmiss_y_pions_normdummysub_data_cut_all.Draw("E1")
#pmiss_y_pions_simc_cut_all.Scale(match_scale)
pmiss_y_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.4, 0.4)
pmiss_y_pions_simc_cut_all.GetYaxis().CenterTitle()
pmiss_y_pions_simc_cut_all.SetLineColor(ROOT.kRed)
pmiss_y_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
pmiss_y_pions_simc_cut_all.Draw("E1same")
# Legend for Pmiss Y Plot
legend_pmiss_y = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_pmiss_y.AddEntry(pmiss_y_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_pmiss_y.AddEntry(pmiss_y_pions_simc_cut_all, "SIMC", "lep")
legend_pmiss_y.Draw()
c4_kin.cd(2)
pmiss_z_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-1, -0.2)
pmiss_z_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
pmiss_z_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
pmiss_z_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
pmiss_z_pions_normdummysub_data_cut_all.Draw("E1")
#pmiss_z_pions_simc_cut_all.Scale(match_scale)
pmiss_z_pions_simc_cut_all.GetXaxis().SetRangeUser(-1, -0.2)
pmiss_z_pions_simc_cut_all.GetYaxis().CenterTitle()
pmiss_z_pions_simc_cut_all.SetLineColor(ROOT.kRed)
pmiss_z_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
pmiss_z_pions_simc_cut_all.Draw("E1same")
# Legend for Pmiss Z Plot
legend_pmiss_z = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_pmiss_z.AddEntry(pmiss_z_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_pmiss_z.AddEntry(pmiss_z_pions_simc_cut_all, "SIMC", "lep")
legend_pmiss_z.Draw()
c4_kin.cd(3)
pmiss_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.3, 1.0)
pmiss_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
pmiss_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
pmiss_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
pmiss_pions_normdummysub_data_cut_all.Draw("E1")
#pmiss_pions_simc_cut_all.Scale(match_scale)
pmiss_pions_simc_cut_all.GetXaxis().SetRangeUser(0.3, 1.0)
pmiss_pions_simc_cut_all.GetYaxis().CenterTitle()
pmiss_pions_simc_cut_all.SetLineColor(ROOT.kRed)
pmiss_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
pmiss_pions_simc_cut_all.Draw("E1same")
# Legend for Pmiss Plot
legend_pmiss = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_pmiss.AddEntry(pmiss_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_pmiss.AddEntry(pmiss_pions_simc_cut_all, "SIMC", "lep")
legend_pmiss.Draw()
c4_kin.cd(4)
emiss_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.9, 1.5)
emiss_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
emiss_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
emiss_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
emiss_pions_normdummysub_data_cut_all.Draw("E1")
#emiss_pions_simc_cut_all.Scale(match_scale)
emiss_pions_simc_cut_all.GetXaxis().SetRangeUser(0.9, 1.5)
emiss_pions_simc_cut_all.GetYaxis().CenterTitle()
emiss_pions_simc_cut_all.SetLineColor(ROOT.kRed)
emiss_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
emiss_pions_simc_cut_all.Draw("E1same")
# Legend for Emiss Plot
legend_emiss = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_emiss.AddEntry(emiss_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_emiss.AddEntry(emiss_pions_simc_cut_all, "SIMC", "lep")
legend_emiss.Draw()
c4_kin.Print(Pion_Analysis_Distributions)

c5_kin = TCanvas("c5_kin", "Variables Distributions", 100, 0, 1200, 1200)
c5_kin.Divide(2,2)
c5_kin.cd(1)
H_dc_x_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-50, 50)
H_dc_x_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_dc_x_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_dc_x_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_dc_x_fp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsxfp_pions_simc_cut_all.Scale(match_scale)
H_hsxfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-50, 50)
H_hsxfp_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsxfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsxfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsxfp_pions_simc_cut_all.Draw("E1same")
# Legend for HMS Xfp Plot
legend_hsxfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsxfp.AddEntry(H_dc_x_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsxfp.AddEntry(H_hsxfp_pions_simc_cut_all, "SIMC", "lep")
legend_hsxfp.Draw()
c5_kin.cd(2)
H_dc_y_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-30, 30)
H_dc_y_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_dc_y_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_dc_y_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_dc_y_fp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsyfp_pions_simc_cut_all.Scale(match_scale)
H_hsyfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-30, 30)
H_hsyfp_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsyfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsyfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsyfp_pions_simc_cut_all.Draw("E1same")
# Legend for HMS Yfp Plot
legend_hsyfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsyfp.AddEntry(H_dc_y_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsyfp.AddEntry(H_hsyfp_pions_simc_cut_all, "SIMC", "lep")
legend_hsyfp.Draw()
c5_kin.cd(3)
H_dc_xp_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
H_dc_xp_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_dc_xp_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_dc_xp_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_dc_xp_fp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsxpfp_pions_simc_cut_all.Scale(match_scale)
H_hsxpfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
H_hsxpfp_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsxpfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsxpfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsxpfp_pions_simc_cut_all.Draw("E1same")
# Legend for HMS XpFp Plot
legend_hsxpfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsxpfp.AddEntry(H_dc_xp_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsxpfp.AddEntry(H_hsxpfp_pions_simc_cut_all, "SIMC", "lep")
legend_hsxpfp.Draw()
c5_kin.cd(4)
H_dc_yp_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
H_dc_yp_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
H_dc_yp_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
H_dc_yp_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
H_dc_yp_fp_pions_normdummysub_data_cut_all.Draw("E1")
#H_hsypfp_pions_simc_cut_all.Scale(match_scale)
H_hsypfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.05, 0.05)
H_hsypfp_pions_simc_cut_all.GetYaxis().CenterTitle()
H_hsypfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
H_hsypfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
H_hsypfp_pions_simc_cut_all.Draw("E1same")
# Legend for HMS YpFp Plot
legend_hsyfpfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_hsyfpfp.AddEntry(H_dc_yp_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_hsyfpfp.AddEntry(H_hsypfp_pions_simc_cut_all, "SIMC", "lep")
legend_hsyfpfp.Draw()
c5_kin.Print(Pion_Analysis_Distributions)

c6_kin = TCanvas("c6_kin", "Variables Distributions", 100, 0, 1200, 1200)
c6_kin.Divide(2,2)
c6_kin.cd(1)
P_dc_x_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-20, 20)
P_dc_x_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_dc_x_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_dc_x_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_dc_x_fp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssxfp_pions_simc_cut_all.Scale(match_scale)
P_ssxfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-20, 20)
P_ssxfp_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssxfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssxfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssxfp_pions_simc_cut_all.Draw("E1same")
# Legend for P Xfp Plot
legend_ssxfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssxfp.AddEntry(P_dc_x_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssxfp.AddEntry(P_ssxfp_pions_simc_cut_all, "SIMC", "lep")
legend_ssxfp.Draw()
c6_kin.cd(2)
P_dc_y_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-30, 30)
P_dc_y_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_dc_y_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_dc_y_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_dc_y_fp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssyfp_pions_simc_cut_all.Scale(match_scale)
P_ssyfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-30, 30)
P_ssyfp_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssyfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssyfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssyfp_pions_simc_cut_all.Draw("E1same")
# Legend for P Yfp Plot
legend_ssyfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssyfp.AddEntry(P_dc_y_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssyfp.AddEntry(P_ssyfp_pions_simc_cut_all, "SIMC", "lep")
legend_ssyfp.Draw()
c6_kin.cd(3)
P_dc_xp_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_dc_xp_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_dc_xp_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_dc_xp_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_dc_xp_fp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssxpfp_pions_simc_cut_all.Scale(match_scale)
P_ssxpfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_ssxpfp_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssxpfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssxpfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssxpfp_pions_simc_cut_all.Draw("E1same")
# Legend for P XpFp Plot
legend_ssxpfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssxpfp.AddEntry(P_dc_xp_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssxpfp.AddEntry(P_ssxpfp_pions_simc_cut_all, "SIMC", "lep")
legend_ssxpfp.Draw()
c6_kin.cd(4)
P_dc_yp_fp_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_dc_yp_fp_pions_normdummysub_data_cut_all.GetYaxis().CenterTitle()
P_dc_yp_fp_pions_normdummysub_data_cut_all.SetLineColor(ROOT.kBlue)
P_dc_yp_fp_pions_normdummysub_data_cut_all.SetMarkerColor(ROOT.kBlue)
P_dc_yp_fp_pions_normdummysub_data_cut_all.Draw("E1")
#P_ssypfp_pions_simc_cut_all.Scale(match_scale)
P_ssypfp_pions_simc_cut_all.GetXaxis().SetRangeUser(-0.1, 0.1)
P_ssypfp_pions_simc_cut_all.GetYaxis().CenterTitle()
P_ssypfp_pions_simc_cut_all.SetLineColor(ROOT.kRed)
P_ssypfp_pions_simc_cut_all.SetMarkerColor(ROOT.kRed)
P_ssypfp_pions_simc_cut_all.Draw("E1same")
# Legend for P YpFp Plot
legend_ssypfp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
legend_ssypfp.AddEntry(P_dc_yp_fp_pions_normdummysub_data_cut_all, "DATA", "lep")
legend_ssypfp.AddEntry(P_ssypfp_pions_simc_cut_all, "SIMC", "lep")
legend_ssypfp.Draw()
c6_kin.Print(Pion_Analysis_Distributions + ')')

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/OUTPUT/Analysis/PionLT/ratios/%s/%s_%s_Ratio_Comp_Data.root" % (VOLATILEPATH, physet_dir_name, PHY_SETTING, MaxEvent) , "RECREATE")
d_Cut_Pion_Events_Accpt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Accpt_Data_Cut_All")
d_Cut_Pion_Events_Prompt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Prompt_Data_Cut_All")
d_Cut_Pion_Events_Random_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Random_Data_Cut_All")
d_Cut_Pion_Events_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_RandomSub_Data")
#d_Cut_Pion_Events_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_DummySub_RandomSub_Data")
d_Cut_Pion_Events_Norm_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_Norm_DummySub_RandomSub_Data")
d_Cut_Pion_Events_Norm_SIMC_Data = outHistFile.mkdir("Cut_Pion_Events_Norm_SIMC_Data")


# Writing Histograms for pions
d_Cut_Pion_Events_Accpt_Data_Cut_All.cd()
H_gtr_beta_pions_data_accpt_cut_all.Write()
H_gtr_xp_pions_data_accpt_cut_all.Write()
H_gtr_yp_pions_data_accpt_cut_all.Write()
H_gtr_dp_pions_data_accpt_cut_all.Write()
H_gtr_p_pions_data_accpt_cut_all.Write()
H_dc_x_fp_pions_data_accpt_cut_all.Write()
H_dc_y_fp_pions_data_accpt_cut_all.Write()
H_dc_xp_fp_pions_data_accpt_cut_all.Write()
H_dc_yp_fp_pions_data_accpt_cut_all.Write()
H_cal_etotnorm_pions_data_accpt_cut_all.Write()
H_cal_etottracknorm_pions_data_accpt_cut_all.Write()
H_cer_npeSum_pions_data_accpt_cut_all.Write()
H_RFTime_Dist_pions_data_accpt_cut_all.Write()
P_gtr_beta_pions_data_accpt_cut_all.Write()
P_gtr_xp_pions_data_accpt_cut_all.Write()
P_gtr_yp_pions_data_accpt_cut_all.Write()
P_gtr_dp_pions_data_accpt_cut_all.Write()
P_gtr_p_pions_data_accpt_cut_all.Write()
P_dc_x_fp_pions_data_accpt_cut_all.Write()
P_dc_y_fp_pions_data_accpt_cut_all.Write()
P_dc_xp_fp_pions_data_accpt_cut_all.Write()
P_dc_yp_fp_pions_data_accpt_cut_all.Write()
P_cal_etotnorm_pions_data_accpt_cut_all.Write()
P_cal_etottracknorm_pions_data_accpt_cut_all.Write()
P_hgcer_npeSum_pions_data_accpt_cut_all.Write()
P_hgcer_xAtCer_pions_data_accpt_cut_all.Write()
P_hgcer_yAtCer_pions_data_accpt_cut_all.Write()
P_ngcer_npeSum_pions_data_accpt_cut_all.Write()
P_ngcer_xAtCer_pions_data_accpt_cut_all.Write()
P_ngcer_yAtCer_pions_data_accpt_cut_all.Write()
P_aero_npeSum_pions_data_accpt_cut_all.Write()
P_aero_xAtAero_pions_data_accpt_cut_all.Write()
P_aero_yAtAero_pions_data_accpt_cut_all.Write()
MMpi_pions_data_accpt_cut_all.Write()
P_RFTime_Dist_pions_data_accpt_cut_all.Write()
CTime_ePiCoinTime_ROC2_pions_data_accpt_cut_all.Write()
pmiss_pions_data_accpt_cut_all.Write()
pmiss_x_pions_data_accpt_cut_all.Write()
pmiss_y_pions_data_accpt_cut_all.Write()
pmiss_z_pions_data_accpt_cut_all.Write()
emiss_pions_data_accpt_cut_all.Write()
W_pions_data_accpt_cut_all.Write()
Q2_pions_data_accpt_cut_all.Write()
Epsilon_pions_data_accpt_cut_all.Write()
t_pions_data_accpt_cut_all.Write()
ph_q_pions_data_accpt_cut_all.Write()


d_Cut_Pion_Events_Prompt_Data_Cut_All.cd()
H_gtr_beta_pions_data_prompt_cut_all.Write()
H_gtr_xp_pions_data_prompt_cut_all.Write()
H_gtr_yp_pions_data_prompt_cut_all.Write()
H_gtr_dp_pions_data_prompt_cut_all.Write()
H_gtr_p_pions_data_prompt_cut_all.Write()
H_dc_x_fp_pions_data_prompt_cut_all.Write()
H_dc_y_fp_pions_data_prompt_cut_all.Write()
H_dc_xp_fp_pions_data_prompt_cut_all.Write()
H_dc_yp_fp_pions_data_prompt_cut_all.Write()
H_cal_etotnorm_pions_data_prompt_cut_all.Write()
H_cal_etottracknorm_pions_data_prompt_cut_all.Write()
H_cer_npeSum_pions_data_prompt_cut_all.Write()
H_RFTime_Dist_pions_data_prompt_cut_all.Write()
P_gtr_beta_pions_data_prompt_cut_all.Write()
P_gtr_xp_pions_data_prompt_cut_all.Write()
P_gtr_yp_pions_data_prompt_cut_all.Write()
P_gtr_dp_pions_data_prompt_cut_all.Write()
P_gtr_p_pions_data_prompt_cut_all.Write()
P_dc_x_fp_pions_data_prompt_cut_all.Write()
P_dc_y_fp_pions_data_prompt_cut_all.Write()
P_dc_xp_fp_pions_data_prompt_cut_all.Write()
P_dc_yp_fp_pions_data_prompt_cut_all.Write()
P_cal_etotnorm_pions_data_prompt_cut_all.Write()
P_cal_etottracknorm_pions_data_prompt_cut_all.Write()
P_hgcer_npeSum_pions_data_prompt_cut_all.Write()
P_hgcer_xAtCer_pions_data_prompt_cut_all.Write()
P_hgcer_yAtCer_pions_data_prompt_cut_all.Write()
P_ngcer_npeSum_pions_data_prompt_cut_all.Write()
P_ngcer_xAtCer_pions_data_prompt_cut_all.Write()
P_ngcer_yAtCer_pions_data_prompt_cut_all.Write()
P_aero_npeSum_pions_data_prompt_cut_all.Write()
P_aero_xAtAero_pions_data_prompt_cut_all.Write()
P_aero_yAtAero_pions_data_prompt_cut_all.Write()
MMpi_pions_data_prompt_cut_all.Write()
P_RFTime_Dist_pions_data_prompt_cut_all.Write()
CTime_ePiCoinTime_ROC2_pions_data_prompt_cut_all.Write()
pmiss_pions_data_prompt_cut_all.Write()
pmiss_x_pions_data_prompt_cut_all.Write()
pmiss_y_pions_data_prompt_cut_all.Write()
pmiss_z_pions_data_prompt_cut_all.Write()
emiss_pions_data_prompt_cut_all.Write()
W_pions_data_prompt_cut_all.Write()
Q2_pions_data_prompt_cut_all.Write()
Epsilon_pions_data_prompt_cut_all.Write()
t_pions_data_prompt_cut_all.Write()
ph_q_pions_data_prompt_cut_all.Write()

d_Cut_Pion_Events_Random_Data_Cut_All.cd()
H_gtr_beta_pions_data_random_cut_all.Write()
H_gtr_xp_pions_data_random_cut_all.Write()
H_gtr_yp_pions_data_random_cut_all.Write()
H_gtr_dp_pions_data_random_cut_all.Write()
H_gtr_p_pions_data_random_cut_all.Write()
H_dc_x_fp_pions_data_random_cut_all.Write()
H_dc_y_fp_pions_data_random_cut_all.Write()
H_dc_xp_fp_pions_data_random_cut_all.Write()
H_dc_yp_fp_pions_data_random_cut_all.Write()
H_cal_etotnorm_pions_data_random_cut_all.Write()
H_cal_etottracknorm_pions_data_random_cut_all.Write()
H_cer_npeSum_pions_data_random_cut_all.Write()
H_RFTime_Dist_pions_data_random_cut_all.Write()
P_gtr_beta_pions_data_random_cut_all.Write()
P_gtr_xp_pions_data_random_cut_all.Write()
P_gtr_yp_pions_data_random_cut_all.Write()
P_gtr_dp_pions_data_random_cut_all.Write()
P_gtr_p_pions_data_random_cut_all.Write()
P_dc_x_fp_pions_data_random_cut_all.Write()
P_dc_y_fp_pions_data_random_cut_all.Write()
P_dc_xp_fp_pions_data_random_cut_all.Write()
P_dc_yp_fp_pions_data_random_cut_all.Write()
P_cal_etotnorm_pions_data_random_cut_all.Write()
P_cal_etottracknorm_pions_data_random_cut_all.Write()
P_hgcer_npeSum_pions_data_random_cut_all.Write()
P_hgcer_xAtCer_pions_data_random_cut_all.Write()
P_hgcer_yAtCer_pions_data_random_cut_all.Write()
P_ngcer_npeSum_pions_data_random_cut_all.Write()
P_ngcer_xAtCer_pions_data_random_cut_all.Write()
P_ngcer_yAtCer_pions_data_random_cut_all.Write()
P_aero_npeSum_pions_data_random_cut_all.Write()
P_aero_xAtAero_pions_data_random_cut_all.Write()
P_aero_yAtAero_pions_data_random_cut_all.Write()
MMpi_pions_data_random_cut_all.Write()
P_RFTime_Dist_pions_data_random_cut_all.Write()
CTime_ePiCoinTime_ROC2_pions_data_random_cut_all.Write()
pmiss_pions_data_random_cut_all.Write()
pmiss_x_pions_data_random_cut_all.Write()
pmiss_y_pions_data_random_cut_all.Write()
pmiss_z_pions_data_random_cut_all.Write()
emiss_pions_data_random_cut_all.Write()
W_pions_data_random_cut_all.Write()
Q2_pions_data_random_cut_all.Write()
Epsilon_pions_data_random_cut_all.Write()
t_pions_data_random_cut_all.Write()
ph_q_pions_data_random_cut_all.Write()

d_Cut_Pion_Events_RandomSub_Data.cd()
H_gtr_beta_pions_randsub_data_cut_all.Write()
H_gtr_xp_pions_randsub_data_cut_all.Write()
H_gtr_yp_pions_randsub_data_cut_all.Write()
H_gtr_dp_pions_randsub_data_cut_all.Write()
H_gtr_p_pions_randsub_data_cut_all.Write()
H_dc_x_fp_pions_randsub_data_cut_all.Write()
H_dc_y_fp_pions_randsub_data_cut_all.Write()
H_dc_xp_fp_pions_randsub_data_cut_all.Write()
H_dc_yp_fp_pions_randsub_data_cut_all.Write()
H_cal_etotnorm_pions_randsub_data_cut_all.Write()
H_cal_etottracknorm_pions_randsub_data_cut_all.Write()
H_cer_npeSum_pions_randsub_data_cut_all.Write()
H_RFTime_Dist_pions_randsub_data_cut_all.Write()
P_gtr_beta_pions_randsub_data_cut_all.Write()
P_gtr_xp_pions_randsub_data_cut_all.Write()
P_gtr_yp_pions_randsub_data_cut_all.Write()
P_gtr_dp_pions_randsub_data_cut_all.Write()
P_gtr_p_pions_randsub_data_cut_all.Write()
P_dc_x_fp_pions_randsub_data_cut_all.Write()
P_dc_y_fp_pions_randsub_data_cut_all.Write()
P_dc_xp_fp_pions_randsub_data_cut_all.Write()
P_dc_yp_fp_pions_randsub_data_cut_all.Write()
P_cal_etotnorm_pions_randsub_data_cut_all.Write()
P_cal_etottracknorm_pions_randsub_data_cut_all.Write()
P_hgcer_npeSum_pions_randsub_data_cut_all.Write()
P_hgcer_xAtCer_pions_randsub_data_cut_all.Write()
P_hgcer_yAtCer_pions_randsub_data_cut_all.Write()
P_ngcer_npeSum_pions_randsub_data_cut_all.Write()
P_ngcer_xAtCer_pions_randsub_data_cut_all.Write()
P_ngcer_yAtCer_pions_randsub_data_cut_all.Write()
P_aero_npeSum_pions_randsub_data_cut_all.Write()
P_aero_xAtAero_pions_randsub_data_cut_all.Write()
P_aero_yAtAero_pions_randsub_data_cut_all.Write()
MMpi_pions_randsub_data_cut_all.Write()
P_RFTime_Dist_pions_randsub_data_cut_all.Write()
CTime_ePiCoinTime_ROC2_pions_randsub_data_cut_all.Write()
pmiss_pions_randsub_data_cut_all.Write()
pmiss_x_pions_randsub_data_cut_all.Write()
pmiss_y_pions_randsub_data_cut_all.Write()
pmiss_z_pions_randsub_data_cut_all.Write()
emiss_pions_randsub_data_cut_all.Write()
W_pions_randsub_data_cut_all.Write()
Q2_pions_randsub_data_cut_all.Write()
Epsilon_pions_randsub_data_cut_all.Write()
t_pions_randsub_data_cut_all.Write()
ph_q_pions_randsub_data_cut_all.Write()

d_Cut_Pion_Events_Norm_DummySub_RandomSub_Data.cd()
H_gtr_beta_pions_normdummysub_data_cut_all.Write()
H_gtr_xp_pions_normdummysub_data_cut_all.Write()
H_gtr_yp_pions_normdummysub_data_cut_all.Write()
H_gtr_dp_pions_normdummysub_data_cut_all.Write()
H_gtr_p_pions_normdummysub_data_cut_all.Write()
H_dc_x_fp_pions_normdummysub_data_cut_all.Write()
H_dc_y_fp_pions_normdummysub_data_cut_all.Write()
H_dc_xp_fp_pions_normdummysub_data_cut_all.Write()
H_dc_yp_fp_pions_normdummysub_data_cut_all.Write()
H_cal_etotnorm_pions_normdummysub_data_cut_all.Write()
H_cal_etottracknorm_pions_normdummysub_data_cut_all.Write()
H_cer_npeSum_pions_normdummysub_data_cut_all.Write()
H_RFTime_Dist_pions_normdummysub_data_cut_all.Write()
P_gtr_beta_pions_normdummysub_data_cut_all.Write()
P_gtr_xp_pions_normdummysub_data_cut_all.Write()
P_gtr_yp_pions_normdummysub_data_cut_all.Write()
P_gtr_dp_pions_normdummysub_data_cut_all.Write()
P_gtr_p_pions_normdummysub_data_cut_all.Write()
P_dc_x_fp_pions_normdummysub_data_cut_all.Write()
P_dc_y_fp_pions_normdummysub_data_cut_all.Write()
P_dc_xp_fp_pions_normdummysub_data_cut_all.Write()
P_dc_yp_fp_pions_normdummysub_data_cut_all.Write()
P_cal_etotnorm_pions_normdummysub_data_cut_all.Write()
P_cal_etottracknorm_pions_normdummysub_data_cut_all.Write()
P_hgcer_npeSum_pions_normdummysub_data_cut_all.Write()
P_hgcer_xAtCer_pions_normdummysub_data_cut_all.Write()
P_hgcer_yAtCer_pions_normdummysub_data_cut_all.Write()
P_ngcer_npeSum_pions_normdummysub_data_cut_all.Write()
P_ngcer_xAtCer_pions_normdummysub_data_cut_all.Write()
P_ngcer_yAtCer_pions_normdummysub_data_cut_all.Write()
P_aero_npeSum_pions_normdummysub_data_cut_all.Write()
P_aero_xAtAero_pions_normdummysub_data_cut_all.Write()
P_aero_yAtAero_pions_normdummysub_data_cut_all.Write()
MMpi_pions_normdummysub_data_cut_all.Write()
P_RFTime_Dist_pions_normdummysub_data_cut_all.Write()
CTime_ePiCoinTime_ROC2_pions_normdummysub_data_cut_all.Write()
pmiss_pions_normdummysub_data_cut_all.Write()
pmiss_x_pions_normdummysub_data_cut_all.Write()
pmiss_y_pions_normdummysub_data_cut_all.Write()
pmiss_z_pions_normdummysub_data_cut_all.Write()
emiss_pions_normdummysub_data_cut_all.Write()
W_pions_normdummysub_data_cut_all.Write()
Q2_pions_normdummysub_data_cut_all.Write()
Epsilon_pions_normdummysub_data_cut_all.Write()
t_pions_normdummysub_data_cut_all.Write()

d_Cut_Pion_Events_Norm_SIMC_Data.cd()
H_hsdelta_pions_simc_cut_all.Write()
H_hsxptar_pions_simc_cut_all.Write()
H_hsyptar_pions_simc_cut_all.Write()
H_hsytar_pions_simc_cut_all.Write()
H_hsxfp_pions_simc_cut_all.Write()
H_hsyfp_pions_simc_cut_all.Write()
H_hsxpfp_pions_simc_cut_all.Write()
H_hsypfp_pions_simc_cut_all.Write()
P_ssdelta_pions_simc_cut_all.Write()
P_ssxptar_pions_simc_cut_all.Write()
P_ssyptar_pions_simc_cut_all.Write()
P_ssytar_pions_simc_cut_all.Write()
P_ssxfp_pions_simc_cut_all.Write()
P_ssyfp_pions_simc_cut_all.Write()
P_ssxpfp_pions_simc_cut_all.Write()
P_ssypfp_pions_simc_cut_all.Write()
q_pions_simc_cut_all.Write()
nu_pions_simc_cut_all.Write()
Q2_pions_simc_cut_all.Write()
epsilon_pions_simc_cut_all.Write()
thetapq_pions_simc_cut_all.Write()
phipq_pions_simc_cut_all.Write()
pmiss_pions_simc_cut_all.Write()
pmiss_x_pions_simc_cut_all.Write()
pmiss_y_pions_simc_cut_all.Write()
pmiss_z_pions_simc_cut_all.Write()
emiss_pions_simc_cut_all.Write()
W_pions_simc_cut_all.Write()
t_pions_simc_cut_all.Write()
MMpi_pions_simc_cut_all.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_DATA.Close() 
infile_DUMMY.Close()
infile_SIMC.Close()
outHistFile.Close()

print ("Processing Complete")