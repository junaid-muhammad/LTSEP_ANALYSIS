#! /usr/bin/python

#
# Description: Adapted from fortran code wt28_3.f
# ================================================================
# Time-stamp: "2025-02-04 16:18:51 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import root_numpy as rnp
import ROOT
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
from array import array
import sys, math, os, subprocess

# Import the dynamic script
import importlib.util

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import open_root_file, run_fortran

##################################################################################################################################################
# Importing param model for weight iteration

sys.path.append("models")
from param_active import set_val, iterWeight

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def iter_weight(param_file, simc_root, inpDict, phi_setting):
    
    formatted_date  = inpDict["formatted_date"]
    iter_num = inpDict["iter_num"]
    ParticleType  = inpDict["ParticleType"]
    q2_set = inpDict["Q2"]
    Q2 = q2_set.replace("p",".")
    w_set = inpDict["W"]
    W = w_set.replace("p",".")
    POL = inpDict["POL"]
    
    if int(POL) == 1:
        pol_str = "pl"
    elif int(POL) == -1:
        pol_str = "mn"
    else:
        print("ERROR: Invalid polarity...must be +1 or -1")
        sys.exit(2)
        
    param_arr = []
    with open(param_file, 'r') as f:
        for i, line in enumerate(f):
            columns = line.split()
            param_arr.append(str(columns[0]))

    if not os.path.isfile(simc_root):
        print("\n\nERROR: No simc file found called {}\n\n".format(simc_root))        
        
    InFile_SIMC = open_root_file(simc_root, "READ")
    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    if iter_num > 1:
        # Create a new ROOT file for writing
        new_InFile_SIMC = open_root_file(simc_root.replace("iter_{}.root".format(iter_num-1),"iter_{}.root".format(iter_num)), "UPDATE")
        new_TBRANCH_SIMC = ROOT.TTree("h10", "Iteration {}".format(iter_num))
    else:
        # Create a new ROOT file for writing
        new_InFile_SIMC = open_root_file(simc_root.replace(".root","_iter_{}.root".format(iter_num)), "UPDATE")
        new_TBRANCH_SIMC = ROOT.TTree("h10", "Iteration {}".format(iter_num))        

    # Grab branches from previous iteration
    hsdelta_array = array( 'f', [0])
    hsyptar_array = array( 'f', [0])
    hsxptar_array = array( 'f', [0])
    hsytar_array = array( 'f', [0])
    hsxfp_array = array( 'f', [0])
    hsxpfp_array = array( 'f', [0])
    hsyfp_array = array( 'f', [0])
    hsypfp_array = array( 'f', [0])
    hsdeltai_array = array( 'f', [0])
    hsyptari_array = array( 'f', [0])
    hsxptari_array = array( 'f', [0])
    hsytari_array = array( 'f', [0])
    ssdelta_array = array( 'f', [0])
    ssyptar_array = array( 'f', [0])
    ssxptar_array = array( 'f', [0])
    ssytar_array = array( 'f', [0])
    ssxfp_array = array( 'f', [0])
    ssxpfp_array = array( 'f', [0])
    ssyfp_array = array( 'f', [0])
    ssypfp_array = array( 'f', [0])
    ssdeltai_array = array( 'f', [0])
    ssyptari_array = array( 'f', [0])
    ssxptari_array = array( 'f', [0])
    ssytari_array = array( 'f', [0])
    q_array = array( 'f', [0])
    nu_array = array( 'f', [0])
    Q2_array = array( 'f', [0])
    W_array = array( 'f', [0])
    epsilon_array = array( 'f', [0])
    epscm_array = array( 'f', [0])
    Em_array = array( 'f', [0])
    Pm_array = array( 'f', [0])
    thetapq_array = array( 'f', [0])
    thetacm_array = array( 'f', [0])
    phipq_array = array( 'f', [0])
    missmass_array = array( 'f', [0])
    missmass_shift_array = array( 'f', [0])    
    mmnuc_array = array( 'f', [0])
    phad_array = array( 'f', [0])
    t_array = array( 'f', [0])
    pmpar_array = array( 'f', [0])
    pmper_array = array( 'f', [0])
    pmoop_array = array( 'f', [0])
    fry_array = array( 'f', [0])
    radphot_array = array( 'f', [0])
    pfermi_array = array( 'f', [0])
    siglab_array = array( 'f', [0])
    decdist_array = array( 'f', [0])
    Mhadron_array = array( 'f', [0])
    pdotqhat_array = array( 'f', [0])
    Q2i_array = array( 'f', [0])
    Wi_array = array( 'f', [0])
    ti_array = array( 'f', [0])
    phipqi_array = array( 'f', [0])
    saghai_array = array( 'f', [0])
    factor_array = array( 'f', [0])
    paero_z_det_array = array( 'f', [0])
    paero_x_det_array = array( 'f', [0])
    paero_y_det_array = array( 'f', [0])
    phgcer_z_det_array = array( 'f', [0])
    phgcer_x_det_array = array( 'f', [0])
    phgcer_y_det_array = array( 'f', [0])
    pend_z_det_array = array( 'f', [0])
    pend_x_det_array = array( 'f', [0])
    pend_y_det_array = array( 'f', [0])
    TBRANCH_SIMC.SetBranchAddress("hsdelta", hsdelta_array);
    TBRANCH_SIMC.SetBranchAddress("hsyptar", hsyptar_array);
    TBRANCH_SIMC.SetBranchAddress("hsxptar", hsxptar_array);
    TBRANCH_SIMC.SetBranchAddress("hsytar", hsytar_array);
    TBRANCH_SIMC.SetBranchAddress("hsxfp", hsxfp_array);
    TBRANCH_SIMC.SetBranchAddress("hsxpfp", hsxpfp_array);
    TBRANCH_SIMC.SetBranchAddress("hsyfp", hsyfp_array);
    TBRANCH_SIMC.SetBranchAddress("hsypfp", hsypfp_array);
    TBRANCH_SIMC.SetBranchAddress("hsdeltai", hsdeltai_array);
    TBRANCH_SIMC.SetBranchAddress("hsyptari", hsyptari_array);
    TBRANCH_SIMC.SetBranchAddress("hsxptari", hsxptari_array);
    TBRANCH_SIMC.SetBranchAddress("hsytari", hsytari_array);
    TBRANCH_SIMC.SetBranchAddress("ssdelta", ssdelta_array);
    TBRANCH_SIMC.SetBranchAddress("ssyptar", ssyptar_array);
    TBRANCH_SIMC.SetBranchAddress("ssxptar", ssxptar_array);
    TBRANCH_SIMC.SetBranchAddress("ssytar", ssytar_array);
    TBRANCH_SIMC.SetBranchAddress("ssxfp", ssxfp_array);
    TBRANCH_SIMC.SetBranchAddress("ssxpfp", ssxpfp_array);
    TBRANCH_SIMC.SetBranchAddress("ssyfp", ssyfp_array);
    TBRANCH_SIMC.SetBranchAddress("ssypfp", ssypfp_array);
    TBRANCH_SIMC.SetBranchAddress("ssdeltai", ssdeltai_array);
    TBRANCH_SIMC.SetBranchAddress("ssyptari", ssyptari_array);
    TBRANCH_SIMC.SetBranchAddress("ssxptari", ssxptari_array);
    TBRANCH_SIMC.SetBranchAddress("ssytari", ssytari_array);
    TBRANCH_SIMC.SetBranchAddress("q", q_array);
    TBRANCH_SIMC.SetBranchAddress("nu", nu_array);
    TBRANCH_SIMC.SetBranchAddress("Q2", Q2_array);
    TBRANCH_SIMC.SetBranchAddress("W", W_array);
    TBRANCH_SIMC.SetBranchAddress("epsilon", epsilon_array);
    TBRANCH_SIMC.SetBranchAddress("epscm", epscm_array);  
    TBRANCH_SIMC.SetBranchAddress("Em", Em_array);
    TBRANCH_SIMC.SetBranchAddress("Pm", Pm_array);
    TBRANCH_SIMC.SetBranchAddress("thetapq", thetapq_array);
    TBRANCH_SIMC.SetBranchAddress("thetacm", thetacm_array);
    TBRANCH_SIMC.SetBranchAddress("phipq", phipq_array);
    TBRANCH_SIMC.SetBranchAddress("missmass", missmass_array);
    TBRANCH_SIMC.SetBranchAddress("missmass_shift", missmass_shift_array);    
    TBRANCH_SIMC.SetBranchAddress("mmnuc", mmnuc_array);
    TBRANCH_SIMC.SetBranchAddress("mmnuc", mmnuc_array);    
    TBRANCH_SIMC.SetBranchAddress("phad", phad_array);
    TBRANCH_SIMC.SetBranchAddress("t", t_array);
    TBRANCH_SIMC.SetBranchAddress("pmpar", pmpar_array);
    TBRANCH_SIMC.SetBranchAddress("pmper", pmper_array);
    TBRANCH_SIMC.SetBranchAddress("pmoop", pmoop_array);
    TBRANCH_SIMC.SetBranchAddress("fry", fry_array);
    TBRANCH_SIMC.SetBranchAddress("radphot", radphot_array);
    TBRANCH_SIMC.SetBranchAddress("pfermi", pfermi_array);
    TBRANCH_SIMC.SetBranchAddress("siglab", siglab_array);
    TBRANCH_SIMC.SetBranchAddress("decdist", decdist_array);
    TBRANCH_SIMC.SetBranchAddress("Mhadron", Mhadron_array);
    TBRANCH_SIMC.SetBranchAddress("pdotqhat", pdotqhat_array);
    TBRANCH_SIMC.SetBranchAddress("Q2i", Q2i_array);
    TBRANCH_SIMC.SetBranchAddress("Wi", Wi_array);
    TBRANCH_SIMC.SetBranchAddress("ti", ti_array);
    TBRANCH_SIMC.SetBranchAddress("phipqi", phipqi_array);
    TBRANCH_SIMC.SetBranchAddress("saghai", saghai_array);
    TBRANCH_SIMC.SetBranchAddress("factor", factor_array);
    TBRANCH_SIMC.SetBranchAddress("paero_z_det", paero_z_det_array)
    TBRANCH_SIMC.SetBranchAddress("paero_x_det", paero_x_det_array)
    TBRANCH_SIMC.SetBranchAddress("paero_y_det", paero_y_det_array)
    TBRANCH_SIMC.SetBranchAddress("phgcer_z_det", phgcer_z_det_array)
    TBRANCH_SIMC.SetBranchAddress("phgcer_x_det", phgcer_x_det_array)
    TBRANCH_SIMC.SetBranchAddress("phgcer_y_det", phgcer_y_det_array)
    TBRANCH_SIMC.SetBranchAddress("pend_z_det", pend_z_det_array)
    TBRANCH_SIMC.SetBranchAddress("pend_x_det", pend_x_det_array)
    TBRANCH_SIMC.SetBranchAddress("pend_y_det", pend_y_det_array)

    # Define new iteration
    Weight_array = array( 'f', [0])        
    sigcm_array = array( 'f', [0])
    iter_weight_array = array( 'f', [0])
    iter_sig_array = array( 'f', [0])
    TBRANCH_SIMC.SetBranchAddress("Weight", Weight_array)
    TBRANCH_SIMC.SetBranchAddress("sigcm", sigcm_array)
    if iter_num > 1:
        TBRANCH_SIMC.SetBranchAddress("iter_weight", iter_weight_array)
        TBRANCH_SIMC.SetBranchAddress("iter_sig", iter_sig_array)

    # Create branches for current iteration
    new_TBRANCH_SIMC.Branch("hsdelta", hsdelta_array, "hsdelta/F");
    new_TBRANCH_SIMC.Branch("hsyptar", hsyptar_array, "hsyptar/F");
    new_TBRANCH_SIMC.Branch("hsxptar", hsxptar_array, "hsxptar/F");
    new_TBRANCH_SIMC.Branch("hsytar", hsytar_array, "hsytar/F");
    new_TBRANCH_SIMC.Branch("hsxfp", hsxfp_array, "hsxfp/F");
    new_TBRANCH_SIMC.Branch("hsxpfp", hsxpfp_array, "hsxpfp/F");
    new_TBRANCH_SIMC.Branch("hsyfp", hsyfp_array, "hsyfp/F");
    new_TBRANCH_SIMC.Branch("hsypfp", hsypfp_array, "hsypfp/F");
    new_TBRANCH_SIMC.Branch("hsdeltai", hsdeltai_array, "hsdeltai/F");
    new_TBRANCH_SIMC.Branch("hsyptari", hsyptari_array, "hsyptari/F");
    new_TBRANCH_SIMC.Branch("hsxptari", hsxptari_array, "hsxptari/F");
    new_TBRANCH_SIMC.Branch("hsytari", hsytari_array, "hsytari/F");
    new_TBRANCH_SIMC.Branch("ssdelta", ssdelta_array, "ssdelta/F");
    new_TBRANCH_SIMC.Branch("ssyptar", ssyptar_array, "ssyptar/F");
    new_TBRANCH_SIMC.Branch("ssxptar", ssxptar_array, "ssxptar/F");
    new_TBRANCH_SIMC.Branch("ssytar", ssytar_array, "ssytar/F");
    new_TBRANCH_SIMC.Branch("ssxfp", ssxfp_array, "ssxfp/F");
    new_TBRANCH_SIMC.Branch("ssxpfp", ssxpfp_array, "ssxpfp/F");
    new_TBRANCH_SIMC.Branch("ssyfp", ssyfp_array, "ssyfp/F");
    new_TBRANCH_SIMC.Branch("ssypfp", ssypfp_array, "ssypfp/F");
    new_TBRANCH_SIMC.Branch("ssdeltai", ssdeltai_array, "ssdeltai/F");
    new_TBRANCH_SIMC.Branch("ssyptari", ssyptari_array, "ssyptari/F");
    new_TBRANCH_SIMC.Branch("ssxptari", ssxptari_array, "ssxptari/F");
    new_TBRANCH_SIMC.Branch("ssytari", ssytari_array, "ssytari/F");
    new_TBRANCH_SIMC.Branch("q", q_array, "q/F");
    new_TBRANCH_SIMC.Branch("nu", nu_array, "nu/F");
    new_TBRANCH_SIMC.Branch("Q2", Q2_array, "Q2/F");
    new_TBRANCH_SIMC.Branch("W", W_array, "W/F");
    new_TBRANCH_SIMC.Branch("epsilon", epsilon_array, "epsilon/F");
    new_TBRANCH_SIMC.Branch("epscm", epscm_array, "epscm/F");  
    new_TBRANCH_SIMC.Branch("Em", Em_array, "Em/F");
    new_TBRANCH_SIMC.Branch("Pm", Pm_array, "Pm/F");
    new_TBRANCH_SIMC.Branch("thetapq", thetapq_array, "thetapq/F");
    new_TBRANCH_SIMC.Branch("thetacm", thetacm_array, "thetacm/F");
    new_TBRANCH_SIMC.Branch("phipq", phipq_array, "phipq/F");
    new_TBRANCH_SIMC.Branch("missmass", missmass_array, "missmass/F");
    new_TBRANCH_SIMC.Branch("missmass_shift", missmass_shift_array, "missmass_shift/F");    
    new_TBRANCH_SIMC.Branch("mmnuc", mmnuc_array, "mmnuc/F");
    new_TBRANCH_SIMC.Branch("phad", phad_array, "phad/F");
    new_TBRANCH_SIMC.Branch("t", t_array, "t/F");
    new_TBRANCH_SIMC.Branch("pmpar", pmpar_array, "pmpar/F");
    new_TBRANCH_SIMC.Branch("pmper", pmper_array, "pmper/F");
    new_TBRANCH_SIMC.Branch("pmoop", pmoop_array, "pmoop/F");
    new_TBRANCH_SIMC.Branch("fry", fry_array, "fry/F");
    new_TBRANCH_SIMC.Branch("radphot", radphot_array, "radphot/F");
    new_TBRANCH_SIMC.Branch("pfermi", pfermi_array, "pfermi/F");
    new_TBRANCH_SIMC.Branch("siglab", siglab_array, "siglab/F");
    new_TBRANCH_SIMC.Branch("decdist", decdist_array, "decdist/F");
    new_TBRANCH_SIMC.Branch("Mhadron", Mhadron_array, "Mhadron/F");
    new_TBRANCH_SIMC.Branch("pdotqhat", pdotqhat_array, "pdotqhat/F");
    new_TBRANCH_SIMC.Branch("Q2i", Q2i_array, "Q2i/F");
    new_TBRANCH_SIMC.Branch("Wi", Wi_array, "Wi/F");
    new_TBRANCH_SIMC.Branch("ti", ti_array, "ti/F");
    new_TBRANCH_SIMC.Branch("phipqi", phipqi_array, "phipqi/F");
    new_TBRANCH_SIMC.Branch("saghai", saghai_array, "saghai/F");
    new_TBRANCH_SIMC.Branch("factor", factor_array, "factor/F");        
    new_TBRANCH_SIMC.Branch("Weight", Weight_array, "Weight/F")
    new_TBRANCH_SIMC.Branch("sigcm", sigcm_array, "sigcm/F")
    new_TBRANCH_SIMC.Branch("iter_weight", iter_weight_array, "iter_weight/F")
    new_TBRANCH_SIMC.Branch("iter_sig", iter_sig_array, "iter_sig/F")
    new_TBRANCH_SIMC.Branch("paero_z_det", paero_z_det_array, "paero_z_det/F")
    new_TBRANCH_SIMC.Branch("paero_x_det", paero_x_det_array, "paero_x_det/F")
    new_TBRANCH_SIMC.Branch("paero_y_det", paero_y_det_array, "paero_y_det/F")
    new_TBRANCH_SIMC.Branch("phgcer_z_det", phgcer_z_det_array, "phgcer_z_det/F")
    new_TBRANCH_SIMC.Branch("phgcer_x_det", phgcer_x_det_array, "phgcer_x_det/F")
    new_TBRANCH_SIMC.Branch("phgcer_y_det", phgcer_y_det_array, "phgcer_y_det/F")
    new_TBRANCH_SIMC.Branch("pend_z_det", pend_z_det_array, "pend_z_det/F")
    new_TBRANCH_SIMC.Branch("pend_x_det", pend_x_det_array, "pend_x_det/F")
    new_TBRANCH_SIMC.Branch("pend_y_det", pend_y_det_array, "pend_y_det/F")

    # Set pol_str, q2_set, w_set for param model script
    set_val(pol_str, q2_set, w_set)
    
    ################################################################################################################################################
    # Run over simc root branch to determine new weight

    # Keep track of bad events
    bad_events = []
    total_events = TBRANCH_SIMC.GetEntries()
    
    print("\nRecalculating weight for %s simc..." % phi_setting)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, total_events,bar_length=25)

      TBRANCH_SIMC.GetEntry(i)

      if iter_num > 1:

          #evt.Weight = evt.iter_weight
          #evt.sigcm = evt.iter_sig

          # Note: ti is used instead of t, ti = main%t which matches its calculation in simc
          #       while t is calculated in recon_hcana (but should be invariant?? Not sure the issue)
          #       This goes for Q2i, Wi, and phiqpi as well
          #inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #print("-"*25,"\n",i,"\n",inp_param)
          #inp_param = '{} {} {} {} {} {} {} {} {} '\
              #.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.iter_sig, evt.iter_weight)+' '.join(param_arr)
          #.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          ###################
          # Before 8/28/2024
          #inp_param = '{} {} {} {} {} {} {} {} {} {} '\
                      #.format(Q2, W, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetapq, evt.phipq, evt.iter_sig, evt.iter_weight)+' '.join(param_arr)
          # After 8/28/2024
          inp_param = '{} {} {} {} {} {} {} {} {} {} '\
                      .format(Q2, W, evt.Q2i, evt.Wi, evt.ti, evt.epsilon, evt.thetapq, evt.phipqi, evt.iter_sig, evt.iter_weight)+' '.join(param_arr)

          iter_lst = iterWeight(inp_param)

          # Check for bad events
          if iter_lst[0] == 0.0:
              bad_events.append(i)
              continue
          
          Weight_array[0] = evt.iter_weight
          sigcm_array[0] = evt.iter_sig
          
          iter_weight_array[0] = iter_lst[0]
          iter_sig_array[0] = iter_lst[1]

          new_TBRANCH_SIMC.Fill()
          
      else:

          # Note: ti is used instead of t, ti = main%t which matches its calculation in simc
          #       while t is calculated in recon_hcana (but should be invariant?? Not sure the issue)
          #       This goes for Q2i, Wi, and phiqpi as well
          #inp_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          #print("-"*25,"\n",i,"\n",inp_param)
          #inp_param = '{} {} {} {} {} {} {} {} {} '\
              #.format(Q2, evt.Q2i, evt.Wi, evt.ti, evt.epscm, evt.thetacm, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)
          ###################
          # Before 8/28/2024
          #inp_param = '{} {} {} {} {} {} {} {} {} {} '\
                      #.format(Q2, W, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetapq, evt.phipq, evt.sigcm, evt.Weight)+' '.join(param_arr)
          # After 8/28/2024
          inp_param = '{} {} {} {} {} {} {} {} {} {} '\
                      .format(Q2, W, evt.Q2i, evt.Wi, evt.ti, evt.epsilon, evt.thetapq, evt.phipqi, evt.sigcm, evt.Weight)+' '.join(param_arr)

          iter_lst = iterWeight(inp_param)

          # Check for bad events
          if iter_lst[0] == 0.0:
              bad_events.append(i)
              continue
          
          Weight_array[0] = evt.Weight
          sigcm_array[0] = evt.sigcm
          
          iter_weight_array[0] = iter_lst[0]
          iter_sig_array[0] = iter_lst[1]

          new_TBRANCH_SIMC.Fill()

    new_TBRANCH_SIMC.Write("h10",ROOT.TObject.kOverwrite)
    new_InFile_SIMC.Close()
    InFile_SIMC.Close()

    print(f"\n\nThere were {len(bad_events)}/{total_events} bad events skipped...")
