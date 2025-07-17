#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-04-14 10:43:03 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math, sys
import logging

logging.basicConfig(level=logging.DEBUG)

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations, prepare_equations

###############################################################################################################################################
# Need to grab polarity Q2 and W string values from xfit script

# Check output equations
DEBUG=False

# First, define empty strings
pol_str = ""
Q2 = ""
W = ""
equations = ""

# Then, set global variables which is called with arguments defined in xfit script
def set_val(inp_pol_str, inp_Q2, inp_W):
    global pol_str, Q2, W, equations
    pol_str = inp_pol_str
    Q2 = inp_Q2
    W = inp_W
    # Load equations from model input file of given setting
    equations = load_equations(f"Q{Q2}W{W}.model")
    if DEBUG:    
        logging.debug(f"Loaded equations: {equations}")
        
###############################################################################################################################################

def iterWeight(arg_str):
    
    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, w_set, qq, ww, tt, eps, theta_cm, phi_cm, sig_prev_iter, weight_prev_iter, *params = args
    par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12, par13, par14, par15, par16 = params
    
    # Grab functional forms from model input file
    fun_Sig_L_optimized = prepare_equations(equations, 'sig_L')
    fun_Sig_T_optimized = prepare_equations(equations, 'sig_T')
    fun_Sig_LT_optimized = prepare_equations(equations, 'sig_LT')
    fun_Sig_TT_optimized = prepare_equations(equations, 'sig_TT')
    fun_wfactor_optimized = prepare_equations(equations, 'wfactor')

    # Calculate SigL, SigT, SigLT, SigTT
    sig_L = fun_Sig_L_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par1, par2, par3, par4)
    sig_T = fun_Sig_T_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par5, par6, par7, par8)
    sig_LT = fun_Sig_LT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par9, par10, par11, par12)
    sig_TT = fun_Sig_TT_optimized(q2_set, w_set, qq, ww, tt, theta_cm, par13, par14, par15, par16)
    
    # Calculate W-factor
    wfactor = fun_wfactor_optimized(q2_set, w_set, qq, ww, tt)    

    # Convert degrees to radians
    theta_cm = theta_cm * math.pi/180
    phi_cm = phi_cm * math.pi/180

    sig = (sig_T + eps * sig_L + eps * math.cos(2. * phi_cm) * sig_TT +
             math.sqrt(2.0 * eps * (1. + eps)) * math.cos(phi_cm) * sig_LT)

    sig = sig * wfactor

    sig = sig / 2.0 / math.pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    
    try:
        wtn = weight_prev_iter * (sig / sig_prev_iter)
    except ZeroDivisionError:
        wtn = 0.0

    #print("sig",sig)
    #print("sigcm",sig_prev_iter)
    #print("wtn",wtn)
    #print("weight_prev_iter",weight_prev_iter)
    
    if isinstance(wtn, complex) or wtn <= 0.0:
        return [0.0, 0.0]
    else:
        return [float(wtn),float(sig)]
