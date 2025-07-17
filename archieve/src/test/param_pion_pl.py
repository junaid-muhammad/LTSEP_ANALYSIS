#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-08 08:38:16 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import math

def iterWeight(arg_str):

    # Define constants
    pi = 3.14159
    mtar_gev = 0.93827231
    mpipl=0.139570
    mkpl=0.493677

    # Split and convert the input string into a list of floats
    args = list(map(float, arg_str.split()))

    # Extract individual values from the list
    q2_set, q2_sim, w_sim, t_sim, eps_sim, thetacm_sim, phicm_sim, sigcm_sim, wt_sim, *params = args
    
    q2_gev = q2_sim # Already GeV
    t_gev = t_sim  # Already GeV, issue here!!! t_sim makes no sense
    s = w_sim**2
    s_gev = s # Already GeV

    ##############
    # HARD CODED #
    ##############
    if q2_set == 2.1:
        q2_set = 2.115
    ##############
    ##############
    ##############

    # Calculate tav, ftav, ft
    #tav = (0.0735 + 0.028 * math.log(q2_set)) * q2_set
    # RLT (10/8/2023): Testing new tav parameterization
    tav=(0.1112 + 0.0066*math.log(q2_set))*q2_set
    ftav = (abs(t_gev) - tav) / tav

    # RLT (7/11/2024): Moved below for Q2dep func form
    #ft = abs(t_gev) / (abs(t_gev) + mkpl**2)**2
    
    # Calculate sigl, sigt, siglt, sigtt, sig219, sig
    # RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for the
    #                 xfit_in_t.py script to work. LT/TT are zeros
    #p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = params
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16 = params

    # RLT (4/5/2024): Added try statments to all sep xsect to suppress any blowing up for some events.
    #                 Setting bad events to zero
    
    try:
        #sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev) - 0.2))
        # RLT (10/12/2023): Removed 0.2 to keep things as simple as possible for initial start parameterization
        # RLT (2/19/2024): Adding a 0.2 term to t dependence to bring down the extreme slope at high t
        # RLT (3/09/2024): Removing +0.2 term for better parameterization of Q2=3.0, W=2.32
        #
        #sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev)))
        #sigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (abs(t_gev)+0.2))
        # RLT (4/23/2024): Marco's thesis functional forms
        #sigl = p1 * math.exp(-p2*abs(t_gev)) * (1.0 / (1.0 + p3*q2_gev))
        # RLT (6/04/2024): Testing simplier exp form for L+T
        ##
        ##sigl = (p1 + p2 * math.log(q2_gev)) * math.exp(p3 * (abs(t_gev)))
        #sigl = (p1 * ((abs(t_gev)/q2_gev)-1)) * math.exp(p2 * (abs(t_gev)))
        ##                
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        ft = abs(t_gev) / (abs(t_gev) + mkpl**2)**2 # pole term
        Qdep_L=q2_gev/(1.0+(1.77*q2_gev)+0.12*(q2_gev**2))
        ##sigl=(p1*Qdep_L*ft)*math.exp(-p2*(abs(t_gev)))
        sigl=(p1*Qdep_L*ft)*math.exp(-p2*(abs(t_gev)))

    except OverflowError:
        sigl = -1000.0
        #print("WARNING: Overflowerror on sigL, setting to zero for this event...")

    try:
        # RLT (2/15/2024): Removing t dependence from sigT because it seems
        #                  to be driving poor sep xsects results
        # RLT (2/20/2024): Added 1/Q^4 term to dampen sigT
        # RLT (2/21/2024): Using global analysis sig T model and params (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.018202)
        #
        #sigt = p5 + p6 * math.log(q2_gev) + (p7 + p8 * math.log(q2_gev)) * ftav
        #sigt = p5 + p6 * math.log(q2_gev)
        #sigt = p5 * math.log(q2_gev) + p6 / (q2_gev**2)
        #sigt = p5 / (1 + p6*q2_gev)
        # RLT (4/20/2024): Adding in t-dependence
        #sigt = (p5 / (1 + p6*q2_gev)) * ftav
        #sigt = (p5 / (1 + p6*q2_gev)) * abs(t_gev)
        # RLT (4/23/2024): Exponential t-dependence
        #sigt = (p5 / (1 + p6*q2_gev)) * math.exp(p7*abs(t_gev))
        # RLT (4/23/2024): Marco's thesis functional forms
        #sigt = p5 * math.exp(-p6*abs(t_gev)) * (1.0 / (1.0 + p7*q2_gev))
        # RLT (6/04/2024): Testing simplier exp form for L+T
        ##
        ##sigt = (p5 * ((abs(t_gev)/q2_gev)-1)) * math.exp(p6 * (abs(t_gev)))
        #sigt = (p5 + p6 * math.log(q2_gev)) * math.exp(p7 * (abs(t_gev)))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        #sigt=(p5/q2_gev)*math.exp(-p6*(q2_gev**2))
        Qdep_T=(math.exp(-q2_gev**2))/q2_gev
        #sigt=p5*(p6+math.exp(-p7*(abs(t_gev))))*(Qdep_T**p8)
        sigt=(p5*math.exp(-p6*(abs(t_gev)))+p7*(abs(t_gev)))*(Qdep_T**p8)

    except OverflowError:        
        sigt = -1000.0
        #print("WARNING: Overflowerror on sigT, setting to zero for this event...")

    try:
        ##
        siglt = (p9 * math.exp(p10 * abs(t_gev)) + p11 / abs(t_gev)) * math.sin(thetacm_sim)
        #siglt = (p9 + p11 / abs(t_gev)) * math.sin(thetacm_sim)
        # RLT (4/23/2024): Marco's thesis functional forms
        #siglt = p9 * math.exp(-p10*abs(t_gev)) * (1.0 / (1.0 + (q2_gev**2)*p11))
        ##        
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        ###siglt=(p9/(1+q2_gev))*math.sin(thetacm_sim)*math.exp(-p10*(abs(t_gev)))
        #siglt=(p9 * math.exp(p10 * abs(t_gev)) + p11 / abs(t_gev)) * math.sin(thetacm_sim)

    except OverflowError:
        siglt = -1000.0
        #print("WARNING: Overflowerror on sigLT, setting to zero for this event...\n\n")

    try:
        # RLT (1/2/2024): Need to have 16 parameters (4 for L/T/LT/TT) for the
        #                 xfit_in_t.py script to work. LT/TT are zeros
        #                 Therefore param 12 was also changed to 13
        ##
        sigtt = (p13 * q2_gev * math.exp(-q2_gev)) * ft * math.sin(thetacm_sim)**2
        # RLT (4/23/2024): Marco's thesis functional forms
        #sigtt = p13 * math.exp(-p14*abs(t_gev)) * (1.0 / (1.0 + (q2_gev**2)*p15))
        ##
        # RLT (7/11/2024): Redefined functional forms of L, T, LT, TT
        #                  that incorporates Q2-dep based of pi FF
        ##sigtt=(-p13/(1+q2_gev))*(math.sin(thetacm_sim)**2)*math.exp(-p14*abs(t_gev))
        ###sigtt=(p13/(1+q2_gev))*(math.sin(thetacm_sim)**2)*ft*math.exp(-p14*(q2_gev))
        
    except OverflowError:
        sigtt = -1000.0
        #print("WARNING: Overflowerror on sigTT, setting to zero for this event...\n\n")

    # RLT (9/25/2023): There are two tav parameterizations in here.
    #                  I am only using the one above, for now.    
    # tav = (-0.178 + 0.315 * math.log(q2_gev)) * q2_gev

    wfactor = 1.0 / (s_gev - mtar_gev**2)**2
    sigl = sigl*wfactor
    sigt = sigt*wfactor
    sigtt = sigtt*wfactor
    siglt = siglt*wfactor

    sig = (sigt + eps_sim * sigl + eps_sim * math.cos(2. * phicm_sim) * sigtt +
             math.sqrt(2.0 * eps_sim * (1. + eps_sim)) * math.cos(phicm_sim) * siglt)
    
    sig = sig / 2.0 / pi / 1e6  # dsig/dtdphicm in microbarns/MeV**2/rad
    #sig = sig / 2.0 / pi  # dsig/dtdphicm in microbarns/GeV**2/rad

    wtn = wt_sim * sig / sigcm_sim

    #print("sig",sig)
    #print("sigcm",sigcm_sim)
    #print("wtn",wtn)
    #print("wt_sim",wt_sim)
    
    return [float(wtn),float(sig)]
