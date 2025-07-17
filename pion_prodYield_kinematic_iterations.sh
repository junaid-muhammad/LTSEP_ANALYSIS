#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2025-03-11 12:12:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################

while getopts 'hsra' flag; do
    case "${flag}" in
        h)
        echo "-------------------------------------------------------------------"
        echo "./pion_prodYield_kinematic_analysis.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the physics analysis..."         
        echo "    -h, help" 
        echo "    -s, run simc yield script to determine simc yields"
        echo "        simc -> Flag=-s Iteration=arg1 RunList=arg2 MaxEvents=arg3-(Optional)"
        echo "    -r, run plotting script to compare Data/SIMC and determine ratios"
        echo "        comp -> Flag=-p Iteration=arg1 RunList=arg2 MaxEvents=arg3-(Optional)"
        echo "    -a, run script to calculate average kinematics and yields for LTSep analysis "
        echo "        comp -> Flag=-a Iteration=arg1 RunList=arg2 MaxEvents=arg3-(Optional)"
        exit 0
        ;;
        s) s_flag='true' ;;
        r) r_flag='true' ;;
        a) a_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

echo "Starting analysis of Pion events"
echo "I take as arguments the -flag, physics setting, runlist and number of events!"

# Input params - beam energy, run list and max number of events
Iteration=$2
if [[ -z "$2" || "$2" -le 00 ]]; then
    echo "Please provide Iteration number. It should be greater than 00!"
    exit 1
fi

RunList=$3
if [[ -z "$3" ]]; then
    echo "I need a runlist"
    echo "Please provide a run list as input"
fi

MAXEVENTS=$4
if [[ -z "$4" ]]; then
    echo "Only Run list entered...I'll assume -1 (all) events!" 
    MAXEVENTS=-1 
fi

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
#    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
    PATHFILE_INFO=`python3 /u/group/c-pionlt/USERS/${USER}/replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"qcd"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd $REPLAYPATH
# Input Arguments for Iteration Scripts
inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/PionLT_2021_2022/${RunList}"
PHY_SETTING=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}')
SIMC_SETTING=$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3 "_" $4 $5}')
DATA_Suffix=ProdCoin_Analysed_Data
DUMMY_Suffix=ProdCoin_Analysed_Dummy_Data
SIMC_Suffix="Prod_Coin_${SIMC_SETTING}"
DATA_RUN_LIST=${PHY_SETTING}
DUMMY_RUN_LIST=${PHY_SETTING}_dummy
CSV_FILE=PionLT_coin_production_Prod_efficiency_data_2025_03_08

# Input Arguments for avergae kinematics and yields calculation Script
PHY_SETTING_C=$(echo "${RunList}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
SIMC_Suffix_C=$(echo "${RunList}" | awk -F'_' '{print $1 $2 $3}')
RUN_LIST_C=${PHY_SETTING_C}_Runlist
DATA_YIELD_CSV=Physics_Data_Yield
DATA_AVG_KIN_CSV=Physics_Avg_Data_Kinematics
SIMC_YIELD_CSV=Physics_SIMC_Yield
SIMC_AVG_KIN_CSV=Physics_Avg_SIMC_Kinematics

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Moving SIMC ROOT files to required standard directory
SIMC_Suffix_File=Prod_Coin_${SIMC_Suffix_C}

if [ ! -d "${VOLATILEPATH}/worksim/${PHY_SETTING_C}_iter${Iteration}" ]; then
    mkdir -p "${VOLATILEPATH}/worksim/${PHY_SETTING_C}_iter${Iteration}"
fi
# Only move if destination does NOT already have files
if ! ls ${VOLATILEPATH}/worksim/${PHY_SETTING_C}_iter${Iteration}/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
    if ls ${VOLATILEPATH}/worksim/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
        mv ${VOLATILEPATH}/worksim/${SIMC_Suffix_File}* "${VOLATILEPATH}/worksim/${PHY_SETTING_C}_iter${Iteration}/"
    fi
fi

if [ ! -d "${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${PHY_SETTING_C}_iter${Iteration}" ]; then
    mkdir -p "${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${PHY_SETTING_C}_iter${Iteration}"
fi
# Only move if destination does NOT already have files
if ! ls ${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${PHY_SETTING_C}_iter${Iteration}/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
    if ls ${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${SIMC_Suffix_File}* 1> /dev/null 2>&1; then
        mv ${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${SIMC_Suffix_File}* "${VOLATILEPATH}/OUTPUT/Analysis/SIMC/${PHY_SETTING_C}_iter${Iteration}/"
    fi
fi

#################################################################################################################################

# Section for phi-bining and simc yield calculation Script
if [[ $s_flag == "true" ]]; then
    if [ ! -d "${UTILPATH}/LTSep_CSVs/simc_yields_csv/${PHY_SETTING_C}_iter${Iteration}" ]; then
    mkdir -p "${UTILPATH}/LTSep_CSVs/simc_yields_csv/${PHY_SETTING_C}_iter${Iteration}"
    fi
    while true; do
        read -p "Do you wish to do t & phi binning and calculate yields for simc setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion simc ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix} ${Iteration}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_simc_yield.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix} ${Iteration}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                simc_dir_yield="simc_yield"
                simc_dir_iter="${PHY_SETTING_C}_iter${Iteration}"
                if [ ! -d "$simc_dir_yield/$simc_dir_iter" ]; then
                    mkdir -p "$simc_dir_yield/$simc_dir_iter"
                    echo "Directory '$simc_dir_yield/$simc_dir_iter' created."
                else
                    echo "Directory '$simc_dir_yield/$simc_dir_iter' already exists."
                fi
                mv *Yield* "$simc_dir_yield/$simc_dir_iter"
                echo "Output files moved to $simc_dir_yield/$simc_dir_iter"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for Data and SIMC comparison Script
elif [[ $r_flag == "true" ]]; then
    if [ ! -d "${UTILPATH}/LTSep_CSVs/datasimc_ratios_csv/${PHY_SETTING_C}_iter${Iteration}" ]; then
        mkdir -p "${UTILPATH}/LTSep_CSVs/datasimc_ratios_csv/${PHY_SETTING_C}_iter${Iteration}"
    fi
    while true; do
        read -p "Do you wish to do Data & SIMC comparison and plotting for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Iteration}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_physics_ratio_comp.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${Iteration}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                ratio_dir="ratios"
                ratio_dir_iter="${PHY_SETTING_C}_iter${Iteration}"
                if [ ! -d "$ratio_dir/$ratio_dir_iter" ]; then
                    mkdir -p "$ratio_dir/$ratio_dir_iter"
                    echo "Directory '$ratio_dir/$ratio_dir_iter' created."
                else
                    echo "Directory '$ratio_dir/$ratio_dir_iter' already exists."
                fi
                mv *Ratio* "$ratio_dir/$ratio_dir_iter"
                echo "Output files moved to $ratio_dir/$ratio_dir_iter"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 3

#################################################################################################################################
# Section for Average Kinematics and Yield Calculation Script
elif [[ $a_flag == "true" ]]; then
    if [ ! -d "${UTILPATH}/LTSep_CSVs/avg_kinematics_csv/${PHY_SETTING_C}_iter${Iteration}" ]; then
        mkdir -p "${UTILPATH}/LTSep_CSVs/avg_kinematics_csv/${PHY_SETTING_C}_iter${Iteration}"
    fi
        if [ ! -d "${UTILPATH}/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_C}_iter${Iteration}" ]; then
        mkdir -p "${UTILPATH}/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_C}_iter${Iteration}"
    fi
    while true; do
        read -p "Do you wish to calculate average kinematics and yields for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${UTILPATH}/scripts/yield/src/pion_yield_avg_kin.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV} ${Iteration}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${UTILPATH}/scripts/yield/src/pion_yield_avg_kin.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV} ${Iteration}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                # Directory setup for ROOT files
                cd $REPLAYPATH/OUTPUT/Analysis/PionLT/
                avgkin_dir="avg_kin"
                avg_kin_dir_iter="${PHY_SETTING_C}_iter${Iteration}"
                if [ ! -d "$avgkin_dir/$avg_kin_dir_iter" ]; then
                    mkdir -p "$avgkin_dir/$avg_kin_dir_iter"
                    echo "Directory '$avgkin_dir/$avg_kin_dir_iter' created."
                else
                    echo "Directory '$avgkin_dir/$avg_kin_dir_iter' already exists."
                fi
                mv *avgkin* "$avgkin_dir/$avg_kin_dir_iter"
                echo "Output files moved to $avgkin_dir/$avg_kin_dir_iter"
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

exit 0