#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2025-05-25 12:12:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################

while getopts 'hxipw' flag; do
    case "${flag}" in
        h)
        echo "-------------------------------------------------------------------"
        echo "./pion_prodYield_ltsep_analysis.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the physics analysis..."         
        echo "    -h, help"
        echo "    -x, run script for cross section calculation"
        echo "        run -> Flag=-x iteration=arg1-(number) Setting=arg2"
        echo "    -i, run script for model iterations"
        echo "        run -> Flag=-i iteration=arg1-(number) Setting=arg2"
        echo "    -p, run script for plotting cross-sections for physics setting"
        echo "        run -> Flag=-p iteration=arg1-(number) Setting=arg2" 
        echo "    -w, run script for weight calculation for physics setting"
        echo "        run -> Flag=-w iteration=arg1-(number) Setting=arg2" 
        exit 0
        ;;
        x) x_flag='true' ;;
        i) i_flag='true' ;;
        p) p_flag='true' ;;
        w) w_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

echo "Starting analysis of Pion events"
echo "I take as arguments the -flag, iteration number, and physics setting"

# Input params - beam energy, run list and max number of events
if [[ $x_flag == "true" ]]; then
    Iteration=$2
    echo "For default cross-section calculation, Provide iteration = 00"
    Iteration=std
elif [[ $p_flag == "true" ]]; then
    Iteration=$2
    if [[ "$Iteration" == "00" ]]; then
        Iteration=std
        echo "Running cross-section calculation for iteration ${Iteration}"
    fi
else
    Iteration=$2
    if [[ -z "$2" || "$2" -le 00 ]]; then
        echo "Please provide Iteration number. It should be greater than 00!"
        echo "Running cross-section calculation for iteration ${Iteration}"
        exit 1
    fi
fi

PHY_SETTING=$3
if [[ -z "$3" ]]; then
    echo "I need a Physics setting as input!"
    echo "Please provide a Physics setting as input"
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

#cd $REPLAYPATH

# Input Arguments for cross section calculation Script
PHY_SETTING_LTSEP=$(echo "${PHY_SETTING}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
SIMC_SETTING=$(echo "${PHY_SETTING}" | awk -F'_' '{print $1 $2 $3 "_" $4 $5}')
SIMC_Suffix="Prod_Coin_${SIMC_SETTING}"

# Input Arguments for plotting cross-sections Script
Q2=$(echo "${PHY_SETTING}" | awk -F'_' '{print $1}')  # Extract Q3p85
Q2_val=$(echo "$Q2" | sed 's/Q\([0-9]\+\)p\([0-9]\+\)/\1\2/')
#echo "$Q2_val"
if [[ "${PHY_SETTING_LTSEP}" == "Q3p85_W2p62_t0p21" ]]; then
    LOWEPS=292
    HIGHEPS=779
else
    echo "Need to update epsilon values for ${PHY_SETTING_LTSEP}"
fi
XSECT_LEPS_DATA="${PHY_SETTING_LTSEP}_x_unsep_pi_${Q2_val}_${LOWEPS}"
XSECT_HEPS_DATA="${PHY_SETTING_LTSEP}_x_unsep_pi_${Q2_val}_${HIGHEPS}"

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd $REPLAYPATH/UTIL_PION/scripts/ltsep_analysis/

# Copying input files to the ltsep_analysis directory fro further processing
BASE_DIR="iterations/${PHY_SETTING_LTSEP}_ltanalysis"
ITER_DIR="${BASE_DIR}/iter_${Iteration}"
INPUT_DIR="${ITER_DIR}/input_files"
OUTPUT_DIR="${ITER_DIR}/output_files"

# Create BASE_DIR if it doesn't exist
if [ ! -d "$BASE_DIR" ]; then
    echo "Creating base directory: $BASE_DIR"
    mkdir -p "$BASE_DIR"
fi
# Create ITER_DIR if it doesn't exist
if [ ! -d "$ITER_DIR" ]; then
    echo "Creating iteration directory: $ITER_DIR"
    mkdir -p "$ITER_DIR"
fi
# Create INPUT_DIR if it doesn't exist
if [ ! -d "$INPUT_DIR" ]; then
    echo "Creating input_files directory: $INPUT_DIR"
    mkdir -p "$INPUT_DIR"
fi
# Create OUTPUT_DIR if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output_files directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

###########################################################################################################################################################################################                                                                                   

# Section for cross-section script
if [[ $x_flag == "true" ]]; then
    # 1. Copy files first
    lt_inputFile="${REPLAYPATH}/UTIL_PION/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_LTSEP}_${Iteration}"
    avg_inputFile="src/averages"
    setting_inputFile="src/input"
    if compgen -G "${lt_inputFile}/*" > /dev/null; then
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_"* "${INPUT_DIR}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_Eb"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_t_bin"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_list_setting"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_ave"* "${avg_inputFile}/"
        echo "All necessary files and directories are ready for analysis."

    else
        echo "No files found in ${lt_inputFile} to copy."
    fi

    # 2. Compile Fortran code (once)
    cd $REPLAYPATH/UTIL_PION/scripts/ltsep_analysis/src/
    read -p "Do you want to compile the default calc_xsect script? (y/n): " compile_std
    if [[ "$compile_std" =~ ^[Yy]$ ]]; then
        gfortran calc_xsect_std.f -o calc_xsect
        echo "Default calc_xsect compiled."
    else
        echo "Skipping compilation of default calc_xsect."
    fi

    # 3. Prompt and run
    while true; do
        read -p "Do you wish to calculate the cross-sections for physics setting ${PHY_SETTING_LTSEP}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                if ./calc_xsect "${PHY_SETTING_LTSEP}" | tee "output/${PHY_SETTING_LTSEP}_calc_xsect_${Iteration}_$(date +%d%b%Y).out"
                then
                    echo "Cross-section calculation completed for physics setting ${PHY_SETTING_LTSEP}."
                else
                    echo "Error: Calculation failed for physics setting ${PHY_SETTING_LTSEP}."
                fi
                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        cp output/${PHY_SETTING_LTSEP}_*std* "../${OUTPUT_DIR}/"
                        cp xsects/${PHY_SETTING_LTSEP}_* "../${OUTPUT_DIR}/"
                        echo "Output files copied to ${OUTPUT_DIR}."
                    else
                        echo "Copy operation cancelled."
                    fi
                else
                    echo "No output files found to copy."
                fi
                break
                ;;
            [Nn]* )
                exit 0
                ;;
            * )
                echo "Please answer yes or no."
                ;;
        esac
    done

sleep 1
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Section for model iterations
elif [[ $i_flag == "true" ]]; then
    # 1. Copy files first
    lt_inputFile="${REPLAYPATH}/UTIL_PION/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_LTSEP}_iter${Iteration}"
    avg_inputFile="src/averages"
    setting_inputFile="src/input"
    if compgen -G "${lt_inputFile}/*" > /dev/null; then
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_"* "${INPUT_DIR}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_Eb"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_t_bin"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_list_setting"* "${setting_inputFile}/"
        cp "${lt_inputFile}/${PHY_SETTING_LTSEP}_ave"* "${avg_inputFile}/"
        echo "All necessary files and directories are ready for analysis."

    else
        echo "No files found in ${lt_inputFile} to copy."
    fi

    # 2. Compile Fortran code (once)
    cd $REPLAYPATH/UTIL_PION/scripts/ltsep_analysis/src/
    read -p "Do you want to compile the default calc_xsect script? (y/n): " compile_std
    if [[ "$compile_std" =~ ^[Yy]$ ]]; then
        gfortran calc_xsect_iter.f -o calc_xsect
        echo "Iteration calc_xsect compiled."
    else
        echo "Skipping compilation of default calc_xsect."
    fi

    # 3. Prompt and run
    while true; do
        read -p "Do you wish to run iteration ${Iteration} for physics setting ${PHY_SETTING_LTSEP}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                if ./calc_xsect "${Iteration}" "${PHY_SETTING_LTSEP}" | tee "output/${PHY_SETTING_LTSEP}_calc_xsect_iter${Iteration}_$(date +%Y%b%d).out"
                then
                    echo "Physics Iteration ${Iteration} completed for physics setting ${PHY_SETTING_LTSEP}."
                else
                    echo "Error: Physics iteration ${Iteration} failed for physics setting ${PHY_SETTING_LTSEP}."
                    exit 1
                fi
                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        cp output/${PHY_SETTING_LTSEP}_*iter${Iteration}* "../${OUTPUT_DIR}/"
                        cp xsects/${PHY_SETTING_LTSEP}_* "../${OUTPUT_DIR}/"
                        echo "Output files copied to ${OUTPUT_DIR}."
                    else
                        echo "Copy operation cancelled."
                    fi
                else
                    echo "No output files found to copy."
                fi
                break
                ;;
            [Nn]* )
                exit 0
                ;;
            * )
                echo "Please answer yes or no."
                ;;
        esac
    done

sleep 1

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Section for plotting cross-sections
elif [[ $p_flag == "true" ]]; then
    cd $REPLAYPATH/UTIL_PION/scripts/ltsep_analysis/src/
    while true; do
        read -p "Do you wish to plot the cross-sections for physics setting ${PHY_SETTING_LTSEP}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics ltsep analysis
                if [ $Iteration == "std" ]; then
                    python3 plotting_xsec_std.py "${PHY_SETTING_LTSEP}" "${Iteration}" "${XSECT_LEPS_DATA}" "${XSECT_HEPS_DATA}"; 
                else
                    python3 plotting_xsec_iter.py "${PHY_SETTING_LTSEP}" "${Iteration}" "${XSECT_LEPS_DATA}" "${XSECT_HEPS_DATA}";
                fi
                echo "Cross-section plotting completed for physics setting ${PHY_SETTING_LTSEP}."
                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        cp output/${PHY_SETTING_LTSEP}_*iter${Iteration}* "../${OUTPUT_DIR}/"
                        if [ "$Iteration" != "std" ]; then
                            cp fit_params/iter${Iteration}/par* "../${INPUT_DIR}/"
                            # Create next iteration directory in fit_params
                            next_iter=$(printf "%02d" $((10#$Iteration + 1)) 2>/dev/null)
                            if [ ! -d "fit_params/iter${next_iter}" ]; then
                                mkdir -p "fit_params/iter${next_iter}"
                            fi
                            for f in output/new_fitparams_iter${Iteration}_par.pl*; do
                                newname="${f#output/new_fitparams_iter${Iteration}_}"
                                mv "$f" "fit_params/iter${next_iter}/${newname}"
                            done
                        fi
                        echo "Output files copied to ${OUTPUT_DIR}."
                    else
                        echo "Copy operation cancelled."
                    fi
                else
                    echo "No output files found to copy."
                fi
                break
                ;;
            [Nn]* )
                exit 0
                ;;
            * )
                echo "Please answer yes or no."
                ;;
        esac
    done

sleep 1

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Section for plotting cross-sections
elif [[ $w_flag == "true" ]]; then
    cd $REPLAYPATH/UTIL_PION/scripts/ltsep_analysis/src/
    while true; do
        read -p "Do you wish to re-calculate SIMC weights for physics setting ${PHY_SETTING_LTSEP}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                python3 reweight_simc.py "${PHY_SETTING_LTSEP}" "${Iteration}" "${SIMC_Suffix}"
                ITERATION_PREV=$(printf "%02d" $((10#$Iteration - 1)))
                HISTPATH="${VOLATILEPATH}/OUTPUT/Analysis/SIMC"
                cd ${HISTPATH}
                cp  "${HISTPATH}/${PHY_SETTING_LTSEP}_iter${ITERATION_PREV}/${SIMC_Suffix}.hist" .
                break
                ;;
            [Nn]* )
                exit 0
                ;;
            * )
                echo "Please answer yes or no."
                ;;
        esac
    done

fi

exit 0