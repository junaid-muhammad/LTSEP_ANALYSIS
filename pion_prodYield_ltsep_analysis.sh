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

while getopts 'hsraxipw' flag; do
    case "${flag}" in
        h)
        echo "-------------------------------------------------------------------"
        echo "./pion_prodYield_ltsep_analysis.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the physics analysis..."         
        echo "    -h, help"
        echo "    -s, run simc yield script to determine simc yields"
        echo "        run -> Flag=-s ITERATION=arg1 PHY_SETTING_INP=arg2"
        echo "    -r, run plotting script to compare Data/SIMC and determine ratios"
        echo "        run -> Flag=-p ITERATION=arg1 PHY_SETTING_INP=arg2"
        echo "    -a, run script to calculate average kinematics and yields for LTSep analysis "
        echo "        run -> Flag=-a ITERATION=arg1 Setting=arg2"
        echo "    -x, run script for cross section calculation"
        echo "        run -> Flag=-x ITERATION=arg1-(number) Setting=arg2"
        echo "    -i, run script for model ITERATIONs"
        echo "        run -> Flag=-i ITERATION=arg1-(number) Setting=arg2"
        echo "    -p, run script for plotting cross-sections for physics setting"
        echo "        run -> Flag=-p ITERATION=arg1-(number) Setting=arg2" 
        echo "    -w, run script for weight calculation for physics setting"
        echo "        run -> Flag=-w ITERATION=arg1-(number) Setting=arg2" 
        exit 0
        ;;
        s) s_flag='true' ;;
        r) r_flag='true' ;;
        a) a_flag='true' ;;
        x) x_flag='true' ;;
        i) i_flag='true' ;;
        p) p_flag='true' ;;
        w) w_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

echo "Starting analysis of Pion events"
echo "I take as arguments the -flag, ITERATION number, and physics setting"

# Input params - beam energy, run list and max number of events
ITERATION=$2
if [[ -z "$2" ]]; then
    echo "Please provide ITERATION number!"
    exit 1
else
    echo "Running iteration number: ${ITERATION}"
fi

PHY_SETTING_INP=$3
if [[ -z "$3" ]]; then
    echo "I need a Physics setting as input!"
    echo "Please provide a Physics setting as input"
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
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` 
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
LTSEP_DIR_PATH=${REPLAYPATH}/LTSEP_ANALYSIS
PHY_SETTING_INP_PATH=${REPLAYPATH}/UTIL_BATCH/InputPHY_SETTING_INPs/PionLT_2021_2022/${PHY_SETTING_INP}
LTSEP_INPCSV_PATH=${UTILPATH}/LTSep_CSVs
BKUP_DIR=${REPLAYPATH}/LTSEP_ANALYSIS/iterations

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

cd $REPLAYPATH
# Input Arguments for ITERATION Scripts
PHY_SETTING=$(echo "${PHY_SETTING_INP}" | awk -F'_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}')
SIMC_SETTING=$(echo "${PHY_SETTING_INP}" | awk -F'_' '{print $1 $2 $3 "_" $4 $5}')
DATA_Suffix=ProdCoin_Analysed_Data
DUMMY_Suffix=ProdCoin_Analysed_Dummy_Data
SIMC_Suffix="Prod_Coin_${SIMC_SETTING}"
DATA_RUN_LIST=${PHY_SETTING}
DUMMY_RUN_LIST=${PHY_SETTING}_dummy
CSV_FILE=PionLT_coin_production_Prod_efficiency_data_2025_03_08

# Input Arguments for avergae kinematics and yields calculation Script
PHY_SETTING_C=$(echo "${PHY_SETTING_INP}" | awk -F'_' '{print $1 "_" $2 "_" $3}')
SIMC_SETTING_C=$(echo "${PHY_SETTING_INP}" | awk -F'_' '{print $1 $2 $3}')
SIMC_Suffix_C="Prod_Coin_${SIMC_SETTING_C}"
RUN_LIST_C=${PHY_SETTING_C}_PHY_SETTING_INP
DATA_YIELD_CSV=Physics_Data_Yield
DATA_AVG_KIN_CSV=Physics_Avg_Data_Kinematics
SIMC_YIELD_CSV=Physics_SIMC_Yield
SIMC_AVG_KIN_CSV=Physics_Avg_SIMC_Kinematics

# Input Arguments for plotting cross-sections Script
Q2=$(echo "${PHY_SETTING}" | awk -F'_' '{print $1}')  # Extract Q3p85
Q2_val=$(echo "$Q2" | sed 's/Q\([0-9]\+\)p\([0-9]\+\)/\1\2/')
#echo "$Q2_val"
if [[ "${PHY_SETTING_C}" == "Q3p85_W2p62_t0p21" ]]; then
    LOWEPS=292
    HIGHEPS=779
else
    echo "Need to update epsilon values for ${PHY_SETTING_C}"
fi
XSECT_LEPS_DATA="${PHY_SETTING_C}_x_unsep_pi_${Q2_val}_${LOWEPS}"
XSECT_HEPS_DATA="${PHY_SETTING_C}_x_unsep_pi_${Q2_val}_${HIGHEPS}"

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create directories if it doesn't exist
RAW_SIMC_PATH=${VOLATILEPATH}/worksim
RECON_SIMC_PATH=${VOLATILEPATH}/OUTPUT/Analysis/SIMC
AVG_KIN_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/avg_kinematics_csv
RATIO_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/datasimc_ratios_csv
DIAMOND_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/diamond_cut_csv
LTSEP_INPUT_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/ltsep_input_csv
MM_OFFSET_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/mm_offset_cut_csv
PHYSICS_YIELDS_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/physics_yields_csv
SIMC_YIELDS_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/simc_yields_csv
T_BINNING_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/t_binning_csv
T_RESOLUTION_DIR=${LTSEP_DIR_PATH}/LTSep_CSVs/t_resolution_csv
INPUT_DIR=${BKUP_DIR}/${PHY_SETTING_C}_ltanalysis/input_files
OUTPUT_DIR=${BKUP_DIR}/${PHY_SETTING_C}_ltanalysis/output_files

directories_csv=("${RECON_SIMC_PATH}" "${RATIO_DIR}" "${LTSEP_INPUT_DIR}" "${SIMC_YIELDS_DIR}" "${AVG_KIN_DIR}" "${INPUT_DIR}" "${OUTPUT_DIR}")
for dir in "${directories_csv[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Creating directory: $dir"
        mkdir -p "$dir"
    fi
    # Create iteration-specific subdirectory if it doesn't exist
    if [ ! -d "$dir/${PHY_SETTING_C}_iter${ITERATION}" ]; then
        mkdir -p "$dir/${PHY_SETTING_C}_iter${ITERATION}"
    fi
done

# Copying required input CSVs from standard model calculations for further processing
csv_types=("avg_kinematics_csv" "datasimc_ratios_csv" "diamond_cut_csv" "ltsep_input_csv" "mm_offset_cut_csv" "physics_yields_csv" "simc_yields_csv" "t_binning_csv" "t_resolution_csv")
dest_dirs=("${AVG_KIN_DIR}" "${RATIO_DIR}" "${DIAMOND_DIR}" "${INPUT_DIR}" "${MM_OFFSET_DIR}" "${PHYSICS_YIELDS_DIR}" "${SIMC_YIELDS_DIR}" "${T_BINNING_DIR}" "${T_RESOLUTION_DIR}")
for i in "${!csv_types[@]}"; do
    csv_type="${csv_types[$i]}"
    dest_dir="${dest_dirs[$i]}"
    dest_subdir="${dest_dir}/${PHY_SETTING_C}_std"
    if [ ! -d "$dest_subdir" ]; then
        if compgen -G "${LTSEP_INPCSV_PATH}/${csv_type}/${PHY_SETTING_C}_std" > /dev/null; then
            cp -r ${LTSEP_INPCSV_PATH}/${csv_type}/${PHY_SETTING_C}_std "${dest_dir}/"
            echo "Copied ${csv_type} directory ${PHY_SETTING_C}_std to ${dest_dir}"
        fi
    fi
done

# Only move SIMC files for iteration
if [[ "${ITERATION}" == "00" ]]; then
    mkdir -p "${RAW_SIMC_PATH}/${PHY_SETTING_C}_iter${ITERATION}"
    # Define source and destination directories
    simc_dirs=("${RAW_SIMC_PATH}" "${RECON_SIMC_PATH}")
    # Move SIMC files to standard directories
    for base_dir in "${simc_dirs[@]}"; do
        dest_dir="${base_dir}/${PHY_SETTING_C}_iter${ITERATION}"
        # Create destination directory if it doesn't exist
        if [ ! -d "$dest_dir" ]; then
            mkdir -p "$dest_dir"
        fi
        # Only move if destination does NOT already have files
        if ! ls ${dest_dir}/${SIMC_Suffix_C}* 1> /dev/null 2>&1; then
            if ls ${base_dir}/${SIMC_Suffix_C}* 1> /dev/null 2>&1; then
                mv ${base_dir}/${SIMC_Suffix_C}* "$dest_dir/"
                echo "Moved SIMC files to $(basename "$base_dir")/${PHY_SETTING_C}_iter${ITERATION}"
            fi
        fi
    done
fi

#################################################################################################################################

# Section for phi-bining and simc yield calculation Script
if [[ $s_flag == "true" ]]; then
    if [ ! -d "${OUTPATH}/PionLT/simc_yield/${PHY_SETTING_C}_iter${ITERATION}" ]; then
        mkdir -p "${OUTPATH}/PionLT/simc_yield/${PHY_SETTING_C}_iter${ITERATION}"
        echo "Directory '${OUTPATH}/PionLT/simc_yield/${PHY_SETTING_C}_iter${ITERATION}' created."
    fi
    while true; do
        read -p "Do you wish to do t & phi binning and calculate yields for simc setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion simc ploting script
                if [ -f "${OUTPATH}/PionLT/simc_yield/${PHY_SETTING_C}_iter${ITERATION}/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        echo "Reprocessing"
                        python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_simc_yield_iter.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix} ${ITERATION}
                    else
                        echo "Skipping phi-binning and yield calculation script step"
                    fi
                elif [ ! -f  "${OUTPATH}/PionLT/simc_yield/${PHY_SETTING_C}_iter${ITERATION}/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_simc_yield_iter.py ${PHY_SETTING} ${MAXEVENTS} ${SIMC_Suffix} ${ITERATION}
                else echo "Pion coin output plots file already found in ${OUTPATH}/PionLT/simc_yield/ - Skipped python plotting script step"
                fi
	            )
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
    if [ ! -d "${OUTPATH}/PionLT/ratios/${PHY_SETTING_C}_iter${ITERATION}" ]; then
        mkdir -p "${OUTPATH}/PionLT/ratios/${PHY_SETTING_C}_iter${ITERATION}"
        echo "Directory '${OUTPATH}/PionLT/ratios/${PHY_SETTING_C}_iter${ITERATION}' created."
    fi
    while true; do
        read -p "Do you wish to do Data & SIMC comparison and plotting for physics setting ${PHY_SETTING}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/ratios/${PHY_SETTING_C}_iter${ITERATION}/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/ratios/${PHY_SETTING_C}_iter${ITERATION}/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_physics_ratio_comp_iter.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${ITERATION}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/${PHY_SETTING}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                        python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_physics_ratio_comp_iter.py ${PHY_SETTING} ${MAXEVENTS} ${DATA_Suffix} ${DUMMY_Suffix} ${SIMC_Suffix} ${DATA_RUN_LIST} ${DUMMY_RUN_LIST} ${CSV_FILE} ${ITERATION}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
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
    while true; do
        read -p "Do you wish to calculate average kinematics and yields for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                i=-1
                (
                # Section for Pion physics ploting script
                if [ -f "${UTILPATH}/OUTPUT/Analysis/PionLT/avg_kin/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                    read -p "Pion coin output plots file already exists, do you want to reprocess it? <Y/N> " option2
                    if [[ $option2 == "y" || $option2 == "Y" || $option2 == "yes" || $option2 == "Yes" ]]; then
                        rm "${UTILPATH}/OUTPUT/Analysis/PionLT/avg_kin/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root"
                        echo "Reprocessing"
                        python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_yield_avg_kin_iter.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV} ${ITERATION}
                    else
                        echo "Skipping  Data & SIMC comparison and plotting script step"
                    fi
                elif [ ! -f  "${UTILPATH}/OUTPUT/Analysis/PionLT/avg_kin/${PHY_SETTING_C}_${MAXEVENTS}_ProdCoin_Yield_Data.root" ]; then
                       python3 ${REPLAYPATH}/LTSEP_ANALYSIS/src/pion_yield_avg_kin_iter.py ${PHY_SETTING_C} ${MAXEVENTS} ${DATA_YIELD_CSV} ${DATA_AVG_KIN_CSV} ${SIMC_YIELD_CSV} ${SIMC_AVG_KIN_CSV} ${ITERATION}
                else echo "Pion coin output plots file already found in ${UTILPATH}/OUTPUT/Analysis/PionLT/ - Skipped python plotting script step"
                fi
	        )
                break;;
            [Nn]* )
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

sleep 1

###########################################################################################################################################################################################                                                                                   

# Section for cross-section script
elif [[ $x_flag == "true" ]]; then
    ITERATION="std"
    INPUT_STD_DIR="${BKUP_DIR}/${PHY_SETTING_C}_ltanalysis/input_files/${PHY_SETTING_C}_${ITERATION}"
    OUTPUT_STD_DIR="${BKUP_DIR}/${PHY_SETTING_C}_ltanalysis/output_files/${PHY_SETTING_C}_${ITERATION}"
    std_dirs=("${INPUT_STD_DIR}" "${OUTPUT_STD_DIR}")
    for dir in "${std_dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir" # Create the directory if it doesn't exist      
            echo "Directory '$dir' created."
        fi
    done
    cd $REPLAYPATH/LTSEP_ANALYSIS/src/
    # Copy files first
    lt_PHY_SETTING_INP_PATH="${REPLAYPATH}/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_C}_${ITERATION}"
    if compgen -G "${lt_PHY_SETTING_INP_PATH}/*" > /dev/null; then
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_"* "${INPUT_STD_DIR}/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_Eb"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_t_bin"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_list_setting"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_ave"* "averages/"
        echo "All necessary files and directories are ready for analysis."
    else
        echo "No files found in ${lt_PHY_SETTING_INP_PATH} to copy."
    fi

    # Compile Fortran code (once)
    read -p "Do you want to compile the default calc_xsect script? (y/n): " compile_std
    if [[ "$compile_std" =~ ^[Yy]$ ]]; then
        gfortran calc_xsect_std.f -o calc_xsect
        echo "Default calc_xsect compiled."
    else
        echo "Skipping compilation of default calc_xsect."
    fi

    # 3. Prompt and run
    while true; do
        read -p "Do you wish to calculate the cross-sections for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                if ./calc_xsect "${PHY_SETTING_C}" | tee "output/${PHY_SETTING_C}_calc_xsect_${ITERATION}_$(date +%d%b%Y).out"
                then
                    echo "Cross-section calculation completed for physics setting ${PHY_SETTING_C}."
                else
                    echo "Error: Calculation failed for physics setting ${PHY_SETTING_C}."
                fi
                # Ask user if they want to run the plotting script
                read -p "Do you want to run the Rosenbluth fitting and plotting script for cross-sections? (y/n): " run_plotting
                if [[ "$run_plotting" =~ ^[Yy]$ ]]; then
                    python3 plotting_xsec_std.py "${PHY_SETTING_C}" "${ITERATION}" "${XSECT_LEPS_DATA}" "${XSECT_HEPS_DATA}"
                    echo "Cross-section plotting completed."
                else
                    echo "Skipping plotting script."
                fi

                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_STD_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        mv output/${PHY_SETTING_C}_*std* "${OUTPUT_STD_DIR}"
                        mv xsects/${PHY_SETTING_C}_* "${OUTPUT_STD_DIR}"
                        mv plots/${PHY_SETTING_C}_*std* "${OUTPUT_STD_DIR}"
                        echo "Output files copied to ${OUTPUT_STD_DIR}."
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

# Section for model ITERATIONs
elif [[ $i_flag == "true" ]]; then
    cd $REPLAYPATH/LTSEP_ANALYSIS/src/
    # 1. Copy files first
    lt_PHY_SETTING_INP_PATH="${REPLAYPATH}/LTSEP_ANALYSIS/LTSep_CSVs/ltsep_input_csv/${PHY_SETTING_C}_iter${ITERATION}"
    if compgen -G "${lt_PHY_SETTING_INP_PATH}/*" > /dev/null; then
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_"* "${INPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_Eb"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_t_bin"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_list_setting"* "input/"
        cp "${lt_PHY_SETTING_INP_PATH}/${PHY_SETTING_C}_ave"* "averages/"
        echo "All necessary files and directories are ready for analysis."

    else
        echo "No files found in ${lt_PHY_SETTING_INP_PATH} to copy."
    fi

    # 2. Compile Fortran code (once)
    read -p "Do you want to compile the iteration calc_xsect script? (y/n): " compile_iter
    if [[ "$compile_iter" =~ ^[Yy]$ ]]; then
        gfortran calc_xsect_iter.f -o calc_xsect
        echo "Iteration calc_xsect compiled."
    else
        echo "Skipping compilation of iteration calc_xsect."
    fi

    # 3. Prompt and run
    while true; do
        read -p "Do you wish to run ITERATION ${ITERATION} for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                if ./calc_xsect "${ITERATION}" "${PHY_SETTING_C}" | tee "output/${PHY_SETTING_C}_calc_xsect_iter${ITERATION}_$(date +%Y%b%d).out"
                then
                    echo "Physics ITERATION ${ITERATION} completed for physics setting ${PHY_SETTING_C}."
                else
                    echo "Error: Physics ITERATION ${ITERATION} failed for physics setting ${PHY_SETTING_C}."
                    exit 1
                fi
                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        cp output/${PHY_SETTING_C}_*iter${ITERATION}* "${OUTPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
                        cp xsects/${PHY_SETTING_C}_* "${OUTPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
                        cp fit_params/iter${ITERATION}/* "${INPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
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
    cd $REPLAYPATH/LTSEP_ANALYSIS/src/
    while true; do
        read -p "Do you wish to plot the cross-sections for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                # Processing logic for Pion physics ltsep analysis
                echo "Running iteration plotting script for physics setting ${PHY_SETTING_C}..."
                if python3 plotting_xsec_iter.py "${PHY_SETTING_C}" "${ITERATION}" "${XSECT_LEPS_DATA}" "${XSECT_HEPS_DATA}"
                then
                    echo "Cross-section plotting completed for physics setting ${PHY_SETTING_C}."
                else
                    echo "Error: Plotting script failed for physics setting ${PHY_SETTING_C}."
                fi
                
                if compgen -G "output/*" > /dev/null; then
                    read -p "Are you sure you want to copy output files to ${OUTPUT_DIR}? (y/n): " confirm
                    if [[ "$confirm" =~ ^[Yy]$ ]]; then
                        mv xsects/${PHY_SETTING_C}_* "${OUTPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
                        mv output/*iter${ITERATION}* "${OUTPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
                        mv plots/${PHY_SETTING_C}_*iter${ITERATION}* "${OUTPUT_DIR}/${PHY_SETTING_C}_iter${ITERATION}/"
                        next_iter=$(printf "%02d" $((10#$ITERATION + 1)) 2>/dev/null)
                        if [ ! -d "fit_params/iter${next_iter}" ]; then
                            mkdir -p "fit_params/iter${next_iter}"
                        fi
                        for f in output/new_fitparams_iter${next_iter}_par.pl*; do
                            cp "$f" "fit_params/iter${next_iter}/par.pl_${Q2_val}"
                        done
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

# Section for SIMC weight recalculation
elif [[ $w_flag == "true" ]]; then
    cd $REPLAYPATH/LTSEP_ANALYSIS/src/
    while true; do
        read -p "Do you wish to re-calculate SIMC weights for physics setting ${PHY_SETTING_C}? (Please answer yes or no) " yn
        case $yn in
            [Yy]* )
                python3 reweight_simc.py "${PHY_SETTING_C}" "${ITERATION}" "${SIMC_Suffix}"
                if [ "${ITERATION}" == "00" ]; then
                    HISTPATH="${VOLATILEPATH}/OUTPUT/Analysis/SIMC"
                    cp  "${HISTPATH}/${PHY_SETTING_C}_itersimc/${SIMC_Suffix}.hist" "${HISTPATH}/${PHY_SETTING_C}_iter${ITERATION}/"
                else
                    ITERATION_PREV=$(printf "%02d" $((10#$ITERATION - 1)))
                    HISTPATH="${VOLATILEPATH}/OUTPUT/Analysis/SIMC"
                    cp  "${HISTPATH}/${PHY_SETTING_C}_iter${ITERATION_PREV}/${SIMC_Suffix}.hist" "${HISTPATH}/${PHY_SETTING_C}_iter${ITERATION}/"
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

fi

exit 0