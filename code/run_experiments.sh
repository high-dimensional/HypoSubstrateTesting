#########################################################################
#                                      ###     ###    #######   ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                       #########     #######   #########
#                                                                       #
#                                                                       #
#                      HypoSubstrateTesting                             #
# Testing hypothesized neural substrates with lesions                   #
#                                                                       #
#                                                                       #
# This script contains all steps to run the full pipeline               #
# of the project.                                                       #
#                                                                       #
#                                                                       #
# Author: Patrik Bey, patrik.bey@ucl.ac.uk                              #
#                                                                       #
#                                                                       #
# last update: 2026/03/06.                                              #
#                                                                       #
#                                                                       #
#########################################################################

# ---- SET VARIABLES ---- #
CODE_DIR="/path/to/GitHub/HypoSubstrateTesting/code"

cd ${CODE_DIR}

DATA_DIR="/path/to/data"
LESION_DIR="${DATA_DIR}/lesions"
SPM_PATH="/path/to/spm"
MATLAB_CMD="matlab"

# ---- substrate maps ---- #

BASE_IMAGES=(
    "MDN.nii"
    "SM.nii"
    "DAN.nii"
    "VAN.nii"
    "AF.nii"
    "IFOF.nii"
    "ILF.nii"
    "SLF.nii"
)


# ---- simulation parameters ---- #

NCAS=300    # sample size
NITS=5000   # iterations
NLOW=0.3    # lower end noise level
NHIGH=0.6   # higher end noise level


# ---- load utility functions ---- #

source utility.sh



# ---- preprocessing morphology ---- #

run_step1_morphology


# ---- preprocessing features ---- #
run_step2_features 


# ---- lesion principal axis regression ---- #

REGRESSION_TARGETS=(
    "morphology"
)

run_step3_simulations


# ---- full experimental setup ---- #
# combining all models and substrates in a single run. This is the most comprehensive way to run the full pipeline, but also the most time consuming.

REGRESSION_TARGETS=(
    "targetfile"
)

for b in ${BASE_IMAGES[@]}; do
    # local overwrite of BASE_IMAGES and TARGET_FILES for each substrate map
    # for function consistency across models.
    BASE_IMAGES=($b)
    for t in ${BASE_IMAGES[@]}; do
        TARGET_FILES=($t)
        run_step3_simulations
    done
done


# ---- create semi-lollipop plot ---- #
# After running all simulation models, the relevant model-specific results may be concatenated in the SignificanceRatios.csv files for the specific models (i.e. net_morph: morphology based simulations using network neural substrates). 
# These can be used to create the semi-lollipop plot, which visualizes the significance ratios for each model and substrate.

python create_semi_lolli.py /mnt/h/HypoLesion/Data/UCL/net_morph/SignificanceRatios.csv --spoke-color 'darkgrey' --limit 65
python create_semi_lolli.py /mnt/h/HypoLesion/Data/UCL/net_true/SignificanceRatios.csv --spoke-color 'darkgrey' --limit 105
python create_semi_lolli.py /mnt/h/HypoLesion/Data/UCL/tract_morph/SignificanceRatios.csv --spoke-color 'darkgrey' --limit 65
python create_semi_lolli.py /mnt/h/HypoLesion/Data/UCL/tract_true/SignificanceRatios.csv --spoke-color 'darkgrey' --limit 105