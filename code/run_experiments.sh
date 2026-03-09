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
# of the project.                                                      #
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


# ---- true substrate regression ---- #

REGRESSION_TARGETS=(
    "targetfile"
)


for b in ${BASE_IMAGES[@]}; do
    BASE_IMAGES=($b)
    TARGET_FILES=($b)
    run_step3_simulations
done


# ---- false MDN substrate regression ---- #

REGRESSION_TARGETS=(
    "targetfile"
)

BASE_IMAGES=(
    "MDN.nii")

TARGET_FILES=(
    "RF.nii")

run_step3_simulations
