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
# This script contains all steps to run the full pipeline described in  #
# (REF).                                                                #
# Lesions and substrate masks are assumed to be on one common grid.     #
# Steps:                                                                #
# 1. brain check of the lesions (masking if lesion voxels lie outside   #
#    the brain)                                                         #
# 2. lesion features                                                    #
# 3. substrate overlaps                                               #
# 4. simulations                                                        #
# 5. scenario figures                                                   #
# 6. radial summary figure                                              #
#                                                                       #
# input:                                                                #
# LESION_DIR, MASK_DIR, BRAIN_MASK, SPM_PATH (environment variables)    #
# FORCE=1 (optional): recompute steps whose output files exist already  #
# otherwise only gnerate missing output files                           #
#                                                                       #
#                                                                       #
# Author: Patrik Bey, patrik.bey@ucl.ac.uk                              #
#                                                                       #
#                                                                       #
# last update: 2026/11/11.                                              #
#                                                                       #
#                                                                       #
#########################################################################


set -euo pipefail


# ---- SET VARIABLES ---- #

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LESION_DIR="${LESION_DIR:?set LESION_DIR to the folder with the lesion .nii files}"
MASK_DIR="${MASK_DIR:?set MASK_DIR to the folder with the substrate masks (<name>.nii)}"
BRAIN_MASK="${BRAIN_MASK:?set BRAIN_MASK to the brain mask (.nii or .nii.gz)}"
SPM_PATH="${SPM_PATH:?set SPM_PATH to the SPM folder}"
OUT_DIR="${OUT_DIR:-${CODE_DIR}/output}"
MATLAB_CMD="${MATLAB_CMD:-matlab}"
PYTHON="${PYTHON:-python3}"
MAX_JOBS="${MAX_JOBS:-4}"
N_PARTS="${N_PARTS:-${MAX_JOBS}}"
FORCE="${FORCE:-0}"

# ---- simulation parameters ---- #

NITS=5000
NLOW=0.3
NHIGH=0.6

# hypothesized compartments (RF is only used as a generator)
COMPARTMENTS=(MDN SM DAN VAN SLF ILF AF IFOF)
OVERLAP_MASKS=("${COMPARTMENTS[@]}" RF)

CLEAN_DIR="${OUT_DIR}/lesions_brain"
FEATURE_DIR="${OUT_DIR}/features"
RESULT_DIR="${OUT_DIR}/results"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${FEATURE_DIR}" "${RESULT_DIR}" "${LOG_DIR}"


# ---- utility functions ---- #

log_msg() {
    # print out text for logging
    _type=$( echo ${1} | cut -d'|' -f1 )
    _message=${1}
    if [[ ${_type,,} = "start " ]] || [[ ${_type,,} = "finished " ]] || [[ ${_type,,} = "error " ]] || [[ ${_type,,} = "warning " ]]; then
            echo -e "\n$(date) $(basename  -- "$0") | ${_message}" | lolcat
        else
            echo -e "\n$(date) $(basename  -- "$0") | ${_message}"
    fi
}

# true if the output file exists and recomputation is not forced
exists() { [ "${FORCE}" != "1" ] && [ -f "$1" ]; }

run_matlab() {
    "${MATLAB_CMD}" -singleCompThread -batch "addpath('${CODE_DIR}'); try, $1; catch ME, disp(getReport(ME)); exit(1); end"
}

# run in the background, at most MAX_JOBS at a time, staggered to avoid simultaneous MATLAB start-ups
PIDS=()
spawn() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do wait -n || true; done
    sleep 3
    "$@" &
    PIDS+=($!)
}
wait_all() {
    local failed=0
    for pid in "${PIDS[@]:-}"; do [ -n "${pid}" ] && { wait "${pid}" || failed=1; }; done
    PIDS=()
    return ${failed}
}


# ---- preprocessing: brain check and masking of the lesions ---- #

if exists "${LOG_DIR}/brain_mask_status.txt"; then
    log_msg "| UPDATE | Step 1: brain check found in ${LOG_DIR}, skipping |"
else
    log_msg "| START | Step 1: checking lesions for voxels outside the brain |"
    "${PYTHON}" "${CODE_DIR}/mask_lesions_to_brain.py" --lesionDir "${LESION_DIR}" --brainMask "${BRAIN_MASK}" \
        --outDir "${CLEAN_DIR}" --logDir "${LOG_DIR}" --jobs "${MAX_JOBS}" 2>&1 | tee "${LOG_DIR}/brain_mask.log"
fi
case "$(cat "${LOG_DIR}/brain_mask_status.txt")" in
    MASKED) [ -d "${CLEAN_DIR}" ] || { log_msg "| ERROR | masked lesions missing: ${CLEAN_DIR} (rerun with FORCE=1) |"; exit 1; }
            LESIONS_USED="${CLEAN_DIR}"; log_msg "| UPDATE | Lesions were brain-masked; using ${LESIONS_USED} |" ;;
    CLEAN)  LESIONS_USED="${LESION_DIR}"; log_msg "| UPDATE | No lesion voxels outside the brain; masking skipped |" ;;
    *)      log_msg "| ERROR | unreadable brain mask status |"; exit 1 ;;
esac


# ---- lesion features ---- #

if exists "${FEATURE_DIR}/lesion_features.mat"; then
    log_msg "| UPDATE | Step 2: lesion_features.mat found, skipping |"
else
    log_msg "| START | Step 2: lesion features (${N_PARTS} parts) |"
    rm -f "${FEATURE_DIR}"/lesion_features_part*.mat   # stale parts from earlier runs would corrupt the merge
    for k in $(seq 1 "${N_PARTS}"); do
        spawn run_matlab "lesion_features('lesionDir','${LESIONS_USED}','outDir','${FEATURE_DIR}','spmPath','${SPM_PATH}','part',${k},'nParts',${N_PARTS})"
    done
    wait_all || { log_msg "| ERROR | lesion features failed |"; exit 1; }
    # with a single part, lesion_features.m already writes lesion_features.mat; nothing to merge
    if [ "${N_PARTS}" -gt 1 ]; then
        run_matlab "lesion_features('lesionDir','${LESIONS_USED}','outDir','${FEATURE_DIR}','nParts',${N_PARTS},'merge',true)" || { log_msg "| ERROR | merging lesion features failed |"; exit 1; }
    fi
fi


# ---- substrate overlaps ---- #

log_msg "| START | Step 3: substrate overlaps |"
for s in "${OVERLAP_MASKS[@]}"; do
    if exists "${FEATURE_DIR}/overlap_${s}.mat"; then
        log_msg "| UPDATE | overlap_${s}.mat found, skipping |"
        continue
    fi
    spawn run_matlab "substrate_overlap('substrate','${s}','maskDir','${MASK_DIR}','lesionDir','${LESIONS_USED}','outDir','${FEATURE_DIR}','spmPath','${SPM_PATH}')"
done
wait_all || { log_msg "| ERROR | overlap computation failed |"; exit 1; }


# ---- simulations ---- #

# MDN hypothesized: global (false positive), MDN (true positive), RF (false positive)
# other compartment: global (false positive) and themselves (true positive)
log_msg "| START | Step 4: simulations |"
SCENARIOS=("MDN:global" "MDN:MDN" "MDN:RF")
for s in "${COMPARTMENTS[@]}"; do
    [ "${s}" = "MDN" ] && continue
    SCENARIOS+=("${s}:global" "${s}:${s}")
done
for sc in "${SCENARIOS[@]}"; do
    if exists "${RESULT_DIR}/hyp-${sc%%:*}_true-${sc##*:}_summary.csv"; then
        log_msg "| UPDATE | scenario ${sc} found, skipping |"
        continue
    fi
    spawn run_matlab "run_rldr_simulation('hypothesized','${sc%%:*}','generator','${sc##*:}','featureDir','${FEATURE_DIR}','outDir','${RESULT_DIR}','nits',${NITS},'nlow',${NLOW},'nhigh',${NHIGH})"
done
wait_all || { log_msg "| ERROR | simulations failed |"; exit 1; }


# ---- scenario figures ---- #

log_msg "| START | Step 5: scenario figures |"
for f in "${RESULT_DIR}"/hyp-*.mat; do
    if exists "${f%.mat}_figure.pdf"; then
        log_msg "| UPDATE | figure for $(basename "${f}") found, skipping |"
        continue
    fi
    spawn run_matlab "plot_rldr_figure('${f}')"
done
wait_all || { log_msg "| ERROR | figure generation failed |"; exit 1; }


# ---- radial lollipop figure ---- #

log_msg "| START | Step 6: radial lollipop figure |"
if exists "${RESULT_DIR}/Figure5.pdf"; then
    log_msg "| UPDATE | Figure5.pdf found, skipping |"
else
    "${PYTHON}" "${CODE_DIR}/plot_summary_lollipop.py" --summaryDir "${RESULT_DIR}" --out "${RESULT_DIR}/Figure5.pdf"
fi

log_msg "| FINISHED | Done. Results in ${OUT_DIR} |"
