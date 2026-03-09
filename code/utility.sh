#!/usr/bin/env bash
#########################################################################
#                                      ###     ###    #######   ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                       #########     #######   #########
#                                                                       #
#                                                                       #
#                               HypoSubstrateTesting                    #
# Testing hypothesized neural substrates with lesions                   #
# This script contains all utility functions for the                    #
# cluster based running of simulation experiments                       #
#                                                                       #
# Author: Patrik Bey, patrik.bey@ucl.ac.uk                              #
#                                                                       #
#                                                                       #
# last update: 2026/03/06.                                              #
#                                                                       #
#                                                                       #
#########################################################################



# ---- general utility functions ---- #

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

UTILITY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_matlab() {
    local matlab_code="$1"
    log "Running MATLAB code..."
    
    local matlab_script="addpath('$(pwd)'); addpath('$UTILITY_DIR'); try, $matlab_code; catch ME, fprintf('Error: %s\n', ME.message); fprintf('Stack trace:\n'); for k = 1:length(ME.stack), fprintf('  In %s at line %d\n', ME.stack(k).name, ME.stack(k).line); end; exit(1); end; exit(0);"
    
    $MATLAB_CMD -batch "$matlab_script"
    
    if [ $? -eq 0 ]; then
        log "MATLAB execution completed successfully"
    else
        log "ERROR: MATLAB execution failed"
        return 1
    fi
}


# ---- matlab wrapper functions ---- #


run_step1_morphology() {
    log "========================================"
    log "STEP 1: PREPROCESSING MORPHOLOGY"
    log "========================================"
    
    log "Computing morphological features for lesions in: $LESION_DIR"
    log "Output will be saved to: $DATA_DIR/$(basename $LESION_DIR)-morphology.mat"
    
    run_matlab "preprocessingMorphology('lesionDir', '$LESION_DIR', 'dataDir', '$DATA_DIR', 'spmPath', '$SPM_PATH');"
    
    if [ $? -eq 0 ]; then
        log "Step 1 complete: Morphology features computed"
    else
        log "ERROR: Step 1 failed"
        return 1
    fi
}


run_step2_features() {
    log "========================================"
    log "STEP 2: PREPROCESSING FEATURES"
    log "========================================"
    log "Processing ${#BASE_IMAGES[@]} base images"
    
    for base_image in "${BASE_IMAGES[@]}"; do
        log "----------------------------------------"
        log "Processing base image: $base_image"
        
        run_matlab "preprocessingFeatures('baseImageFile', '$base_image', 'lesionDir', '$LESION_DIR', 'dataDir', '$DATA_DIR', 'spmPath', '$SPM_PATH');"
        
        if [ $? -eq 0 ]; then
            log "Features computed for: $base_image"
        else
            log "ERROR: Failed to compute features for: $base_image"
        fi
    done
    
    log "Step 2 complete: All feature preprocessing done"
}



run_step3_simulations() {
    log "========================================"
    log "STEP 3: RUNNING SIMULATIONS"
    log "========================================"
    log "Simulations parameters:"
    log "  Cases per iteration: $NCAS"
    log "  Number of iterations: $NITS"
    log "  Noise range: $NLOW - $NHIGH"
    log ""
    
    # Loop over base images
    for base_image in "${BASE_IMAGES[@]}"; do
        log "----------------------------------------"
        log "Base image: $base_image"
        log "----------------------------------------"
        
        # Loop over regression targets
        for target in "${REGRESSION_TARGETS[@]}"; do
            log "Running simulations with target: $target"
            
            if [ "$target" == "targetfile" ]; then
                # For targetfile mode, loop over target files
                for target_file in "${TARGET_FILES[@]}"; do
                    log "  Using target file: $target_file"
                    
                    run_matlab "simulations('baseImageFile', '$base_image', 'regressionTarget', 'targetfile', 'targetFile', '$target_file', 'lesionDir', '$LESION_DIR', 'dataDir', '$DATA_DIR', 'spmPath', '$SPM_PATH', 'ncas', $NCAS, 'nlow', $NLOW, 'nhigh', $NHIGH, 'nits', $NITS, 'figureFormat', 'pdf');"
                    
                    if [ $? -eq 0 ]; then
                        log "  Simulations complete for $base_image + $target_file"
                    else
                        log "  ERROR: Simulations failed for $base_image + $target_file"
                    fi
                done
            else
                # For morphology or baseimage mode
                run_matlab "simulations('baseImageFile', '$base_image', 'regressionTarget', '$target', 'lesionDir', '$LESION_DIR', 'dataDir', '$DATA_DIR', 'spmPath', '$SPM_PATH', 'ncas', $NCAS, 'nlow', $NLOW, 'nhigh', $NHIGH, 'nits', $NITS, 'figureFormat', 'pdf');"
                
                if [ $? -eq 0 ]; then
                    log "Simulations complete for $base_image + $target"
                else
                    log "ERROR: Simulations failed for $base_image + $target"
                fi
            fi
        done
    done
    
    log "Step 3 complete: All simulations done"
}

