
<p align='center'>
    <img src= 'img/Figure1.png' width="30%">

</p>


# *HypoSubstrateTesting*
Testing hypothesized neural substrates with lesions

## ABOUT

This code has been developed by the High Dimensional Neurology group at the [Queen Square Institute of Neurology UCL](https://www.ucl.ac.uk/brain-sciences/ion).

Authors: <br>
Parashkev Nachev; <br>
[Patrik Bey](https://github.com/PatrikBey)



### INPUT REQUIREMENTS
The framework assumes the following file structure:
```bash
DATA_DIR
    |_lesions
        |_lesion1.nii
        |_lesion2.nii
    |_masks
        |_mask1.nii
        |_mask2.nii
```

### SOFTWARE REQUIREMENTS
1. Matlab
2. SPM
3. BASH (optional for automated run script usage)

### STEPS

#### 1. Computing lesion mask morphological properties

```bash
# single function call
run_step1_morphology
```
#### 2. Computing lesion mask / target mask features
```bash
# single function call
run_step2_features 
```
#### 3. Running simulations
```bash
# ---- true substrate regression ---- #
# example use case for true substrate regression model.

REGRESSION_TARGETS=(
    "targetfile"
)

BASE_IMAGES=(
    "MDN.nii")

TARGET_FILES=(
    "MDN.nii")

run_step3_simulations

# ---- false substrate regression ---- #
# example use case for MDN.nii base image and RF.nii showcasing the false substrate regression model.

REGRESSION_TARGETS=(
    "targetfile"
)

BASE_IMAGES=(
    "MDN.nii")

TARGET_FILES=(
    "RF.nii")

run_step3_simulations

```
#### 4. Visualizations of significance ratios across different models
```bash
# exemplary call for true-substrate regression models
python create_semi_lolli.py /path/to/true-substrate/SignificanceRatios.csv --spoke-color 'darkgrey' --limit 105
```

## REFERENCES
Patrik Bey, Joe Mole, James K. Ruffle, Amy Nelson, Edgar Chan, Parashkev Nachev, and Lisa Cipolotti, Testing hypothesized neural substrates with lesions, (2026), *under review*