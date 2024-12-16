# Whole Liver Phase-Based R2 Mapping in Liver Iron Overload within a Breath-hold

## Authors of Code
- Daiki Tamada
- Ruvini Navaratna

## Contents
- **Figure1 - Figure5**: MATLAB and Python code to generate the plots and T2 maps presented in Figures 1 to 5 of the paper.
- **tools**: MATLAB code used for simulations and image reconstruction.
- **phantom_raw_T2_data**: Phantom and volunteer data/images.
  - **phantom_raw_T2_data/Phantom**: Measured T1 and T2 values of phantoms, and T2 maps acquired using Spin Echo (SE) R2 and Phase-Based (PB) R2 methods. These data are MATLAB array files (.mat).
  - **phantom_raw_T2_data/Volunteer**: R2 and R2* maps of volunteers. R2 maps are measured using PBR2 and saved as MATLAB array files. R2* maps are measured using the product IDEALIQ sequence and saved as DICOM files.

## Data Description

### Phantom Data (phantom_raw_T2_data/Phantom)
- Measured T1 and T2 Values: Includes the relaxation times of the phantoms used in the study.
- T2 Maps: Acquired using SE R2 and PBR2 methods.
- File Format: All data are saved as MATLAB array files (.mat).

### Volunteer Data (phantom_raw_T2_data/Volunteer)
- Raw PB images used for R2 mapping: Measured using the PBR2 method and saved as MATLAB array files.
- R2 Maps: Reconstructed from the raw PB images and saved as DICOM file. Unnecessary DICOM tags are removed using the CTP DICOM Anonymizer.
- R2* Maps: Measured using the product IDEALIQ sequence and saved as DICOM file. Unnecessary DICOM tags are removed using the CTP DICOM Anonymizer.

## Usage

### Generating Figures:
1. Navigate to the folder corresponding to the figure you wish to reproduce (Figure1 to Figure5).
2. Run the MATLAB scripts (Figure1.m-Figure5.m) provided to generate the plots and T2 maps as presented in the paper.

### Simulations and Image Reconstruction:
The tools directory contains MATLAB scripts needed for simulations and reconstructing images. Refer to the comments within each script for detailed instructions on usage.

## Requirements
- **MATLAB**: Version R2022b or newer (ensure compatibility with the scripts)
- **Python**: Version 3.10 with necessary libraries installed (e.g., NumPy, Matplotlib)

## License
Please see LICENSE before downloading the code.

## Contact
PLease let me know if you have any comments or suggestions.

Daiki Tamada (dtamada@wisc.edu)
