#3D segmentation of cell-type specific nuclei from multi-color stainings

The purpose of this script is to automatically obtain ROIs (region of interest selections) of specific cell nuclei from multiple slices of an image stack for further analyses such as intensity measurements.
A general nuclear staining (i. e. DAPI) is combined with a cell type marker to restrict segmentation to nuclei of the cell type of interest.
The input images are two-channel hyperstacks, with the convention that C1 (red) is the nuclear channel and C2 (green) is the cell channel.
First, the cell stainings are segmented to obtain the sub-volume of the image containing the nuclei of interest, then the nuclei are segmented within the respective cellular volumes using the nuclear stainings.
The segmentation procedure for both sets of images involves the following steps:
1. Pre-processing: application of a Gaussian filter with a user-defined width
2. Thresholding with the ‘Auto threshold…’ command using a method of choice
3. Watershed segmentation (only for the nuclei)
4. Selection of ROIs within a user-defined range of size and shape properties with the ‘Analyze Particles’ command
The final output is saved as an ROI object in a directory of choice. To try different segmentation parameters, the script can be run in ‘test mode’, where the results are displayed for each image.
