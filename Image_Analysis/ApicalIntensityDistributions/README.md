#Measurement of radial intensity profiles from multi-color stainings of hair cell apical surfaces

The purpose of this script is to automatically obtain an average radial intensity profile across multiple cells of a protein of interest that is distributed on the apical surface (i.e. a planar cell polarity protein).
An apical surface marker (i. e. actin or spectrin) is used to detect the centroid positions of the apical surfaces of all cells within the 3D stack.
A stack of 60x60 images, each centered at a detected apical surface, is generated and saved together with an average intensity projection over all cells.
Subsequently, the intensity of the protein of interest is measured within a ring of specified radius around the center of the average intensity projection.
It is assumed that the polarity axis of the cells is aligned to the image axis.  
