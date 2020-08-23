# Fracture-Segmentation
Code to segment fracture from images

The code first uses a small and large median filter to segment the voids, i.e. pores and fractures, in the image. The sizes of the two median filters should
preferentially be: (1) small median filter, smaller than the pore size; and (2) large median filter, larger than fracture widths. The resulting filtered image 
is the binarised.

Because we are interested in only the linear features, i.e. the fractures, using a ratio between the maximum axis length and the eccentricity values of the 
connected components, the code filters out the round-shaped features.

The example attached use a CT micrograph of a ceramics sample (kaolin matrix + quartz temper).
