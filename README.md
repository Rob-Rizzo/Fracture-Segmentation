# Fracture-Segmentation
Code to segment fracture from images

The code first uses a small and large median filter to segment the voids, i.e. pores and fractures, in the image. The sizes of the two median filters should
preferentially be: (1) small median filter, smaller than the pore size; and (2) large median filter, larger than fracture widths. The resulting filtered image 
is the binarised.

Because we are interested in only the linear features, i.e. the fractures, using a ratio between the maximum axis length and the eccentricity values of the 
connected components, the code filters out the round-shaped features.

The test image attached here is a CT micrograph of a ceramics sample (kaolin matrix + quartz temper).

The main script calls two functions:

  (1) roseEqualArea.m --> which allows to plot equal area rose digrams of the fractures' trace angles
                          This function has been written by Dr. David Healy (Univerisity of Aberdeen) and it's part of the FracPaQ toolbox
                          (http://fracpaq.com/index.html)
                          
  (2) RDPsimplify.m --> which, by using theRamer-Douglas-Peucker algorithm for curve semplification, reduces the number of vertices in the fracture traces 
                        according to a specified tolerance.
                        This function has been written by Wolfgang Schwanghart on 13 July 2010
