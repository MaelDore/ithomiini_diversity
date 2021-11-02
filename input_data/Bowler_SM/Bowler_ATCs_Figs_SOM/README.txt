Spatial files to figures 4 and 6 of:
https://www.biorxiv.org/content/10.1101/432880v3

Both files can be read into R using the raster package:
library(raster)
fig4 <- raster(Figure_4.tif)
fig6 <- stack(Figure_6.grd)#this is a raster stack and requires the auxillary file Figure_6.gri to be in the same folder.
names(fig6)#to see names of layers - fig6[["Cumulative"]] extracts the main plot in Fig. 6

Data can be used for any purpose, with citation to the preprint (or better the manuscript, when it has been published by a journal).

Any queries to be directed to: diana.bowler@gmail.com


