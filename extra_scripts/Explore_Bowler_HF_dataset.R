library(raster)

fig4 <- raster("./input_data/Bowler_SM/Bowler_ATCs_Figs_SOM/Figure_4.tif")
fig6 <- stack("./input_data/Bowler_SM/Bowler_ATCs_Figs_SOM/Figure_6.grd")

plot(fig4)
plot(fig6)

fig4@crs # Units in m
res(fig4) # Resolution = 25km x 25km

fig6@crs # Units in m
res(fig6) # Resolution = 10km x 10km

names(fig6)

plot(fig6[["Climate_change"]])

# Qualitative values with 5 levels, mixed between terrestrial and marine... Not very useful...
table(readAll(fig6[["Climate_change"]])@data@values)

# Qualitative values with 6 levels, mixed between terrestrial and marine... Not very useful...
table(readAll(fig6[["Human_use"]])@data@values)

# Binary raster with just the top 10% highest...
table(readAll(fig6[["Human_population"]])@data@values)

# Qualitative values with 5 levels, mixed between terrestrial and marine... Not very useful...
table(readAll(fig6[["Pollution"]])@data@values)

# Binary raster with just the top 10% highest...
table(readAll(fig6[["Alien_potential"]])@data@values)

# Qualitative values with 12 levels, mixed between terrestrial and marine... Not very useful...
table(readAll(fig6[["Cumulative"]])@data@values)

test <- calc(subset(fig6, 1:5), sum)

# Fig 6 is the sum of the previous 5 qualitative maps
plot(test)
plot(fig6[["Cumulative"]])
