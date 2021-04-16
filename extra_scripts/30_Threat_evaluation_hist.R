
###### Script 30: Species/mimicry ring-level threat evaluation ######

# Compute %area under threat (based on HII top 25% and 5%) for each OMU/species/ring
# Make histogram and density plot summarizing the results

### Inputs: 
   # List of OMU
   # Probability stack for OMU/species/rings from script 10/12/13
   # Mask for top 5% and 25% threats from script 19

### Outputs:
   # Summary table
   # Histogram and density plot
   # Quantiles values

   # Density curve with transparency. Per type of rings/species/OMUs based on the regions were they are found
       # Make mask per regions and label OMU/ring/species if presence in only 1 or several regions


### 0/ Prepare stuff ####

# Effacer l'environnement
rm(list = ls())

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Create summary table for rings
list.rings <- as.data.frame(table(list.models$Mimicry.model)) # Table of nb of OMU per ring
names(list.rings) <- c("ring", "nb_OMU")

# Load probability stacks
OMU_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"))
nlayers(OMU_proba_stack) # 783 OMU

sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
nlayers(sp_proba_stack) # 388 species

ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))
nlayers(ring_proba_stack) # 44 rings


# Load threat masks
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_quantiles_files.rds") 


##### 1/ Extract regional information for OMU/sp/ring #####


# Add region_type (or get it from automatic procedure based on bioregions masks?)

# Make a summary table per regions ?
  # Nb of OMU/species/ring
  # % under threat
  # Quantiles for threats for OMU/species/ring


##### 2/ Compute area threatened #####

plot(HII_mask_rank_5)
plot(HII_mask_rank_25)

### 2.1/ For OMUs ####
OMU_proba_brick <- readAll(brick(OMU_proba_stack))

# 2.1.1/ Compute full area
area_full <- apply(X = OMU_proba_brick@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_full <- area_full[match(as.character(list.models$Tag.model), names(OMU_proba_stack))] # Reorder in the same order than list.models
list.models$area_full <- area_full # Save in table

# 2.1.2/ Compute area under threat for top 25% threat
threat_25_OMU <- readAll(HII_mask_rank_25*OMU_proba_stack) # Crop distri with threat top 25%

area_threat_25 <- apply(X = threat_25_OMU@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_25 <- area_threat_25[match(as.character(list.models$Tag.model), names(OMU_proba_stack))] # Reorder in the sam eorder than list.models
list.models$area_threat_25 <- area_threat_25 # Save in table
list.models$area_threat_25_perc <- list.models$area_threat_25/list.models$area_full*100 # Compute % of threatened area

# 2.1.3/ Compute area under threat for top 5% threat
threat_5_OMU <- readAll(HII_mask_rank_5*OMU_proba_stack) # Crop distri with threat top 5%

area_threat_5 <- apply(X = threat_5_OMU@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_5 <- area_threat_5[match(as.character(list.models$Tag.model), names(OMU_proba_stack))] # Reorder in the sam eorder than list.models
list.models$area_threat_5 <- area_threat_5 # Save in table
list.models$area_threat_5_perc <- list.models$area_threat_5/list.models$area_full*100 # Compute % of threatened area

rm(OMU_proba_brick, OMU_proba_stack, threat_25_OMU, threat_5_OMU)


### 2.2/ For species ####
sp_proba_brick <- readAll(brick(sp_proba_stack))

# 2.2.1/ Compute full area
area_full <- apply(X = sp_proba_brick@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_full <- area_full[match(as.character(list.sp$Sp_full), names(sp_proba_stack))] # Reorder in the same order than list.sp
list.sp$area_full <- area_full # Save in table

# 2.2.2/ Compute area under threat for top 25% threat
threat_25_sp <- readAll(HII_mask_rank_25*sp_proba_stack) # Crop distri with threat top 25%

area_threat_25 <- apply(X = threat_25_sp@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_25 <- area_threat_25[match(as.character(list.sp$Sp_full), names(sp_proba_stack))] # Reorder in the sam eorder than list.sp
list.sp$area_threat_25 <- area_threat_25 # Save in table
list.sp$area_threat_25_perc <- list.sp$area_threat_25/list.sp$area_full*100 # Compute % of threatened area

# 2.2.3/ Compute area under threat for top 5% threat
threat_5_sp <- readAll(HII_mask_rank_5*sp_proba_stack) # Crop distri with threat top 5%

area_threat_5 <- apply(X = threat_5_sp@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_5 <- area_threat_5[match(as.character(list.sp$Sp_full), names(sp_proba_stack))] # Reorder in the sam eorder than list.sp
list.sp$area_threat_5 <- area_threat_5 # Save in table
list.sp$area_threat_5_perc <- list.sp$area_threat_5/list.sp$area_full*100 # Compute % of threatened area

rm(sp_proba_brick, sp_proba_stack, threat_25_sp, threat_5_sp)


### 2.3/ For rings ####
ring_proba_brick <- brick(ring_proba_stack)

# 2.3.1/ Compute full area
area_full <- apply(X = ring_proba_brick@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_full <- area_full[match(as.character(list.rings$ring), names(ring_proba_stack))] # Reorder in the sam eorder than list.rings
list.rings$area_full <- area_full # Save in table

# 2.3.2/ Compute area under threat for top 25% threat
threat_25_ring <- readAll(HII_mask_rank_25*ring_proba_stack) # Crop distri with threat top 25%

area_threat_25 <- apply(X = threat_25_ring@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_25 <- area_threat_25[match(as.character(list.rings$ring), names(ring_proba_stack))] # Reorder in the sam eorder than list.rings
list.rings$area_threat_25 <- area_threat_25 # Save in table
list.rings$area_threat_25_perc <- list.rings$area_threat_25/list.rings$area_full*100 # Compute % of threatened area

# 2.3.3/ Compute area under threat for top 5% threat
threat_5_ring <- readAll(HII_mask_rank_5*ring_proba_stack) # Crop distri with threat top 5%

area_threat_5 <- apply(X = threat_5_ring@data@values, MARGIN = 2, FUN = sum, na.rm = T)*772.84/1000 # Compute area in 10^3 km²
area_threat_5 <- area_threat_5[match(as.character(list.rings$ring), names(ring_proba_stack))] # Reorder in the same order than list.rings
list.rings$area_threat_5 <- area_threat_5 # Save in table
list.rings$area_threat_5_perc <- list.rings$area_threat_5/list.rings$area_full*100 # Compute % of threatened area

rm(ring_proba_brick, ring_proba_stack, threat_25_ring, threat_5_ring)


##### 3/ Plot histogram #####

### 3.1/ For species. Top 25% as threatened ####

pdf(file = "./graphs/Threat_areas/Species/Perc_threatened_areas_sp_25.pdf", height = 6.5, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = list.sp$area_threat_25_perc, 
     breaks = 20,
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     # main = "Proportion of species distribution range\nunder threat (Top 25%)",
     main = "Species: Top 25% threat   ",
     xlab = "Threatened area (%)",
     axes = T,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.8, lwd = 2)

d <- density(list.sp$area_threat_25_perc, adjust = 0.75)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*1700, col="black", lwd=2) # add a density estimate with defaults
abline(v = 25, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.sp$area_threat_25_perc), 1), nsmall = 1)

nb_over_25 <- sum(list.sp$area_threat_25_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.sp$area_threat_25_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_75 <- sum(list.sp$area_threat_25_perc > 75)
perc_over_75 <- format(round(nb_over_75/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.sp$area_threat_25_perc == 0)
perc_0 <- format(round(nb_0/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.sp$area_threat_25_perc == 100)
perc_100 <- format(round(nb_100/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " sp. (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " sp. (", perc_100,"%)"), 
                  paste0("> 25% = ", nb_over_25, " sp. (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " sp. (", perc_over_50,"%)"),
                  paste0("> 75% = ", nb_over_75, " sp. (", perc_over_75,"%)")), 
       x = "topright", inset = c(0, 0.03), 
       cex = 1.2, bty = "n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

# Names of species with 100% of their range threatened
list.sp$Sp_full[list.sp$area_threat_25_perc == 100]

### 3.2/ For species. Top 5% as threatened ####

pdf(file = "./graphs/Threat_areas/Species/Perc_threatened_areas_sp_5.pdf", height = 6.5, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = list.sp$area_threat_5_perc, 
     breaks = seq(0, 100, 10/3),
     freq = TRUE, col = "gray",
     xlim = c(0, 60),
     # ylim = c(0, 400),
     # main = "Proportion of species distribution range\nunder threat (Top 5%)",
     main = "Species: Top 5% threat     ",
     xlab = "Threatened area (%)",
     axes = T,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.8, lwd = 2)

d <- density(list.sp$area_threat_5_perc, adjust = 1)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*1200, col="black", lwd=2) # add a density estimate with defaults
abline(v = 5, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.sp$area_threat_5_perc), 1), nsmall = 1)

nb_over_5 <- sum(list.sp$area_threat_5_perc > 5)
perc_over_5 <- format(round(nb_over_5/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_25 <- sum(list.sp$area_threat_5_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.sp$area_threat_5_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)


nb_0 <- sum(list.sp$area_threat_5_perc == 0)
perc_0 <- format(round(nb_0/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.sp$area_threat_5_perc == 100)
perc_100 <- format(round(nb_100/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " sp. (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " sp. (", perc_100,"%)"), 
                  paste0("> 5% = ", nb_over_5, " sp. (", perc_over_5,"%)"),
                  paste0("> 25% = ", nb_over_25, " sp. (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " sp. (", perc_over_50,"%)")
                  ), 
       x = "topright", inset = c(0, 0.03), 
       cex = 1.2, bty = "n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

# Names of species with 100% of their range threatened
list.sp$Sp_full[list.sp$area_threat_5_perc == 100]

# Names of species with > 50% of their range threatened
list.sp$Sp_full[list.sp$area_threat_5_perc > 50]
# "Eutresis.dilucida"  "Napeogenes.benigna" "Pteronymia.asopo"   "Pteronymia.nsp1" 

### 3.3/ For rings. Top 25% as threatened ####

pdf(file = "./graphs/Threat_areas/Rings/Perc_threatened_areas_ring_25.pdf", height = 6.5, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = list.rings$area_threat_25_perc, 
     breaks = 20,
     freq = TRUE, col = "gray",
     xlim = c(0, 100),
     # ylim = c(0, 400),
     # main = "Proportion of mimicry ring distribution range\nunder threat (Top 25%)",
     main = "Mimicry rings: Top 25% threat     ",
     xlab = "Threatened area (%)",
     axes = T,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.8, lwd = 2)

d <- density(list.rings$area_threat_25_perc, adjust = 0.75)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*200, col="black", lwd=2) # add a density estimate with defaults
abline(v = 25, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.rings$area_threat_25_perc), 1), nsmall = 1)

nb_over_25 <- sum(list.rings$area_threat_25_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.rings$area_threat_25_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_75 <- sum(list.rings$area_threat_25_perc > 75)
perc_over_75 <- format(round(nb_over_75/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.rings$area_threat_25_perc == 0)
perc_0 <- format(round(nb_0/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.rings$area_threat_25_perc == 100)
perc_100 <- format(round(nb_100/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " ring (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " ring (", perc_100,"%)"), 
                  paste0("> 25% = ", nb_over_25, " rings (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " rings (", perc_over_50,"%)"),
                  paste0("> 75% = ", nb_over_75, " ring (", perc_over_75,"%)")), 
       x = "topright", inset = c(0, 0.03), 
       cex = 1.2, bty = "n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

# Names of mimicry rings with > 75% of their range threatened
list.rings$ring[list.rings$area_threat_25_perc > 75]
# HEMIXANTHE

### 3.4/ For rings. Top 5% as threatened ####

pdf(file = "./graphs/Threat_areas/Rings/Perc_threatened_areas_ring_5.pdf", height = 6.5, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = list.rings$area_threat_5_perc, 
     breaks = seq(0, 20, 1),
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     # main = "Proportion of mimicry ring distribution range\nunder threat (Top 5%)",
     main = "Mimicry rings: Top 5% threat     ",
     xlab = "Threatened area (%)",
     axes = T,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.8, lwd = 2)

d <- density(list.rings$area_threat_5_perc, adjust = 1)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*50, col="black", lwd=2) # add a density estimate with defaults
abline(v = 5, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.rings$area_threat_5_perc), 1), nsmall = 1)

nb_over_5 <- sum(list.rings$area_threat_5_perc > 5)
perc_over_5 <- format(round(nb_over_5/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_15 <- sum(list.rings$area_threat_5_perc > 15)
perc_over_15 <- format(round(nb_over_15/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.rings$area_threat_5_perc == 0)
perc_0 <- format(round(nb_0/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)


legend(legend = c(paste0("0% threat = ", nb_0, " rings (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("> 5% = ", nb_over_5, " rings (", perc_over_5,"%)"),
                  paste0("> 15% = ", nb_over_15, " rings (", perc_over_15,"%)")
                  ), 
x = "topright", inset = c(0, 0.03), 
cex = 1.2, bty = "n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

# Names of mimicry rings with > 15% of their range threatened
list.rings$ring[list.rings$area_threat_5_perc > 15]
# HEMIXANTHE HUMBOLDT LYSIMNIA

##### 4/ Multiple plot with all 4 plots: Species/rings vs. level of threats #####

pdf(file = "./graphs/Threat_areas/Perc_threatened_areas_all_options.pdf", height = 13, width = 15)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.5,7,3,1), mfrow = c(2,2))

# A/ For species. Top 25% as threatened

hist(x = list.sp$area_threat_25_perc, 
     breaks = 20,
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     # main = "Proportion of species distribution range\nunder threat (Top 25%)",
     main = "Species: Top 25% threat   ",
     xlab = "Threatened area (%)",
     ylab = "",
     axes = T,
     cex.axis = 2.0, cex.lab = 2.1, cex.main = 2.0, lwd = 2)

mtext(text = "Frequency", side = 2, cex = 1.9, line = 3.5)

d <- density(list.sp$area_threat_25_perc, adjust = 0.75)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*1700, col="black", lwd=2) # add a density estimate with defaults
abline(v = 25, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.sp$area_threat_25_perc), 1), nsmall = 1)

nb_over_25 <- sum(list.sp$area_threat_25_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.sp$area_threat_25_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_75 <- sum(list.sp$area_threat_25_perc > 75)
perc_over_75 <- format(round(nb_over_75/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.sp$area_threat_25_perc == 0)
perc_0 <- format(round(nb_0/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.sp$area_threat_25_perc == 100)
perc_100 <- format(round(nb_100/length(list.sp$area_threat_25_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " sp. (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " sp. (", perc_100,"%)"), 
                  paste0("> 25% = ", nb_over_25, " sp. (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " sp. (", perc_over_50,"%)"),
                  paste0("> 75% = ", nb_over_75, " sp. (", perc_over_75,"%)")), 
       x = "topright", inset = c(0, 0.03), y.intersp = 1.05, 
       cex = 1.8, bty = "n")

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 2.3, inset=c(0.01, 0.10))

# B/ For species. Top 5% as threatened

hist(x = list.sp$area_threat_5_perc, 
     breaks = seq(0, 100, 10/3),
     freq = TRUE, col = "gray",
     xlim = c(0, 60),
     ylim = c(0, 200),
     # main = "Proportion of species distribution range\nunder threat (Top 5%)",
     main = "Species: Top 5% threat     ",
     xlab = "Threatened area (%)",
     ylab = "",
     axes = T,
     cex.axis = 2.0, cex.lab = 2.1, cex.main = 2.0, lwd = 2)

mtext(text = "Frequency", side = 2, cex = 1.9, line = 3.5)

d <- density(list.sp$area_threat_5_perc, adjust = 1)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*1200, col="black", lwd=2) # add a density estimate with defaults
abline(v = 5, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.sp$area_threat_5_perc), 1), nsmall = 1)

nb_over_5 <- sum(list.sp$area_threat_5_perc > 5)
perc_over_5 <- format(round(nb_over_5/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_25 <- sum(list.sp$area_threat_5_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.sp$area_threat_5_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)


nb_0 <- sum(list.sp$area_threat_5_perc == 0)
perc_0 <- format(round(nb_0/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.sp$area_threat_5_perc == 100)
perc_100 <- format(round(nb_100/length(list.sp$area_threat_5_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " sp. (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " sp. (", perc_100,"%)"), 
                  paste0("> 5% = ", nb_over_5, " sp. (", perc_over_5,"%)"),
                  paste0("> 25% = ", nb_over_25, " sp. (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " sp. (", perc_over_50,"%)")), 
       x = "topright", inset = c(0, 0.03), y.intersp = 1.05, 
       cex = 1.8, bty = "n")

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 2.3, inset=c(0.01, 0.04))

# C/ For rings. Top 25% as threatened

hist(x = list.rings$area_threat_25_perc, 
     breaks = 20,
     freq = TRUE, col = "gray",
     xlim = c(0, 100),
     # ylim = c(0, 400),
     # main = "Proportion of mimicry ring distribution range\nunder threat (Top 25%)",
     main = "Mimicry rings: Top 25% threat     ",
     xlab = "Threatened area (%)",
     ylab = "",
     axes = T,
     cex.axis = 2.0, cex.lab = 2.1, cex.main = 2.0, lwd = 2)

mtext(text = "Frequency", side = 2, cex = 1.9, line = 3.5)

d <- density(list.rings$area_threat_25_perc, adjust = 0.75)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*200, col="black", lwd=2) # add a density estimate with defaults
abline(v = 25, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.rings$area_threat_25_perc), 1), nsmall = 1)

nb_over_25 <- sum(list.rings$area_threat_25_perc > 25)
perc_over_25 <- format(round(nb_over_25/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_50 <- sum(list.rings$area_threat_25_perc > 50)
perc_over_50 <- format(round(nb_over_50/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)
nb_over_75 <- sum(list.rings$area_threat_25_perc > 75)
perc_over_75 <- format(round(nb_over_75/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.rings$area_threat_25_perc == 0)
perc_0 <- format(round(nb_0/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

nb_100 <- sum(list.rings$area_threat_25_perc == 100)
perc_100 <- format(round(nb_100/length(list.rings$area_threat_25_perc)*100, 1), nsmall = 1)

legend(legend = c(paste0("0% threat = ", nb_0, " ring (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("100% threat = ", nb_100, " ring (", perc_100,"%)"), 
                  paste0("> 25% = ", nb_over_25, " rings (", perc_over_25,"%)"),
                  paste0("> 50% = ", nb_over_50, " rings (", perc_over_50,"%)"),
                  paste0("> 75% = ", nb_over_75, " ring (", perc_over_75,"%)")), 
       x = "topright", inset = c(0, 0.03), y.intersp = 1.05, 
       cex = 1.8, bty = "n")

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 2.3, inset=c(0.01, 0.03))

# D/ For rings. Top 5% as threatened

hist(x = list.rings$area_threat_5_perc, 
     breaks = seq(0, 20, 1),
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     # main = "Proportion of mimicry ring distribution range\nunder threat (Top 5%)",
     main = "Mimicry rings: Top 5% threat     ",
     xlab = "Threatened area (%)",
     ylab = "",
     axes = T,
     cex.axis = 2.0, cex.lab = 2.1, cex.main = 2.0, lwd = 2)

mtext(text = "Frequency", side = 2, cex = 1.9, line = 3.5)

d <- density(list.rings$area_threat_5_perc, adjust = 1)
d_x_cropped <- d$x[d$x > 0 & d$x < 100]
d_y_cropped <- d$y[d$x > 0 & d$x < 100]

lines(x = d_x_cropped, y = d_y_cropped*50, col="black", lwd=2) # add a density estimate with defaults
abline(v = 5, col = "red", lwd = 2, lty = 2) # Add vertical line for mean value

mean_threat_lvl <- format(round(mean(list.rings$area_threat_5_perc), 1), nsmall = 1)

nb_over_5 <- sum(list.rings$area_threat_5_perc > 5)
perc_over_5 <- format(round(nb_over_5/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)
nb_over_15 <- sum(list.rings$area_threat_5_perc > 15)
perc_over_15 <- format(round(nb_over_15/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)

nb_0 <- sum(list.rings$area_threat_5_perc == 0)
perc_0 <- format(round(nb_0/length(list.rings$area_threat_5_perc)*100, 1), nsmall = 1)


legend(legend = c(paste0("0% threat = ", nb_0, " rings (", perc_0,"%)"), 
                  paste0("Mean = ", mean_threat_lvl, "%"),
                  paste0("> 5% = ", nb_over_5, " rings (", perc_over_5,"%)"),
                  paste0("> 15% = ", nb_over_15, " rings (", perc_over_15,"%)")), 
       x = "topright", inset = c(0, 0.03), y.intersp = 1.05, 
       cex = 1.8, bty = "n")

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 2.3, inset=c(0.01, 0.12))

par(mar = original_int_margins, oma = original_ext_margins, mfrow = c(1, 1))

dev.off()


##### 5/ Define regions membership #####

# Then do with different colors for different regional types
# Lot of work, for not so much interest since results are quite obvious from other community maps looking at threat location

### 5.1/ Load all shp files
Guyana_Shield_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Guyana_Shield_shp.rds")
Caatinga_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Caatinga_shp.rds")
Cerrado_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Cerrado_shp.rds")
Mata_Atlantica_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp2.rds")
full_CA_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/full_CA_shp.rds")
Northern_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Northern_Andes_shp.rds")
Central_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_Andes_shp.rds")
Western_Lowlands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_shp.rds")
Coastal_desert_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Coastal_desert_shp.rds")
Chacos_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Chacos_shp.rds")
Llanos_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Llanos_shp.rds")
Caribbean_Islands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Caribbean_Islands_shp.rds")
Western_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Amazon_shp.rds")
Lower_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Lower_Amazon_shp.rds")

### 5.2/ Create masks by rasterizing
rasterize(x = sp_obj, y = raster_template, # Provide the grid to fill with CRS, bbox and resolution
          field = c("variable", value), # How to fill non empty cells. With the value of a varibel in the df of the sp_obj, or directly with a fixed value ?
          background = NA, # Value to use to fill empty cells
          fun = c("last", "count", sum, mean, custom_function, ...)) # What to do when multiple sp features overlap with a cell ?


