# Load map for testing
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))


# HEX code:	#00afbe	
# RGB code:	rgb(0, 174, 210, max = 255)

# Deep dark blue
#00021d	RGB code:	rgb(0, 2, 29)

# Dark blue
#013458	RGB code:	rgb(1, 52, 88)

# Mid blue
#05879f	RGB code:	rgb(5, 135, 159)

# Light blue
#00b0b9	RGB code:	rgb(0, 174, 210)

# Green-blueish
#00b38d	RGB code:	rgb(0, 179, 141)

# Green
#8dc63b	RGB code:	rgb(141, 198, 59)

# Yellow
#fddf00	RGB code:	rgb(253, 223, 0)

# Orange
#f67d21	RGB code:	rgb(246, 125, 33)

# Reddish-orange
#ee3a22	RGB code:	rgb(238, 58, 34)

# Red
#ae1418	RGB code:	rgb(174, 20, 24)

# Dark red
#86080c	RGB code:	rgb(134, 8, 12)

# Deep dark red
#560000	RGB code:	rgb(86, 0, 0)

?colorRampPalette

DDblue_to_Dblue <- colorRampPalette(colors = c("#00021d", "#013458"))
Dblue_to_midblue <- colorRampPalette(colors = c("#013458", "#05879f"))
midblue_to_Lightblue <- colorRampPalette(colors = c("#05879f", "#00AED2"))
Lightblue_to_Grnbl <- colorRampPalette(colors = c("#00AED2", "#00b38d"))
Grnbl_to_Grn <- colorRampPalette(colors = c("#00b38d", "#8dc63b"))
Grn_to_Yellow <- colorRampPalette(colors = c("#8dc63b", "#fddf00"))
Yellow_to_Org <- colorRampPalette(colors = c("#fddf00", "#f67d21"))
Org_to_Orgred <- colorRampPalette(colors = c("#f67d21", "#ee3a22"))
Orgred_to_Red <- colorRampPalette(colors = c("#ee3a22", "#ae1418"))
Red_to_Drkred <- colorRampPalette(colors = c("#ae1418", "#86080c"))
Drkred_to_DDRed <- colorRampPalette(colors = c("#86080c", "#560000"))

test <- DDblue_to_Dblue(10)
test <- Dblue_to_midblue(10)
test <- midblue_to_Lightblue(10)
test <- Lightblue_to_Grnbl(10)
test <- Grnbl_to_Grn(10)
test <- Grn_to_Yellow(10)
test <- Yellow_to_Org(10)
test <- Org_to_Orgred(10)
test <- Orgred_to_Red(10)
test <- Red_to_Drkred(10)
test <- Drkred_to_DDRed(10)

pal_bl_red_Mannion <- c(gplots::col2hex("grey93"), DDblue_to_Dblue(20)[9:10], Dblue_to_midblue(12), midblue_to_Lightblue(12), Lightblue_to_Grnbl(18),
                        Grnbl_to_Grn(22), Grn_to_Yellow(24), Yellow_to_Org(23), Org_to_Orgred(22), Orgred_to_Red(22),
                        Red_to_Drkred(22), Drkred_to_DDRed(20))

plot(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion)
plot(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, colNA = "#c3edfd")
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)

plot(tot.sp.richness_Jaccard.80, col = pal_bl_red)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)

# Save the new palette
pal_bl_red_Mannion <- saveRDS(pal_bl_red_Mannion, file = "./maps/pal_bl_red_Mannion.rds", version = "2")
