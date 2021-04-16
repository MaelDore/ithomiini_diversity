##### Filtering out Environmental outliers

# Effacer l'environnement
rm(list = ls())

load(file = "./Ithomiini_DB_final.RData")
load(file = "./list.units.RData")

# Besoin d'un example
# load(file ="testset_for_Mahalanobis.RData")
# Env.table <- testset_mahalanobis

# Real data

load(file = "Env2.5.RData")

Full_env <- env2.5

#### Merger les dataset + variables DB_orig = c(Keith, Jim)

# ACP initiale pour le 3D plot env.
library(ade4)

data_full <- data.frame(rep(NA,10000)) # Need to subsample for only 10000 pixels

index.sample <- sample(x = which(!is.na(Full_env[[1]]@data@values)), size = 10000) # Get index of cells which are not NA

for (i in 1:(length(names(Full_env)))) {
  data_full[,i] <- Full_env[[i]]@data@values[index.sample] # Retrieve 10000 values from each layer
}
names(data_full) <- names(Full_env)[1:length(names(Full_env))]

data_full <- data_full[,-which(names(data_full) == "Elevation")] # Remove Elevation

ACP <-  dudi.pca(df = data_full, center = T, nf = 2, scannf = F) # G?n?re l'ACP
# round(ACP$eig/sum(ACP$eig)*100,3) # % variances expliqu?e par les axes : 1e = 62,8% ; 2e = 24,9% ; 3e = 6.0%
# cumsum(round(ACP$eig/sum(ACP$eig)*100,3)) # Variance expliqu?e cumulative : 1e = 62,8%, 1+2e = 87,7% ; 1+2+3e = 93,7%
# ACP$co # Correlations variables ~ PCAxis

vec <- ACP$c1[,1:2] # Extract des Eigenvectors pour pouvoir projeter des nouvelles occurences sur les 2 premiers axes

# Tableau final pour stocker les infos d'outliers
names(Full_env)

Outliers.index <- data.frame(Tag=character(), ID_obs=character(), bio1=numeric(), bio3=numeric(), bio4=numeric(), bio12=numeric(), bio15=numeric(), Elevation=numeric(), Forests=numeric())
Ithomiini_final$outlier_maha <- F

# Librairies and parameters for 3D plots
library(rgl)
zoom <-  1
userMatrix <- matrix(data = c(0.8154871, 0.5784692, 0.01883241, 0, -0.1358215, 0.1596396, 0.97778738, 0, 0.5626130, -0.7999308, 0.20875259, 0, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow =T)
windowRect <- c(34,  57, 533, 551)
viewport <-  c(0, 0, 499, 494) ; names(viewport) <- c("x","y", "width", "height")

# Function to check for variance in dataset
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

# Run per unit

# j <- 265
# j <- 312

# for (j in 1:nrow(list.units)) {
for (j in 1:nrow(list.units)) {
  Tag <- list.units$Tag[j]
  sp <- list.units$Sp_full[j]
  
  Dataset <- Ithomiini_final[Ithomiini_final$Tag == Tag,] # Extraction des occurrences par unit
  
  if (nrow(Dataset)>2) { # Mahalanobis distance cannot be computed for 2 points
    Index.Keith <- Dataset$DB_orig == "Keith"
    Index.Jim <- Dataset$DB_orig == "Jim"
    
    Env.table <- extract(x = Full_env, y = Dataset[,c("Longitude","Latitude")]) # Extraction des valeurs 
    index.error <- which(!complete.cases(Env.table)) # Extraction des occurrences en erreur
    
    # R?cup?ration d'infos dans les pixels les plus proches de l'occurrence en erreur
    for (i in index.error) {
      Error.occ <- Env.table[i,]
      buffer <- 5000  # Premier buffer = 5km
      while (anyNA(Error.occ)) { # Elargit le buffer tant que l'on pas r?cup?r? des donn?es
        Error.occ <- unlist(extract(x = Full_env, y = Dataset[i,c("Longitude","Latitude")], buffer = buffer, fun = mean))
        buffer <- buffer + 5000 # Nouveau buffer = + 5km
      }
      cat("Buffer for",as.character(list.units$Tag[j]),"=",buffer/1000,"km\n") # Check la largeur finale du buffer
      Env.table[i,] <- Error.occ # Ecrit l'information
    }
    Env.table.orig <- Env.table # Stocke les valeurs non standardis?es
    
    # Standardisation des variables pour ?viter le probl?me d'?chelle lors du calcul de distance
    if(sum(apply(X = Env.table, MARGIN = 2, FUN = zero_range))==0) { # To check for the existence of variance in all variables
      Env.table <- as.data.frame(scale(Env.table, center = T, scale = T))
      
      if(det(as.matrix(cov(Env.table)))>2e-16) { # To avoid error with matrix that cannot be inverted
        # Calculate Mahalanobis Distance
        m_dist <- mahalanobis(Env.table, colMeans(Env.table), cov(Env.table))
        Env.table$m_dist <- round(m_dist, 2)
        
        # hist(Env.table$m_dist)
        # head(Env.table[order(Env.table$m_dist, decreasing = T),])
        
        # Mahalanobis Outliers for p = 7 = nb of env. variables
        p <- 7
        thres <- qchisq(p= 0.9999999, df = p-1) # Follow a Khi? distri for p-1 df (Etherington, 2019 => distri of p df !!!)
        Env.table$outlier_maha <- Env.table$m_dist > thres
        
        ### Tracer la loi de distri et les outliers sur le plot
        # x <- seq(0, max(Env.table$m_dist[Env.table$outlier_maha]*1.1), length=1000) ; y <- dchisq(x, df = p-1)
        # plot(x, y, type="l", lwd=1)
        # points(x = Env.table$m_dist[Env.table$outlier_maha], y = rep(0.01,sum(Env.table$outlier_maha)), pch = 16, cex = 1.2, col = "red")
        # segments(x0 = Env.table$m_dist[Env.table$outlier_maha], y0 = rep(0.01,sum(Env.table$outlier_maha)),
        #          x1 = Env.table$m_dist[Env.table$outlier_maha], y1 = rep(0,sum(Env.table$outlier_maha)), col = "red", lwd = 2)
        # abline(v=thres, col = "red", lty = 2)
        
        ### Liste des outliers
        if (sum(Env.table$outlier_maha) > 0) { # Si il y a des outliers
          if (sum(Env.table$outlier_maha) == 1) {
            New.outliers <- as.data.frame(t(c(as.character(Dataset$Tag[Env.table$outlier_maha]),Dataset$ID_obs[Env.table$outlier_maha], Env.table.orig[Env.table$outlier_maha,]))) # G?n?rer le tableau des outliers
          }else{
            New.outliers <- as.data.frame(cbind(as.character(Dataset$Tag[Env.table$outlier_maha]),Dataset$ID_obs[Env.table$outlier_maha], Env.table.orig[Env.table$outlier_maha,])) # G?n?rer le tableau des outliers
          }
          names(New.outliers) <- names(Outliers.index)
          Outliers.index <- rbind(Outliers.index, New.outliers) # Ajouter au tableau existant
          Ithomiini_final[Ithomiini_final$ID_obs %in% New.outliers[,2],"outlier_maha"] <- T # Stocker l'info dans le dataset
        }
        
        ### 3D plot only if outliers !
        
        if (sum(Env.table$outlier_maha) > 0) {
          coords <- as.matrix(Env.table[,names(data_full)])%*%as.matrix(vec)
          coords <- as.data.frame(cbind(coords,Env.table$Elevation)) ; names(coords) <- c("PC1", "PC2", "Elevation")
          open3d(zoom = 1, userMatrix = userMatrix, windowRect=windowRect, viewport=viewport) # Ouvre la fen?tre 3D avec une orientation d?finie
          plot3d(x = coords, type = "n") # To generate the axis
          # Plot des outliers de la base de Keith
          text3d(x = coords[(Env.table$outlier_maha)&(Index.Keith),], texts = Dataset$ID_obs[(Env.table$outlier_maha)&(Index.Keith)], col = "red", alpha = 1) # To add the labels in the desired color in the 3D plot
          points3d(x = coords[Index.Keith,], col=c("black","red")[Env.table[Index.Keith,]$outlier_maha+1]) # To add points in the desired color in the 3D plot
          # Plot des outliers de la base de Jim
          if (sum(Index.Jim)>0) {
            text3d(x = coords[(Env.table$outlier_maha)&(Index.Jim),], texts = Dataset$ID_obs[(Env.table$outlier_maha)&(Index.Jim)], col = "orange") # To add the labels in the desired color in the 3D plot
            points3d(x = coords[Index.Jim,], col=c("limegreen","orange")[Env.table[Index.Jim,]$outlier_maha+1]) # To add points in the desired color in the 3D plot
          }
          rgl.snapshot(filename = paste0("./Maps/Outliers_env/3D_plot/",Tag,".png"), top =T)
          rgl.close()
          # Dans le dossier par esp?ces
          file.copy(from =  paste0("./Maps/Outliers_env/3D_plot/",Tag,".png"), to = paste0(getwd(),"/Maps/By_species/",sp,"/Outliers_env_3D_",Tag,".png"), overwrite = T)
        }
      }
    }
  }
  if (j%%10 == 0)  { # Print incr?ment seulement tous les 10
    print(j)
  }
}

### Plot outliers in geographical space, per species

load(file = "list.sp.RData")

# Altitude
load(file = "Elevation10.RData")
Map <- Elevation10

# Maps on full study area, by mimicry
for (i in 1:nrow(list.sp)) {
  sp <- list.sp[i,"Sp_full"]
  sp.set <- subset(Ithomiini_final, subset = Ithomiini_final$Sp_full==sp)[,c("Longitude","Latitude","Mimicry_full","outlier_maha","DB_orig")]
  sp.set.outliers <- subset(Ithomiini_final, subset = (Ithomiini_final$Sp_full==sp) & (Ithomiini_final$outlier_maha==T))[,c("Longitude","Latitude","Mimicry_full","DB_orig")]
  
  if (nrow(sp.set.outliers)>0) { # Seulement si il y a des outliers
    # Dans le dossier par types de maps
    jpeg(filename = paste0("./Maps/Outliers_env/By_mimic/Full_extent/",sp,".jpeg"), quality =100) 
    plot(Map, main = sp)
    points(sp.set, col = 1, pch = (-as.numeric(factor(sp.set$DB_orig))+3)+15, cex = 0.6)
    points(sp.set, col = as.numeric(factor(sp.set$Mimicry_full))+1, pch = (-as.numeric(factor(sp.set$DB_orig))+3)+15, cex = 0.3)
    points(sp.set.outliers, col = 1, pch = (-as.numeric(factor(sp.set.outliers$DB_orig))+3)+15, cex = 1.3)
    col.outliers <- as.numeric(factor(sp.set$Mimicry_full))+1
    points(sp.set.outliers, col = col.outliers[sp.set$outlier_maha], pch = (-as.numeric(factor(sp.set.outliers$DB_orig))+3)+15, cex = 0.9)
    legend(legend = levels(sp.set$Mimicry_full), pch = 16, pt.cex = 1.3, col = 2:(1+length(levels(sp.set$Mimicry_full))), x = "bottomleft", cex = 0.9, bty ="o")
    dev.off()
    
    # Dans le dossier par esp?ces
    file.copy(from =  paste0("./Maps/Outliers_env/By_mimic/Full_extent/",sp,".jpeg"), to = paste0(getwd(),"/Maps/By_species/",sp,"/Outliers_env_By_mimic_Full_extent_",sp,".jpeg"), overwrite = T)
  }
  
  if (i%%10 == 0)  { # Print incr?ment seulement tous les 10
    print(i)
  }
}

# Maps on full study area, by sub.species
for (i in 1:nrow(list.sp)){
  sp <- list.sp[i,"Sp_full"]
  sp.set <- subset(Ithomiini_final, subset = Ithomiini_final$Sp_full==sp)[,c("Longitude","Latitude","Sub.species","outlier_maha","DB_orig")]
  sp.set.outliers <- subset(Ithomiini_final, subset = (Ithomiini_final$Sp_full==sp) & (Ithomiini_final$outlier_maha==T))[,c("Longitude","Latitude","Sub.species","outlier_maha","DB_orig")]
  
  if (nrow(sp.set.outliers)>0) { # Seulement si il y a des outliers
    # Dans le dossier par types de maps
    jpeg(filename = paste0("./Maps/Outliers_env/By_sub.species/Full_extent/",sp,".jpeg"), quality =100) 
    plot(Map, main = sp)
    points(sp.set, col = 1, pch = (-as.numeric(factor(sp.set$DB_orig))+3)+15, cex = 0.6)
    points(sp.set, col = as.numeric(factor(sp.set$Sub.species))+1, pch = (-as.numeric(factor(sp.set$DB_orig))+3)+15, cex = 0.3)
    points(sp.set.outliers, col = 1, pch = (-as.numeric(factor(sp.set.outliers$DB_orig))+3)+15, cex = 1.3)
    col.outliers <- as.numeric(factor(sp.set$Sub.species))+1
    points(sp.set.outliers, col = col.outliers[sp.set$outlier_maha], pch = (-as.numeric(factor(sp.set.outliers$DB_orig))+3)+15, cex = 0.9)
    legend(legend = unique(sp.set$Sub.species), pch = 16, pt.cex = 1.3, col = 2:(1+length(unique(sp.set$Sub.species))), x = "bottomleft", cex = 0.9, bty ="o")
    dev.off()
    
    # Dans le dossier par esp?ces
    file.copy(from =  paste0("./Maps/Outliers_env/By_sub.species/Full_extent/",sp,".jpeg"), to = paste0(getwd(),"/Maps/By_species/",sp,"/Outliers_env_By_sub.species_Full_extent_",sp,".jpeg"), overwrite = T)
  }
  if (i%%10 == 0)  { # Print incr?ment seulement tous les 10
    print(i)
  }
}

# G?n?ration du tableau final

Outliers.final <- merge(x = Outliers.index, y = Ithomiini_final, by = "ID_obs")

save(Outliers.final, file = "./Outliers.final.RData")
write.csv2(Outliers.final, file = "./Outliers.final.csv")

list.tag.outliers_env <- as.character(unique(Outliers.final$Tag.x))
save(list.tag.outliers_env, file = "./list.tag.outliers_env.RData")

