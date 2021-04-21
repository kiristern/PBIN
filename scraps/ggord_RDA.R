#https://github.com/fawda123/ggord

library(ggord)

#run initialize script first
rem_rare_vir <- as.data.frame(t(filt_vir))
rem_rare_vir$year <-row.names(rem_rare_vir)

#keep date only (ie. remove everything before first period)
#change _ to .
rem_rare_vir$year <- gsub("_", ".", rem_rare_vir$year)
#keep only between 3rd and 4th period (just for year)
rem_rare_vir$year <- gsub("(?:[^.]+\\.){3}([^.]+).*", "\\1", rem_rare_vir$year)
  #(?:[^.]+\\.) is a group which matches non-period characters and then a single period. 
  #The {3} after the group means that the preceding token (the group) is repeated thrice - that is, 
  #"non-periods, followed by a period, followed by non-periods, followed by a period.". 
  #Then, the final ([^.]+) matches as many non-period characters as it can past the third period, 
  #thereby matching non-periods between the second period and the third period (or the end of the string).

#### PCA ###
dim(rem_rare_vir)
colnames(rem_rare_vir)

ord_vir <- prcomp(rem_rare_vir[,1:845])

p_vir <- ggord(ord_vir, rem_rare_vir$year, poly = FALSE) #transparent ellipses
p_vir

#embelish plot
p_vir + scale_shape_manual('Groups', values = c(1, 2, 3,4,5,6,7,8,9,10)) +
  theme_classic() 
  #theme(legend.position = "top") +
  

# change vector scaling, arrow length, line color, size, and type
p_vir <- ggord(ord_vir, grp_in = rem_rare_vir$year, arrow = 1, vec_ext = 3, veccol = 'red', veclsz = 0.5, vectyp = 'dotted')
p_vir
  
# faceted by group
p_vir <- ggord(ord_vir, rem_rare_vir$year, facet = TRUE, nfac = 1)
p_vir


# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(rem_rare_vir[,1:845])
ggord(ord, rem_rare_vir$year)


# distance-based redundancy analysis
#need to run dbRDA.R script
head(environnement)
ord <- dbrda(species ~ Months + Years + Site + Period + bloom2 +
                Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp+
               Cyano.Abundance + Micro.Abundance + Dolicho.Abundance, 
              data=complete_env_keep, dist="bray")
ggord(ord)

