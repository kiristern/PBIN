#Hellinger transformed and filtered: RDA + MRT, ML

bact_hel
vir_hel

### Cyanobacteria ###
bactps_helli <- transform(bact_physeq, transform = "hellinger", target = "OTU")
(cyano_ps_helli <- subset_taxa(bactps_helli, Phylum == "p__Cyanobacteria"))
(cyano_helli_filt_ps = filter_taxa(cyano_ps_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE))
cyano_helli_filt <- cyano_helli_filt_ps %>% otu_table() 


#### Viral ####
(virps_helli <- transform(viral_physeq, transform = "hellinger", target = "OTU"))
(virps_helli_filt = filter_taxa(virps_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE))
vir_helli_filt <- virps_helli_filt %>% otu_table()
vir_helli_filt <- samp2date(vir_helli_filt)




