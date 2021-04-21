vir_abun_removed
complete_env_keep

species <- vir_abun_removed
environnement <- complete_env_keep

#decide which distance measure to use by looking at the rank correlations between dissimilarity indices and gradient separation (the higher the value the better)
rankindex(environnement, species, indices = c("euc", "man", "gow", "bra", "kul"), stepacross = F, method = "spearman")
#The Kulczynski measure seems like a good distance measure to use. We’ll proceed with the Bray-Curtis distance measure for this example though since it is one we are familiar with and the difference between bra and kul is marginal.

dbRDA <- capscale(species ~ Months + Years + Site + Period + bloom2 +
                    Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp+ 
                    Cyano.Abundance + Micro.Abundance + Dolicho.Abundance
                    , data=complete_env_keep, dist="bray")
plot(dbRDA)
anova(dbRDA) #overall test of model significance
anova(dbRDA, by="axis", perm.max=500) #test axes for significance
anova(dbRDA, by="terms", permu=200) #test for sig env vars





#Transforming negative eigenvalues by
#1) adding a constant
dbRDA_add <- capscale(species ~Months + Years + Site + Period + bloom2 +
                        Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp, data=complete_env_keep, dist="bray", 
                      add = T)
plot(dbRDA_add)
anova(dbRDA_add)
#2) take the sqrt of dissimilarities
dbRDA_sqrt <- capscale(species ~ Months + Years + Site + Period + bloom2 +
                         Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp+cyano_count, data=complete_env_keep, dist="bray",
                       sqrt.dist=T)
plot(dbRDA_sqrt)
anova(dbRDA_sqrt)
#3) Sqrt transformation, Wisconsin db standartdization (this emmplahizes the env vars)
dbRDA_medaMDF <- capscale(species ~Months + Years + Site + Period + bloom2 +
                            Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp+cyano_count, data=complete_env_keep, dist="bray",
                          metaMDF=T, sqrt.dist = T)
plot(dbRDA_medaMDF)
anova(dbRDA_medaMDF)

#Modify db-RDA plot (optional) to change the vector sizes or include only the important axes (since we know the significant env vars)
#getting the scores
scrs_dbRDA <- scores(dbRDA) #species and sites are together
site_scrs <- scrs_dbRDA$sites # separating out the site scoes, get CAP1 and CAP2 scores
fix(site_scrs)







