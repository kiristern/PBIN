# abund_table<-read.table("ASVs_counts_copy.tsv",row.names=1,check.names=FALSE,header=T)
# 
# #Transform the table
# abund_table3 <-decostand(abund_table, method="hellinger")
# 
# meta_table1<-read.csv("metadata_w_cmd.csv",header=T,row.names=1,check.names=FALSE)
# 
# meta_table2=meta_table1
# colnames(meta_table2)
# str(meta_table1)
# meta_table2[,c(7:12)]<-decostand(meta_table1[,c(7:12)], method="standardize", na.rm = TRUE)
# 
# names(meta_table2)[names(meta_table2) == "Dolicho Abundance"] <- "Doli.Abundance"
# names(meta_table2)[names(meta_table2) == "Micro Abundance"] <- "Micro.Abundance"
# names(meta_table2)[names(meta_table2) == "Cyano Abundance"] <- "Cyan.Abundance"

sp.asv.rm
env_keep

#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(sp.asv.rm ~ Dolicho.Abundance + Micro.Abundance + Cyano.Abundance +
                               Mean_temperature_t0_t7+Total_Phosphorus_ug+Total_Nitrogen_mg+Dissolved_P+Dissolved_N+Cumulative_precipitation_t1_t7_mm, 
                               data=env_keep, permutations = 9999)
abund_table.adonis
#Extract the best variables
(bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.05])
#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- rda(sp.asv.rm ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=env_keep)",sep="")))
#You can use the following to use all the environmental variables
#sol<-cca(abund_table ~ ., data=meta_table2)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"),scaling=2)
#Extract site data first
df_sites<-data.frame(scrs$sites,bloom=as.factor(env_keep[,5]),Site=as.factor(env_keep[,3]),Months=as.factor(env_keep[,1]))
colnames(df_sites)<-c("x","y","bloom","Site","Months")
#Draw sites
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=bloom),size=4)
p<-p+labs(x="RDA1(%)",y="RDA2(%)")
p
#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines
# which are the coordinates of the heads of unit length vectors. In plot these are
# scaled by their correlation (square root of the column r2) so that "weak" predictors
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths
# using command scores. The plotted (and scaled) arrows are further adjusted to the
# current graph using a constant multiplier: this will keep the relative r2-scaled
# lengths of the arrows but tries to fill the current plot. You can see the multiplier
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul.
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.6)
p2<-p+geom_text(data=as.data.frame(df_arrows*1.3),aes(x, y, label = rownames(df_arrows),size=6,face="bold"),color="black",alpha=0.6, size=6, fontface="bold")
p2
# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")
# Either choose text or points
#p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species)))
#p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2)
p2<-p2+theme_bw()+geom_text(size=6, fontface="bold")
p3<-p2+scale_colour_manual(values = c("red","blue", "green")) + theme(panel.background = element_rect(fill = "white"))
p4<-p3+theme(axis.text.x  = element_text(vjust=0.5, size=12), axis.text.y  = element_text(vjust=0.5, size=12), axis.title.x = element_text(size = 15, face = "bold", color="black"),axis.title.y = element_text(size=15, face = "bold",color="black"),panel.border = element_rect(colour = "black", fill=NA, size=1))
p4
RsquareAdj (sol)
anova(sol)
