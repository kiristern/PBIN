viral <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
viral <- t(viral)
bacterial <- read.csv("metadata_w_cmd.csv", row.names=1, header=T)
bacterial <- bacterial[,c(13:15)]

#merge dolico with viral
doli <- merge(bacterial[1], viral, by="row.names", all=T)
#transform col1 into row.names
doli2 <- doli[,-1]
rownames(doli2)<-doli[,1]
#transform data
doli_trans <-decostand(doli2, method="hellinger", na.rm = T)

#merge micro with viral
mic <- merge(bacterial[2], viral, by="row.names", all=T)
#transform col1 into row.names
mic2 <- mic[,-1]
rownames(mic2)<-mic[,1]
#transform data
mic_trans <-decostand(mic2, method="hellinger", na.rm = T)
