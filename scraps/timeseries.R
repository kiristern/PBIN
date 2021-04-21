filt_vir
viral <- as.data.frame(t(filt_vir))

#store data in timeseries object
vir20 <- viral %>% select(
  "ASV_2",
  "ASV_1",
  "ASV_3",
  "ASV_5",
  "ASV_4",
  "ASV_7",
  "ASV_8",
  "ASV_12",
  "ASV_6",
  "ASV_21",
  "ASV_15",
  "ASV_33",
  "ASV_11",
  "ASV_18",
  "ASV_10",
  "ASV_17",
  "ASV_19",
  "ASV_48",
  "ASV_13",
  "ASV_25"
)

vir_ts <- ts(vir20)

vir_ts$sample <- row.names(vir_ts)

ts <- reshape2::melt(vir_ts, id="sample")
ggplot(ts) +
  geom_line(aes(x=Var1, y=value, group=Var2, color=Var2))

start(vir_ts)
end(vir_ts)
frequency(vir_ts)
summary(vir_ts)

decompose(vir_ts)

##### Timeseries by season ######
#combine ASV with meta
asv <- as.data.frame(t(filt_vir))
meta_cyano

row.names(asv) %in% row.names(meta_cyano)
df <- merge(asv, meta_cyano, by="row.names")
colnames(df)
df_keep <- df[, c(1:846, 851)]

df_keep$Period


