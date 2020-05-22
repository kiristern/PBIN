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

ts <- melt(vir_ts, id="sample")
ggplot(ts) +
  geom_line(aes(x=Var1, y=value, group=Var2, color=Var2))

