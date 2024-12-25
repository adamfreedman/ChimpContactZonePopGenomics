if (file.exists("~/software/eems/plotting/rEEMSplots")) {
  install.packages("rEEMSplots",repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

coord=read.table("maf05_Phase_3_CHIMP_10k_NEUTRAL.coord")
mcmcpath =c("chimps_neutral_nDemes200-1","chimps_neutral_nDemes200-2","chimps_neutral_nDemes200-3","chimps_neutral_nDemes200-4","chimps_neutral_nDemes200-5")
plotpath = "plots_2024.08.23"

eems.plots(mcmcpath, plotpath, longlat = TRUE,projection.in='+proj=longlat +datum=WGS84',projection.out='+proj=longlat +datum=WGS84',add.map=TRUE,m.plot.xy={ points(coord, col = "black") },q.plot.xy={ points(coord, col = "black") })

