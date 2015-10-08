load(file='~/research/gHuman/aaag_analysis/placa_esp_uab_car_samples.RData')

s1 = sample_list[[3]][[1]]
s1 = spTransform(s1, center=T)
length(s1)
test = brownian.bridge.dyn(s1, raster=100,location.error=rep(10, length(s1)), ext=.5, by.step=T)

pdf('~/research/gHuman/aaag_analysis/placa_esp_uab_car_sample_x_stepwise_dBBMM.pdf')
plot(s1, ylim=c(-8000, 8000))
for(ti in test){
contour(ti, add=TRUE, levels=seq(0,.9, .1))
}
points(s1, col="red", pch=20, cex=2)
lines(s1, col="red", lwd=2)
dev.off()

