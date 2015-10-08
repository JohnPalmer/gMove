load(file='~/research/gHuman/aaag_analysis/placa_esp_uab_car_samples.RData')

s1 = sample_list[[2]][[1]]
s1 = spTransform(s1, center=T)
length(s1)
test = brownian.bridge.dyn(s1, raster=500,location.error=rep(10, length(s1)), ext=2, by.step=T)

plot(s1)
for(ti in test){
contour(ti, add=TRUE)
}

contour(test[[90]], add=TRUE)

sapply(1:nrow(test[[12]]), function(i) sum(test[[12]][i,]))

