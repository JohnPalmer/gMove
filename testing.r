load(file='~/gHuman/aaag_analysis/placa_esp_uab_car_samples.RData')

s1 = sample_list[[3]][[1]]
s1 = spTransform(s1, center=T)
test = brownian.bridge.dyn(s1, raster=5000,location.error=rep(10, length(s1)),ext=1.5, by.step=T)


sum(test[[13]])
