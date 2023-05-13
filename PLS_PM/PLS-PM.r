#install.packages('devtools')
#devtools::install_github('gastonstat/plspm')

library(plspm)

dat <- read.delim('data.txt', sep = '\t',check.names = F)

dat_blocks <- list(
DA = c('DA'), 
AOM = c('TOC', 'TN'), 
pH= c('pH'), 
Ferrous= c('Ferrous'), 
Denitrification=c('Denitrification'),
 Nitrification=c('Nitrification'), 
 Ammonification=c('Ammonification'),
 Anammox=c('hzs'),
 DNRA=c('DNRA'), 
 nitrogen_fixation=c('Nitrogen fixation'))
dat_blocks

 DA<- c(0, 0, 0, 0, 0, 0, 0, 0, 0,0)
 AOM <- c(1, 0, 0, 0, 0, 0, 0, 0, 0,0)
 pH<- c(1, 1, 0, 0, 0, 0, 0, 0, 0,0)
 Ferrous<- c(1, 1, 1, 0, 0, 0, 0, 0, 0,0)
Denitrification<- c(1, 1, 1, 1, 0, 0, 0, 0, 0,0)
Nitrification<- c(1, 1, 0, 0, 0, 0, 0, 0, 0,0)
Ammonification<- c(1, 1, 0, 0, 0, 0, 0, 0, 0,0)
 Anammox<- c(1, 0, 1, 0, 0, 1, 0, 0, 0,0)
 DNRA<- c(1, 0, 0, 0, 0, 0, 0,1, 0,0)
 nitrogen_fixation<- c(1, 0, 0, 0, 0, 0, 1,0, 1,0)

dat_path <- rbind(DA, AOM, pH,Ferrous, Denitrification, Nitrification,Ammonification, Anammox,DNRA,   nitrogen_fixation )
colnames(dat_path) <- rownames(dat_path)
dat_path

dat_modes <- rep('A', 10)
dat_modes

dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes, scaled = T)
dat_pls
summary(dat_pls)
dat_pls$path_coefs
dat_pls$inner_model

innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

dat_pls$inner_summary

dat_pls$effects

dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

dat_pls$gof

dat_pls$scores

write.table(dat_pls$path_coefs,'Path_coefs.tsv',sep = '\t')
write.table(dat_pls$inner_summary,'inner_summary.tsv',sep = '\t')
write.table(dat_pls$effects,'effects.tsv',sep = '\t')

