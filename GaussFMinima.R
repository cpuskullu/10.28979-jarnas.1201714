#########################################################################################
#                      -- GaussFMinima.R - 08.01.2020 --                                #
# GIRDI:: *.tbl, *.dat verisi (gecis, isik egrisi)                                      #
# ISLEM:: Gecis egrisi verisine Gauss fiti yapar. "output_KW" gerekli.                  #
# CIZIM:: Gecis ortasi                                                                  #
# CIKTI:: *.MIN: Tek bir dosyada min. verileri                                          #
# KAYNAK:: 
# 1) https://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian
# 2) https://stats.stackexchange.com/questions/220109/fit-a-gaussian-to-data-with-r-with-optim-and-nls
# 3) https://earlglynn.github.io/RNotes/package/minpack.lm/index.html
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR.      ###
### Kucuk Harflerle yazilan basliklar altinda islemler yapilir.                   ###

####################
### Kutuphaneler ###
library(minpack.lm)

###################
### --AYARLAR-- ###
setwd("D:/Academia/Arastirma Projeleri/2019 - BAP (Ortak)/Veri/5_KOI196 (Kepler41)/Kepler41SC/Extracted Data")
dosyalar <- list.files(path = ".", pattern = "*.tbl", all.files = FALSE)
output_KW <- read.table(file="../MIN-KW/0_output_KW.MIN", skip=0, header=FALSE, sep="\t", dec = ".", numerals=c("no.loss"), stringsAsFactors=FALSE)

ciktidosyasi <- "0_output_GF.MIN"
ciktidizini <- "../MIN-GF/"

yuvarla <- 6 # JD kesirli yuvarlama hanesi (LC verisi icin '5', SC verisi icin '7')
verinoktasayisisiniri <- 10  # veride yer almasi gereken en az nokta sayisi (min civari nokta sayisina bakilacak)

# Initial Values
A <- 0.0120 # alfa   (amplitude)
B <- 1.000 # beta    (normalization shift)
sigm <- 0.030  #     (width)

### --AYARLAR-- ###
###################

####################
### --ISLEMLER-- ###
#dosyalar <- "kepler412-002.ctbl"; i=2;
for(i in 1:length(dosyalar)){
# Veri Sutunlari: ObsTime ObsFlux/Mag
veri <- read.table(dosyalar[i], skip=0, header=FALSE, sep="\t", dec = ".", numerals=c("no.loss"), stringsAsFactors=FALSE)

N <- nrow(veri)
if(N > verinoktasayisisiniri){
# Determine Normalisation Value: ilk 3 ve son uc nokta kullaniliyor
normV <- mean(max(veri[,2][1:3]):max(veri[,2][nrow(veri)-3]:veri[,2][nrow(veri)]))

# Determine initial value: mu (T0)
Ti <- veri[1,1]
Ts <- veri[N,1]

# Read inital mu from output_KW
mu <- as.double(output_KW[,2][i])
#mu <- Ti + (Ts - Ti)/2

# Her nokta arasi zaman farki bulunuyor
dT <- c()
for(n in 1:N-1){
  dT[n] <- veri[n+1,1] - veri[n,1]
}
dT <- mean(dT)

# Add when v2: sT func. to find initial mu

# Veriler yeni degiskene aktariliyor
Time <- veri[,1]
g <- veri[,2]/normV


##########################
##### --GRAFIK CIZ-- #####
### -- Min grafigi --  ###
if(F){
plot(Time, g, type="l")
title(main = dosyalar[i])
#lines(Time, B-gf, col="red")
abline(v=mu, col="green", lwd="2")
}
##### --GRAFIK CIZ-- #####
##########################

##########################
### --GRAFIGI YAZDIR-- ###
# pdf(), png(), postscript(): saydamligi desteklemez: pdf'i eps'ye cevir.
yazdir <- 1
if(yazdir){
png(paste0(ciktidizini,strsplit(basename(dosyalar[i]), "\\.")[[1]][1],".png"))#, width = 40, height = 30)
# En altta dev.off komutunu kapatmali; bu nedenle 'yazdir' degiseni ile kontrol ediliyor
}
### --GRAFIGI YAZDIR-- ###
##########################

# Gauss fonksiyonu hesaplaniyor, gecis ortasi zamani 'mu' bulunuyor.
# YONTEM 1: NLS
nlinModGaussF <- nls(g ~ B-A*exp(-((Time-mu)^2)/(2*sigm^2)), start=c(mu=mu,sigm=sigm,A=A,B=B), data=veri, model=T, trace=T)
residuals(nlinModGaussF)
fitted.values(nlinModGaussF)
fitc <- paste0(dosyalar[i],"\t", residuals(nlinModGaussF),"\t", fitted.values(nlinModGaussF))
write(fitc, file=paste0(ciktidizini,ciktidosyasi), append=TRUE)

v <- summary(nlinModGaussF)$parameters[,"Estimate"] # coef(nlinModGaussF)
v_err <- summary(nlinModGaussF)$parameters[,"Std. Error"]

##########################
##### --GRAFIK CIZ-- #####
### GF sonrasi grafigi ###
plot(g~Time, data=veri, ylim=c(range(g)[1]-0.005,range(g)[2]))
title(main = dosyalar[i])
plot(function(Time) v[4]-v[3]*exp(-1/2*(Time-v[1])^2/v[2]^2), add=T, xlim=range(Time), col="blue", lwd="3")
abline(v=v[1], col="red", lwd="2")
##### --GRAFIK CIZ-- #####
##########################
# YONTEM 1 sonu

if(F){
# YONTEM 2: NLSLM (munipack.lm) * Fark: Y2-Y1: 10^-9
nlinModGaussFLM <- nlsLM(g ~ B-A*exp(-((Time-mu)^2)/(2*sigm^2)), start=c(mu=mu,sigm=sigm,A=A,B=B), data=veri, model=T, trace=T)
residuals(nlinModGaussF)
fitted.values(nlinModGaussF)
vLM <- summary(nlinModGaussFLM)$parameters[,"Estimate"] # coef(nlinModGaussF)
v_errLM <- summary(nlinModGaussFLM)$parameters[,"Std. Error"]

##########################
##### --GRAFIK CIZ-- #####
### GF sonrasi grafigi ###
if(F){
yenig <- v[4]-v[3]*exp(-1/2*(Time-v[1])^2/v[2]^2)
dev.new()
plot(Time, g, type="p")
title(main = dosyalar[i])
lines(Time, predict(fitted.values(nlinModGaussF)), col="red", lwd="1")
#lines(Time, yenig, col="red")
abline(v[1], col="blue", lwd="2")
}
##### --GRAFIK CIZ-- #####
##########################
} # YONTEM 2 sonu

###########################
##### --DOSYAYA YAZ-- #####
yaz <- paste0(dosyalar[i],"\t", v[1],"\t", v_err[1])
write(yaz, file=paste0(ciktidizini,ciktidosyasi), append=TRUE)
##### --DOSYAYA YAZ-- #####
###########################

if(yazdir){
dev.off() }

} # if(N > verinoktasayisisiniri) SONU
else {
###########################
##### --DOSYAYA YAZ-- #####
yaz <- paste0(dosyalar[i],"\t", "Nokta sayisi yetersiz","\t")
write(yaz, file=paste0(ciktidizini,ciktidosyasi), append=TRUE)
##### --DOSYAYA YAZ-- #####
###########################
}


} # DONGU SONU: i

### --ISLEMLER-- ###
####################

####################
### --KOD SONU-- ###
winDialog("ok", paste0("Kod Sonlandi!","\n\nBiLGi: ",
  "\ninlenen Dosya Top.Sayisi = ", length(dosyalar)))
  #"\nAt�lan Sat�r Top.Say�s� = ", GecersizSayaci))
### --KOD SONU-- ###
####################

#########################################################################################
#                        GaussFMinima.R - 05.01.2021                                    #
#                         ** Gelistirme Notlari **                                      #
#                                                                                       #
#  s#: ...                                                                              #
#########################################################################################




