#########################################################################################
#                         -- KWMinima.R - 10.08.2017 --                                 #
# GIRDI:: *.tbl, *.dat dosyalari                                                        #
# ISLEM:: Kwee - van Woerden minima detection, 1956                                     #
# CIZIM:: 1) sT fonksiyonunu, polinom fiti ile birlikte cizer                           #
#         2) Minimum eðrisini çizer                                                     #
# CIKTI:: *.MIN: Tek bir dosyada min. verileri                                          #
# KAYNAK:: https://www.r-bloggers.com/fitting-polynomial-regression-in-r/               #
#          https://rpubs.com/kikihatzistavrou/80124                                     #
#          http://www.talkstats.com/threads/error-in-parabolic-regression-of-data-points.14476/
#          https://stackoverflow.com/questions/29999900/poly-in-lm-difference-between-raw-vs-orthogonal
#          https://courses.lumenlearning.com/boundless-algebra/chapter/graphs-of-quadratic-functions/
#                                                                                       #
#########################################################################################

### --BUYUK HARFLERLE-- YAZILAN BASLIKLAR ICINDEKI GIRDILER DUZENLENMELIDIR. Kucuk Harflerle yazilan basliklar altinda islemler yapilir. ###

####################
### Kutuphaneler ###
library(stats)
library(quantmod)

###################
### --AYARLAR-- ###
setwd("D:/Academia/Araþtýrma Projeleri/2019 - BAP (Ortak)/Veri/1_KOI202 (K412)/Kepler412SC/Extracted Data")
# D:/Academia/Araþtýrma Projeleri/2019 - BAP (Ortak)/Veri/KOI22 (K422)/Kepler422LC/MIN-AVE
# D:/Academia/Araþtýrma Projeleri/2019 - BAP (Ortak)/Veri/KOI192 (K427)/Kepler427SC//Extracted Data
# Girdi Sutunlari: [1]ObsTime [2]ObsFlux/Mag 
dosyalar <- list.files(path = ".", pattern = "*.tbl", all.files = FALSE)
ciktidosyasi <- "0_output_KW.MIN"
ciktidizini <- "../MIN-KW/"

Intp <- F # branch'larý ayri ayri interpole et
yuvarla <- 5 # JD kesirli yuvarlama hanesi
minnoktasayisi <- 16 # Minimumdaki nokta sayisi; LC: 12(+4), SC: 
verinoktasayisisiniri <- 5  # veride yer almasi gereken en az nokta sayisi #LC icin en az 5 #SC icin en az 20
### --AYARLAR-- ###
###################

####################
### Fonksiyonlar ###
# Constructing Quadratic Formula
QuadRoot <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
        x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
        x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
        QuadRoot = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
        x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
      b^2-4*a*c
} 
### Fonksiyonlar ###
####################

####################
### --ISLEMLER-- ###  
for(i in 1:length(dosyalar)){
if(file.exists(dosyalar[i])){
veri <- read.table(dosyalar[i], header=FALSE, skip=0, sep = "\t", stringsAsFactors=FALSE)

Tbas_tam <- as.integer(veri[1,1])
veri[,1] <- veri[,1] - as.integer(Tbas_tam)

N <- nrow(veri)
if(N > verinoktasayisisiniri){

# Her nokta arasi zaman farki bulunuyor
# ...
##############
# min Kontrolu
# if(veri[which.min(veri[,2]),1] < T1){ T1 <- veri[which.min(veri[,2]),1] }
# min Kontrolu
##############

#sN <- (T1-Ti)/(dT/2)

const_Z <- 1/4*N # In the case of linear interp.: Z =1/4*N, mag. were already equdistant in time, making interp. unnecesary, Z=1/2*N (Sec.2)
Nn <- 2*(2*const_Z)+1
#Nn <- 3600 # AVE'nin tahmini degeri: 3640
const_Zn <- 1/4*Nn

mag_b1 <- c(); mag_b2 <- c(); dmag <- c()
# verideki noktalar, Nn kadar interpole ediliyor.
veri <- as.data.frame(approx(veri, method="linear", n=Nn))

# Her nokta arasi zaman farki bulunuyor
dT <- c()
for(n in 1:Nn-1){
  dT[n] <- veri[n+1,1] - veri[n,1]
}
dT <- mean(dT)
                         
Ti <- veri[1,1] 
Ts <- veri[Nn,1]
T1 <- Ti + (Ts - Ti)/2

############
# sT Dongusu
repeat{
branch1 <- veri[round(veri[,1],yuvarla) <= round(T1,yuvarla), ]
branch2 <- veri[round(veri[,1],yuvarla) >= round(T1,yuvarla), ]

if(Intp){ n = 5000 # AVE'nin tahmini degeri: 3640
branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
}

# Min. Egrisi Kirpma ---
# Nn ile belirlenen kaba T1 degeri uzerinden min. egrisini esit iki parcaya ayirmak icin kirpma uygulaniyor (v2 ekle: elle T1 belirleme icin grafik cizdirilecek) 
if(nrow(branch1) != nrow(branch2)){
  
  if(nrow(branch1) > nrow(branch2)){
    n <- nrow(branch2)
    b1 <- branch1[order(-branch1[,1]), ]
    b1 <- as.data.frame(lapply(b1, "[", c(1:n)))
    b2 <- branch2
  }
  
  if(nrow(branch2) > nrow(branch1)){
    n <- nrow(branch1)
    b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
    b2 <- b2[order(-b2[,1]), ]
    b1 <- branch1  
  }                                    
   
} else {
  
  n <- nrow(branch1)
  b1 <- branch1
  b2 <- branch2[order(-branch2[,1]), ]

} 
# --- Min. Egrisi Kirpma

  for(j in 1:n){
    mag_b1[j] <- b1[j,2]
    mag_b2[j] <- b2[j,2]
    dmag[j] = mag_b1[j] - mag_b2[j]
  }
  
  sT1 <- sum(dmag^2)

  #Shifting dT for test
  # T1 + 0.5dT  
  T1art <- T1 + 0.5*dT
  branch1 <- veri[round(veri[,1],yuvarla) <= round(T1art,yuvarla), ]
  branch2 <- veri[round(veri[,1],yuvarla) >= round(T1art,yuvarla), ]

  if(Intp){ 
  branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
  branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
  }

  if(nrow(branch1) != nrow(branch2)){
    
    if(nrow(branch1) > nrow(branch2)){
      n <- nrow(branch2)
      b1 <- branch1[order(-branch1[,1]), ]
      b1 <- as.data.frame(lapply(b1, "[", c(1:n)))
      b2 <- branch2
    }
    
    if(nrow(branch2) > nrow(branch1)){
      n <- nrow(branch1)
      b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
      b2 <- b2[order(-b2[,1]), ]
      b1 <- branch1  
    }                                    
      
  } else {
  
    n <- nrow(branch1)
    b1 <- branch1
    b2 <- branch2[order(-branch2[,1]), ]
  
  }
  
  for(j in 1:n){
    mag_b1[j] <- b1[j,2]
    mag_b2[j] <- b2[j,2]
    dmag[j] = mag_b1[j] - mag_b2[j]
  }    
  sT1art <- sum(dmag^2)
  
  if(sT1art < sT1){
    T1 <- T1art
  } else {
  
    # T1 - 0.5dT           
    T1eks <- T1 - 0.5*dT
    branch1 <- veri[round(veri[,1],yuvarla) <= round(T1eks,yuvarla), ]
    branch2 <- veri[round(veri[,1],yuvarla) >= round(T1eks,yuvarla), ] 

  if(Intp){ 
  branch1 <- as.data.frame(approx(branch1, method="linear", n=n))
  branch2 <- as.data.frame(approx(branch2, method="linear", n=n))
  }

  if(nrow(branch1) != nrow(branch2)){
    
    if(nrow(branch1) > nrow(branch2)){
      n <- nrow(branch2)
      b1 <- branch1[order(-branch1[,1]), ]
      b1 <- as.data.frame(lapply(branch1, "[", c(1:n)))
      b2 <- branch2
    }
    
    if(nrow(branch2) > nrow(branch1)){
      n <- nrow(branch1)
      b2 <- as.data.frame(lapply(branch2, "[", c(1:n)))
      b2 <- b2[order(-b2[,1]), ]
      b1 <- branch1  
    }                                    
      
  } else {
  
    n <- nrow(branch1)
    b1 <- branch1
    b2 <- branch2[order(-branch2[,1]), ]
  
  }

    for(j in 1:n){
      mag_b1[j] <- b1[j,2]
      mag_b2[j] <- b2[j,2]
      dmag[j] = mag_b1[j] - mag_b2[j]
    }                                     
    sT1eks <- sum(dmag^2)
  
    if(sT1eks < sT1){
      T1 <- T1eks
    } else {
    break
    }
  }
}
# sT Dongusu  
############

Time <- c(T1eks, T1, T1art)
sT <- c(sT1eks, sT1, sT1art)

##############
# Polinom Fiti
# y = ax2 + bx + c, The y-intercept is the constant term, c. In every polynomial, the y-intercept is the constant term because the constant term is the value of y when x = 0.

# NOT
# R polinom fit uygulamasinda asagidaki iki yaklasim ayni olmakla birlikte 
# degisenlerin 'correlated' olmasini dikkate alir. 
# sTpoly <- lm(sT ~ poly(Time, 2, raw=F)) # Meth1
# sTpoly <- lm(sT ~ Time + I(Time^2))     # Meth2
# Note that when the explanatory variables are strongly correlated,
# the individual confidence intervals will usually underestimate
# the uncertainty in the parameter estimates, as indicated by the confidence region
# The rules of error propagation are different for correlated and uncorrelated errors
# -----------------------------------------------------------------------------------

sTfun <- function(Time) { Time + I(Time^2) }
optimize(sTfun, interval=c(Ti,Ts), maximum=F)

# Polinom Fiti: tum egri [ uncorrelated variables (orthogonal polynomials) ]
sTpoly <- lm(sT ~ poly(Time, 2))
sTa <- summary(sTpoly)$coefficients[3,1]
sTb <- summary(sTpoly)$coefficients[2,1]
sTc <- summary(sTpoly)$coefficients[1,1]
delta(sTa,sTb,sTc)
sTpolyroot <- QuadRoot(sTa,sTb,sTc)

# Minimuma yakin secilen iki nokta arasinda interpolasyon icin 10^-yuvarla kadar nokta turetme
yeniTime <- round(seq(Ti, Ts, 10^-yuvarla),yuvarla)
yenidf <- data.frame(Time=yeniTime)
yenisT <- predict(sTpoly, newdata=yenidf)

grafikverisi <- data.frame(yeniTime,yenisT)
minima <- grafikverisi[findValleys(grafikverisi[,2], 0),]

# sT fonk. grafigi
plot(yeniTime, yenisT, type='l', col="grey", cex=0.5)#, ylim=c(-0.02,0.02)) #, xlim=c(1.04,1.06))
points(Time, sT, col="blue", pch=16)
lines(Time, predict(sTpoly), col="red")

# Polinom Fiti: minumuma yakin noktalarin egrisi [ correlated variables ]
sTpoly2 <- lm(yenisT ~ yeniTime + I(yeniTime^2))
sTa2 <- summary(sTpoly2)$coefficients[3,1]
sTb2 <- summary(sTpoly2)$coefficients[2,1]
sTc2 <- summary(sTpoly2)$coefficients[1,1]
#approx(x=Time, y=sT, xout=T0)$y
cf = coef(sTpoly2)

#y <- sTa2^2*T0 + sTb2*T0 + sTc2

delta(sTa2,sTb2,sTc2)
sTpoly2root <- QuadRoot(sTa2,sTb2,sTc2)
sTpoly2root

# Minimum time
T0 <- (-sTb2/(2*sTa2))

# Minima mean error over Z
Z <- const_Z

# Minima mean error: poly(Time, 2)
# Tum Egri uzerinden hata hesabi: Z=1/4*N(asil verideki nokta sayisi)
sigT0 = sqrt((4*sTa*sTc-(sTb^2))/(4*(sTa^2)*(Z-1)))
sigT02 = sigT0^2

# Minumum civari uretilen egrinin hata hesabi
sigT02_2 = sqrt(abs((4*sTa2*sTc2-(sTb2^2)))/(4*(sTa2^2)*(Z-1)))
sigT02_2 = sigT02_2^2

# Tum Egri uzerinden hata hesabi, Z=1/4*Nn(interpolasyon icin uret. nokta sayisi)
Zn <- const_Zn
sigT02_n <- abs((4*sTa*sTc-(sTb^2)))/(4*(sTa^2)*(Zn-1))

#Get our covariance matrix
v <- vcov(sTpoly2)
b <- coef(sTpoly2)
#use delta method to calculate variance
xminvar <- (1/(2*b[3]))^2*v[2,2] + (b[2]/(2*b[3]^2))^2*v[3,3] - (b[2]/(2*b[3]^3))*v[2,3]

# Minima mean error: Time + I(Time^2)
#sTpolyMeth2 <- lm(sT ~ Time + I(Time^2))
#sTaMeth2 <- summary(sTpolyMeth2)$coefficients[3,1]
#sTbMeth2 <- summary(sTpolyMeth2)$coefficients[2,1]
#sTcMeth2 <- summary(sTpolyMeth2)$coefficients[1,1]
#sigT02Meth2 = sqrt(abs(4*sTaMeth2*sTcMeth2-(sTbMeth2^2))/(4*(sTaMeth2^2)*(Z-1)))

# Polinom Fiti
##############

T0 #+Tbas_tam
minima
sigT02
n
Time
sT
sTpoly

###########################
##### --DOSYAYA YAZ-- #####
yaz <- paste0(dosyalar[i],"\t", T0+Tbas_tam,"\t", sigT02)
write(yaz, file=paste0(ciktidizini,ciktidosyasi), append=TRUE)
##### --DOSYAYA YAZ-- #####
###########################

##########################
### --GRAFIGI YAZDIR-- ###
# pdf(), png(), postscript(): saydamligi desteklemez: pdf'i eps'ye cevir.
yazdir <- 1
if(yazdir){
dev.off()
png(paste0(ciktidizini,strsplit(basename(dosyalar[i]), "\\.")[[1]][1],".png"))#, width = 40, height = 30)
# En altta dev.off komutunu kapatmali; bu nedenle 'yazdir' degiseni ile kontrol ediliyor
}
### --GRAFIGI YAZDIR-- ###
##########################

##########################
##### --GRAFIK CIZ-- #####
### Min verisi grafigi ###
#dev.new()
dev.set(dev.prev()) 
plot(veri[,1]+Tbas_tam,veri[,2]) #, ylim=c(0.990, 1.002))
#lines(approx(veri[,1], veri[,2], method="linear", n=n))
#lines(aradeger)
title(main = dosyalar[i])
lines(b1)
lines(b2)
abline(v=T0+Tbas_tam, col="red", lwd="2")
#lines(spline(branch1[,1], branch1[,2]), df=70, method = "natural") 
##### --GRAFIK CIZ-- #####
##########################
if(yazdir){
dev.off() }

} # if(N > verinoktasayisisiniri) SONU
else {
###########################
##### --DOSYAYA YAZ-- #####
yaz <- paste0(dosyalar[i],"\t", "Nokta sayisi yetersiz", "\t")
write(yaz, file=paste0(ciktidizini,ciktidosyasi), append=TRUE)
##### --DOSYAYA YAZ-- #####
###########################
}

}} # DONGU SONU: i
                                    
### --ISLEMLER-- ###  
####################

####################
### --KOD SONU-- ###
winDialog("ok", paste0("Kod Sonlandý!","\n\nBÝLGÝ: ",
  "\nÝþlenen Dosya Top.Sayýsý = ", length(dosyalar))) 
  #"\nAtýlan Satýr Top.Sayýsý = ", GecersizSayaci))
### --KOD SONU-- ###
####################


#########################################################################################
#                          KWMinima.R - 10.08.2017                                      #
#                         ** Geliþtirme Notlari **                                      #
#                                                                                       #
#  s#: ...                                                                              #
#########################################################################################