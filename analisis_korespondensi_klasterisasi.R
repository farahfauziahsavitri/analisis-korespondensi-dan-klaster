#Data
covid=read.csv("C:/Users/farah/Dropbox/PC/Downloads/C-Download/github/11_21.csv",sep=';')
covid
data=covid[,-8]
data=data[,-1:-2]
kecamatan=covid[,1]
rownames(data)=kecamatan
data

#Statistik deskriptif
summary(data)
positif=data$positif
suspek=data$suspek
kontak_erat=data$kontak_erat
meninggal=data$meninggal
sembuh=data$sembuh
sd(positif)
sd(suspek)
sd(kontak_erat)
sd(meninggal)
sd(sembuh)
data<-as.table(as.matrix(data))
data


#Uji chi square #mengetahui ada tidaknya korelasi antar variabel
test<- chisq.test(data)
test

#Analisis korespondensi
ankor=function(data){
  total=sum(data)
  totBaris=apply(data, 1,sum)
  totKolom=apply(data, 2,sum)
  p=data/total
  r=totBaris/total
  R=diag(r)
  c=totKolom/total
  C=diag(c)
  S=solve(R)^(0.5)%*%(p-t(t(r))%*%c)%*%solve(C)^(0.5)
  lamda=eigen(S%*%t(S))$values[1:(min(nrow(data),ncol(data))-
                                    1)]
  D=diag(sqrt(lamda))
  prop=c()
  propcum=c()
  for(i in (1:length(lamda))){
    prop=c(prop,lamda[i]/sum(lamda))
    propcum=c(propcum,sum(prop[1:i]))
  }
  if(nrow(data)<ncol(data)){
    U=(eigen(S%*%t(S))$vectors)[,1:(min(nrow(data),ncol(data))-
                                      1)]
    V=solve(D)%*%U%*%S
    Y=solve(R)^(0.5)%*%U%*%D
    Z=solve(C)^(0.5)%*%t(V)%*%D
  }else {
    V=(eigen(S%*%t(S))$vectors)[,1:(min(nrow(data),ncol(data))-
                                      1)]
    U=t(S)%*%V%*%solve(D)
    Y=solve(R)^(0.5)%*%V%*%D
    Z=solve(C)^(0.5)%*%U%*%D
  }
  hasil=list("Inersia"=lamda,"Proporsi
Inersia"=prop,"Proporsi Kumulatif Inersia"=propcum,"Estimasi
Koordinat Utama dari Baris"=Y,"Estimasi Koordinat Utama dari
Kolom"=Z,"Matriks_residual_standar"=S)
  library(ca)
  plot=plot(ca(data))
  return(c(hasil,plot))
}
P=ankor(data)
P

#Analisis Klaster


#Identifikasi Outlier
library(MVN)
par(mfrow=c(1,1))
mvn(X,mvnTest = "mardia",multivariateOutlierMethod = "quan")

X=P$rows[,1:2]
X
multicol=function(X){
  VIF=diag(solve(cor(X)))
  result=ifelse(VIF>10,"mulicolinearity", "non multicolinearity")
  data1=data.frame(VIF,result)
  return(data1)
}
multicol(X)

#matriks korelasi antar dimensi
library(car)
multikol<-cor(X)
multikol
#terjadinya korelasi negatif antar dimensi

#kecukupan sampel
library(psych)
kmo <- function(x)
{
  x <- subset(x, complete.cases(x))       # menghilangkan data kosong (NA)
  r <- cor(x)                             # Membuat matrix korelasi
  r2 <- r^2                               # nilai koefisien untuk r squared
  i <- solve(r)                           # Inverse matrix dari matrix korelasi
  d <- diag(i)                            # element diagonal dari inverse matrix
  p2 <- (-i/sqrt(outer(d, d)))^2          # koefisien korelasi Parsial kuadrat
  diag(r2) <- diag(p2) <- 0               # menghapus element diagonal 
  KMO <- sum(r2)/(sum(r2)+sum(p2))
  MSA <- colSums(r2)/(colSums(r2)+colSums(p2))
  return(list(KMO=KMO, MSA=MSA))
}
kmo(X) #sampel cukup untuk digunakan jika kmo>0.5

#Jarak Euclidean
D<-dist(X,method="euclidian")
D

#penentuan jml klaster
library(pvclust)
library(factoextra)
library(cluster)

gap_stat<-clusGap(X,FUN=hcut,K.max=10,B=1000)
fviz_gap_stat(gap_stat)


#----------------------------METODE AVERAGE------------------------------------
hierarkiave<-hclust(dist(data), method="ave")
hierarkiave
windows()
plot(hierarkiave) #dendogram

rect.hclust(hierarkiave,7) 		#plot mengelompokkan data
anggotaave<-cutree(hierarkiave,k=7) #hasil kelompok data
anggotaave

tabulasiave<-data.frame(X,anggotaave) # hasil kelompok data dalam bentuk data frame
View(tabulasiave)

cophenetic(hierarkiave) #jarak cophenetic average
#korelasi cophenetic
d1 <- dist(X)
hc <- hclust(d1, "ave")
d2 <- cophenetic(hc)
corave=cor(d1, d2)
corave

#-----------------------------METODE COMPLETE----------------------------------
hierarkicomp<-hclust(dist(scale(data)), method="complete")
hierarkicomp
windows()
plot(hierarkicomp) #dendogram

rect.hclust(hierarkicomp,7) 		  #plot mengelompokkan data
anggotacomp<-cutree(hierarkicomp,k=7) #hasil kelompok data
anggotacomp

tabulasicomp<-data.frame(data,anggotacomp) # hasil kelompok data dalam bentuk data frame
View(tabulasicomp)

cophenetic(hierarkicomp) #jarak cophenetic complete 
#korelasi cophenetic
d1 <- dist(data)
hc <- hclust(d1, "complete")
d2 <- cophenetic(hc)
corcomp=cor(d1, d2)
corcomp

#-------------------------------METODE SINGLE----------------------------------
hierarkising<-hclust(dist(scale(data)), method="single")
hierarkising
windows()
plot(hierarkising) #dendogram

rect.hclust(hierarkising,7) 		  #plot mengelompokkan data
anggotasing<-cutree(hierarkising,k=7) #hasil kelompok data
anggotasing

tabulasising<-data.frame(data,anggotasing) # hasil kelompok data dalam bentuk data frame
View(tabulasising)

cophenetic(hierarkising) #jarak cophenetic single
#korelasi cophenetic
d1 <- dist(data)
hc <- hclust(d1, "single")
d2 <- cophenetic(hc)
corsing=cor(d1, d2)
corsing

#---------------------------------METODE WARD----------------------------------
hierarkiward<-hclust(dist(scale(data)), method="ward.D")
hierarkiward
windows()
plot(hierarkiward) #dendogram

rect.hclust(hierarkiward,7) 		  #plot mengelompokkan data
anggotaward<-cutree(hierarkiward,k=7) #hasil kelompok data
anggotaward

tabulasiward<-data.frame(data,anggotaward) # hasil kelompok data dalam bentuk data frame
View(tabulasiward)

cophenetic(hierarkiward) #jarak cophenetic ward
#korelasi cophenetic
d1 <- dist(data)
hc <- hclust(d1, "ward.D")
d2 <- cophenetic(hc)
corward=cor(d1, d2)
corward

#-------------------------------METODE CENTROID--------------------------------
hierarkicent<-hclust(dist(scale(data)), method="centroid")
hierarkicent
windows()
plot(hierarkicent) #dendogram

rect.hclust(hierarkicent,7) 		  #plot mengelompokkan data
anggotacent<-cutree(hierarkicent,k=7) #hasil kelompok data 
anggotacent

tabulasicent<-data.frame(data,anggotacent) # hasil kelompok data dalam bentuk data frame
View(tabulasicent)

cophenetic(hierarkicent) #jarak cophenetic centroid
#korelasi cophenetic
d1 <- dist(data)
hc <- hclust(d1, "centroid")
d2 <- cophenetic(hc)
corcent=cor(d1, d2) 
corcent