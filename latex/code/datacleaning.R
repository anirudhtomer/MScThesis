blood_donor=read.table(file="F:/docs/Dropbox/MSc Stats/Thesis/MScThesis/latex/code/dataset.txt", header = T)

blood_donor=blood_donor[,-c(8,9)]
blood_donor$Season<-ifelse(blood_donor$Season==3 | blood_donor$Season==0,"Cold","Hot")
blood_donor$Donate<-ifelse(blood_donor$Donate==0, "NO","YES")
#Age, Id, (Donation date, TSFV), Season, Donate are all available

#Hb is response but missing in 370 places, but there donation is also NO and 0 is the volume of the blood taken. This is like a rejection due to some other reason...not low Hb
#Every first observation has no blood donation, only 10ml is taken
#Volume is NA at other places, but that's coz no donation was done, yet we know HB. So these are true rejections
#We are anyway not interested in volume right now

nsubjects = length(unique(blood_donor$Id))
#checking why Hb is NA
hbNaBloodDonor = blood_donor[is.na(blood_donor$Hb),]
nrow(hbNaBloodDonor)
length(unique(hbNaBloodDonor$Id))
#the last two are not equal so some people didn't donate twice
lastVisitMap=tapply(blood_donor$Visit, blood_donor$Id, max)
for(i in hbNaBloodDonor$Id){
  
}
