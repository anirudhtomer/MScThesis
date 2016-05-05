blood_donor=read.table(file="dataset.txt", header = T)

blood_donor = blood_donor[,-c(8,9)]
blood_donor$Season = as.factor(ifelse(blood_donor$Season==3 | blood_donor$Season==0,"Cold","Hot"))
blood_donor$Donate = ifelse(blood_donor$Donate==0, FALSE,TRUE)
blood_donor$firstAge = rep(0, nrow(blood_donor))
blood_donor$donationLast2Years = rep(0, nrow(blood_donor))

#Age, Id, (Donation date, TSFV), Season, Donate are all available

#Hb is response but missing in 370 places, but there donation is also NO and 0 is the volume of the blood taken. This is like a rejection due to some other reason...not low Hb
#Every first observation has no blood donation, only 10ml is taken
#Volume is NA at other places, but that's coz no donation was done, yet we know HB. So these are true rejections
#We are anyway not interested in volume right now

#select all those subjects who don't have a NA in Hb
rejectIdList = c()
for(i in unique(blood_donor$Id)){
  if(sum(is.na(blood_donor[blood_donor$Id == i,]$Hb))>0){
    rejectIdList = c(rejectIdList, i) 
  }
}

ds = blood_donor[!(blood_donor$Id %in% rejectIdList),]

#remove rows with volume 10
#Process the selected blood donors
uniqueIds = unique(ds$Id)

randomSample = sample(1:1595, size = 250, replace = F)
uniqueIds = uniqueIds[randomSample]

ds = ds[ds$Id %in% uniqueIds, ]

timePerSubject = list()
agePerSubject = list()
donationLast2YearsPerSubject = list()
nsubjects = length(uniqueIds)
nrepPerSubject = rep(0, nsubjects)
tspdpersubject = list()

for(i in 1:nsubjects){
  id = uniqueIds[i]
  
  temp = ds[ds$Id==id,]
  timePerSubject[[paste(id)]] = temp$TSPD
  agePerSubject[[paste(id)]] = temp$Age
  ds[ds$Id==id,]$firstAge = rep(min(temp$Age), length(temp$Age))
  
  donationLast2years = rep(0, nrow(temp))
  for(j in 2:nrow(temp)){
    for(k in j:1){
      if(sum(temp$TSPV[j:k])>365*2){
        break
      }
    }
    if(k<j){
      donationLast2years[j] = sum(temp$Donate[(j-1):k])
    }
  }
  ds[ds$Id==id,]$donationLast2Years = donationLast2years/10
  
  donationLast2YearsPerSubject[[paste(id)]] = donationLast2years/10
  tspdpersubject[[paste(id)]] = temp$TSPD
}

ds$TSPDdonate = ds$TSPD*ds$Donate
ds$ageSeason = ds$Age*(as.numeric(ds$Season)-1)
ds$TSPDSeason = ds$TSPD*(as.numeric(ds$Season)-1)
ds$donateLast2TSPD = ds$TSPD*ds$donationLast2Years
ds$donateLast2Donate = ds$Donate*ds$donationLast2Years
ds$donateLast2Season = (as.numeric(ds$Season)-1)*ds$donationLast2Years
ds$donateLast2Square = ds$donationLast2Years^2


numHb = table(ds$Id)
cumsumHb = c(0, cumsum(numHb))

write.csv(ds[,-2], file = "dataset.csv", row.names = F)

rm(rejectIdList)
rm(blood_donor)
rm(i)
rm(j)
rm(k)
rm(id)
rm(donationLast2years)
rm(temp)