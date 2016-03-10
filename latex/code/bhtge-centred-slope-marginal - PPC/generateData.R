library(MASS)

extractRandomComp = function(viaReg=F){
  temp = ds$weight
  if(viaReg==T){
    reg=lm(weight~gender+by+time,data=ds)
    temp = temp-(as.numeric(ds$gender)-1)*reg$coefficients[2]
    temp = temp-(as.numeric(ds$by)-1)*reg$coefficients[3]
    #temp = temp-ds$time*reg$coefficients[4]
  }else{
    temp = temp-(as.numeric(ds$gender)-1)*40
    temp = temp-(as.numeric(ds$by)-1)*30
    #temp = temp-ds$time*10
  }
  temp

  randIntercept = numeric()
  randSlope = numeric()
  for(i in 1:nsubjects){
    startIndex = (i-1)*nrep + 1
    endIndex = i*nrep
    reg=lm(temp[startIndex:endIndex]~time)   
    randIntercept[i] = reg$coefficients[1]
    randSlope[i] = reg$coefficients[2]
  }
  
  return(cbind(randIntercept, randSlope))
}

nrep = 10
time = 1:nrep
weightgen=function(gender, by, diet){
  weight = rep(150, length(time));
  
  if(gender=="M"){
    weight = weight + 40
  }
  
  if(by=="97"){
    weight = weight + 30
  }
  
  weight = weight + time*10
  
  varcovMatrix=matrix(c(8,4.2,4.2,18),2,2)
  
  switch(diet,
         poor={
           randeff=mvrnorm(1,c(-20,-30), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
         },
         ok={
           randeff=mvrnorm(1,c(0,0), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
         },
         good={
           randeff=mvrnorm(1,c(20, 30), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
         })
  
  weight = weight + rnorm(length(time), 0, 6)
}

ncomponents = 3

if(ncomponents == 3){
  diets = c("poor", "ok", "good")
}else{
  diets = c("poor", "good", "ok", "baam", "kaam")
}

sublist = matrix(nrow = 0, ncol = 6)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
    for(diet in diets){
      for(k in 1:15){
        id = rep(paste(gender, by, diet, k, sep=""), length(time))
        weight = weightgen(gender, by, diet)
        subject = cbind(id, rep(gender, length(time)),
                        rep(by, length(time)), time,
                        rep(diet, length(time)), weight)
        sublist=rbind(sublist, subject)
      }
    }
  }
}


ds = data.frame(sublist)
colnames(ds) = c("id", "gender", "by", "time", "diet", "weight")
ds$time =  as.numeric(as.character(ds$time))
ds$weight = as.numeric(as.character(ds$weight))
ds = ds[order(ds$id, ds$time),]

dswide=data.frame(dcast(ds, id+gender+by ~ time, value.var = "weight"))
dswide_y=as.matrix(dswide[,c(-1,-2,-3)])

nobs = nrow(ds)
nsubjects = nobs/nrep

rm(subject)
rm(by)
rm(diet)
rm(gender)
rm(id)
rm(k)
rm(weight)
