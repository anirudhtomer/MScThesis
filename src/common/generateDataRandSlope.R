
extractRandomComp = function(viaReg=F){
  temp = ds$weight
  if(viaReg==T){
    reg=lm(weight~gender+by+age+time,data=ds)
    temp = temp-(as.numeric(ds$gender)-1)*reg$coefficients[2]
    temp = temp-(as.numeric(ds$by)-1)*reg$coefficients[3]
    temp = temp-ds$age*reg$coefficients[4]
    #temp = temp-ds$time*reg$coefficients[4]
  }else{
    temp = temp-(as.numeric(ds$gender)-1)*40
    temp = temp-(as.numeric(ds$by)-1)*30
    temp = temp-ds$age*reg$coefficients[4]
    #temp = temp-ds$time*10
  }

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
weightgen=function(gender, by, diet, age){
  weight = rep(150, length(time));
  
  if(gender=="M"){
    weight = weight + 20
  }
  
  if(by=="97"){
    weight = weight + 10
  }
  
  weight = weight + age*6
  
  weight = weight + time*10
  
  varcovMatrix=matrix(c(12,4.2,4.2,20),2,2)
  
  switch(diet,
         poor={
           randeff=mvrnorm(1,c(-17,-15), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
           splaash1 <<- rbind(splaash1, randeff)
         },
         ok={
           randeff=mvrnorm(1,c(0,0), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
           splaash2 <<- rbind(splaash2, randeff)
         },
         good={
           randeff=mvrnorm(1,c(17, 15), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
           splaash3 <<- rbind(splaash3, randeff)
         },
         baam={
           randeff=mvrnorm(1,c(0, 30), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
           splaash4 <<- rbind(splaash4, randeff)
         },
         kaam={
           randeff=mvrnorm(1,c(-30, 0), varcovMatrix)
           weight = weight + randeff[1] + randeff[2]*time
           splaash5 <<- rbind(splaash5, randeff)
         })
  
  weight = weight + rnorm(length(time), 0, 2)
}

splaash1 = matrix(nrow=0, ncol=2)
splaash2 = matrix(nrow=0, ncol=2)
splaash3 = matrix(nrow=0, ncol=2)
splaash4 = matrix(nrow=0, ncol=2)
splaash5 = matrix(nrow=0, ncol=2)

ncomponents = 3

if(ncomponents == 3){
  diets = c("poor", "ok", "good")
}else if(ncomponents==5){
  diets = c("poor", "good", "ok", "baam", "kaam")
}else if(ncomponents == 1){
  diets = c("ok")
}

dietSubjects = c("poor"=1,"ok"=10,"good"=10, "baam"=12, "kaam"=4)

sublist = matrix(nrow = 0, ncol = 7)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
    for(diet in diets){
      for(k in 1:dietSubjects[diet]){
        id = rep(paste(gender, by, diet, k, sep=""), length(time))
        age = round(runif(1, min = 0, max=16), 2)
        weight = weightgen(gender, by, diet, age)
        subject = cbind(id, rep(age, length(time)),rep(gender, length(time)),
                        rep(by, length(time)), time,
                        rep(diet, length(time)), weight)
        sublist=rbind(sublist, subject)
      }
    }
  }
}


ds = data.frame(sublist)
colnames(ds) = c("id", "age" ,"gender", "by", "time", "diet", "weight")
ds$time =  as.numeric(as.character(ds$time))
ds$weight = as.numeric(as.character(ds$weight))
ds$age = as.numeric(as.character(ds$age))
ds = ds[order(ds$id, ds$time),]

dswide=data.frame(dcast(ds, id+gender+by+age ~ time, value.var = "weight"))
dswide_y=as.matrix(dswide[,c(-1,-2,-3,-4)])

nobs = nrow(ds)
nsubjects = nobs/nrep

rm(subject)
rm(by)
rm(diet)
rm(gender)
rm(id)
rm(k)
rm(weight)
rm(sublist)
