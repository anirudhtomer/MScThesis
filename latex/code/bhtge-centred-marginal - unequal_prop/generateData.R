extractRandomComp = function(viaReg=F){
  temp = ds$weight
  if(viaReg==T){
    reg=lm(weight~gender+by+time,data=ds)
    temp = temp-(as.numeric(ds$gender)-1)*reg$coefficients[2]
    temp = temp-(as.numeric(ds$by)-1)*reg$coefficients[3]
    temp = temp-ds$time*reg$coefficients[4]
  }else{
    temp = temp-(as.numeric(ds$gender)-1)*40
    temp = temp-(as.numeric(ds$by)-1)*30
    temp = temp-ds$time*10
  }
  temp
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
  
  precision = 0.0625
  sd = sqrt(1/precision)
  
  switch(diet,
         poor={
           weight = weight + rnorm(1, -30, sd)
         },
         ok={
           weight = weight + rnorm(1, 0, sd)
         },
         good={
           weight = weight + rnorm(1, 30, sd)
         })
  
  weight = weight + rnorm(length(time), 0, 3)
}

ncomponents = 3

if(ncomponents == 3){
  diets = c("poor", "ok", "good")
}else{
  diets = c("poor", "good")
}

dietProp = list("poor"=0.6,"ok"=0.3,"good"=0.1)

sublist = matrix(nrow = 0, ncol = 6)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
    for(diet in diets){
      for(k in 1:round(45*dietProp[[diet]])){
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
