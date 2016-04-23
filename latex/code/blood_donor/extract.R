############ extracting random components to plot them ##############
extractRandomComp = function(){
  temp = ds$Hb
  reg=lm(Hb~Age+Season+TSPD+Donate+donationLastTwoYears,data=ds)
  temp = temp-ds$firstAge*reg$coefficients[2]
  temp = temp-(as.numeric(ds$Season)-1)*reg$coefficients[3]
  temp = temp-ds$TSPD*reg$coefficients[4]
  temp = temp-ds$Hb.PD*reg$coefficients[5]
  
  randIntercept = numeric()
  randSlope = numeric()
  prevEnd = 0
  k=0
  for(i in unique(ds$Id)){
    startIndex = prevEnd + 1
    endIndex = prevEnd + length(timePerSubject[[paste(i)]])
    
    reg=lm(temp[startIndex:endIndex]~I(donationLast2YearsPerSubject[[paste(i)]]/10))   
    randIntercept[k] = reg$coefficients[1]
    randSlope[k] = reg$coefficients[2]
    k=k+1
    
    prevEnd = endIndex
  }
  return(cbind(randIntercept, randSlope))
}

randomComp=extractRandomComp()
randomCompDf = data.frame(randomComp[complete.cases(randomComp)==TRUE,])
qplot(x=randIntercept, y=randSlope, data=randomCompDf)
var(randomCompDf)
#plot(density(extractRandomComp()[,1]))
