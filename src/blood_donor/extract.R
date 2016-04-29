############ extracting random components to plot them ##############
extractRandomComp = function(){
  temp = ds$Hb
  reg=lm(Hb~Age*Season + TSPD*Donate+ TSPD*Season + 
           donationLast2Years*TSPD + donationLast2Years*Donate + 
           donationLast2Years*Season + I(donationLast2Years^2),data=ds)
  #reg=lm(Hb~Age+Donate + Season+TSPD + donationLast2Years,data=ds)
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

# reg=lm(Hb~Age*Donate+ Donate*Season+TSPD*donationLast2Years + donationLast2Years*Donate,data=ds)
# modelMatrix = model.matrix(reg)[,-1]
# 
# corrplot(round(cor(modelMatrix),2), order="original", method="circle", tl.pos="lt", type="lower",
#          tl.col="black", tl.cex=0.6, tl.srt=45,
#          addCoef.col="black", addCoefasPercent = F,
#          p.mat = 1-abs(cor(modelMatrix)), sig.level=0.3, insig = "blank")
