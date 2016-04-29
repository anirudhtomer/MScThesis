resultDf = data.frame(run=integer(0),ncomp=integer(0),DIC1=numeric(0), DIC2=numeric(0), DIC3=numeric(0), DIC7=numeric(0))

saveResult=function(ncomp,dic1,dic2,dic3,dic7){
  index = nrow(resultDf) + 1
  resultDf[index,] = c(index,ncomp,dic1,dic2,dic3,dic7)
  assign('resultDf',resultDf,envir=.GlobalEnv)
  write.csv(resultDf,file = "result.csv",row.names = F)  
}