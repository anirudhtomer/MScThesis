variableList=c("Age","Donate", "TSPD", "Season","donationLast2Years")
interactions = vector();
k=1;
for(i in 1:length(variableList)){
  for(j in i:length(variableList)){
    interactions[k] = paste(variableList[i], "*", variableList[j])
    k=k+1
  }
}
write.csv(data.frame(interactions), file.choose(),quote = F)