# TMN
tmn_stage = function(t,m,n){
  tmn_table = read.table("../data/clinical/TMN_breast.tsv.csv",header = T,sep = "\t",stringsAsFactors = F)
  if(is.na(m)|is.na(t)|is.na(n))
    return(NA)
  if(m==1)
    return(4)
  else{
    return(tmn_table$stage[ which(tmn_table$T==t & tmn_table$N==n)])
  }
}
