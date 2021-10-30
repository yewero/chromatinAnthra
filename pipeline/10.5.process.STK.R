library(GEOquery)

if(file.exists("../data/clinical/GSE1456_raw.RData")) {
  load(file = "../data/clinical/GSE1456_raw.RData")
} else {
  rawsetL = getGEO("GSE1456")
  dir.create("../data/clinical", showWarnings = F, recursive = T)
  save(rawsetL,file = "../data/clinical/GSE1456_raw.RData")
}

if(! file.exists("../data/clinical/Sweden_Clinical.RData")) {
  download.file(url = "https://www.meb.ki.se/sites/yudpaw/wp-content/uploads/sites/5/papers/Intrinsic_Signature.zip", 
                destfile = "../data/clinical/Intrinsic_Signature.zip")
  unzip(zipfile = "../data/clinical/Intrinsic_Signature.zip", files = "Sweden_Clinical.RData", exdir = "../data/clinical/")
  file.remove("../data/clinical/Intrinsic_Signature.zip")
}
load("../data/clinical/Sweden_Clinical.RData")

if(! file.exists("../data/clinical/stockholm-U133A.RData")) {
  download.file(url = "https://www.meb.ki.se/sites/yudpaw/wp-content/uploads/sites/5/papers/stockholm-U133A.RData", 
                destfile = "../data/clinical/stockholm-U133A.RData")
}
load("../data/clinical/stockholm-U133A.RData")

expOrig = exprs(rawsetL[[1]])
prExp = expOrig["208305_at",]
her2Exp = expOrig["216836_s_at",]
erExp = expOrig["205225_at",]

getThreshold = function(exp){
  intersect <- function(m1, s1, m2, s2, prop1, prop2){
    # fix order
    if(m1 > m2) {
      mx <- m1
      m1 <- m2
      m2 <- mx
      sx <- s1
      s1 <- s2
      s2 <- sx
      propx <- prop1
      prop1 <- prop2
      prop2 <- propx
    }
    
    B <- (m1/s1^2 - m2/s2^2)
    A <- 0.5*(1/s2^2 - 1/s1^2)
    C <- 0.5*(m2^2/s2^2 - m1^2/s1^2) - log((s1/s2)*(prop2/prop1))
    
    (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
  }
  library(flexmix)
  kk = flexmix(exp~1,k=2,model = list(FLXMRglm(exp ~ ., family = "gaussian"), FLXMRglm(exp ~ ., family = "gaussian")))
  mu1 = parameters(kk)[[1]][1,1]
  mu2 = parameters(kk)[[1]][1,2]
  sigma1 = parameters(kk)[[1]][2,1]
  sigma2 = parameters(kk)[[1]][2,2]
  hist(exp,breaks = 100,col = kk@cluster)
  intersect(mu1,sigma1,mu2,sigma2,table(clusters(kk))[1]/length(exp),table(clusters(kk))[2]/length(exp))[2]
  
}

getThreshold(prExp)
thresholdPR = 4.579567
getThreshold(her2Exp)
thresholdher2 = 8.896044
getThreshold(erExp)
thresholder =8.464946
pr_probe = prExp>thresholdPR
her2_probe = her2Exp>thresholdher2
er_probe = erExp>thresholder

phenoData = pData(rawsetL[[1]])

# sample ID mapping
probeTest = intersect(rownames(exprs(rawsetL[[1]])), rownames(stk.x))
exp_geo = exprs(rawsetL[[1]])[probeTest,]
exp_ext = stk.x[probeTest,]

cor_exp <- cor(exp_geo, exp_ext)
dic_samples <- lapply(colnames(cor_exp), function(x) {
  y <- data.frame(ext_id = x, geo_id = rownames(cor_exp)[which.max(cor_exp[, x])], cor_value = max(cor_exp[, x]), stringsAsFactors = F)
  return(y)
})
dic_samples <- do.call("rbind", dic_samples)
dic_samples <- subset(dic_samples, cor_value > 0.99)
all(as.numeric(gsub("GSM", "", dic_samples[, 2])) == setdiff(seq(107072, 107231), 107116))

# fix IDs
id_mu <- data.frame(ori_id = c(rownames(subset(SwedenClinical, cohort == "Uppsala")), colnames(stk.x)), stringsAsFactors = F)
id_mu$cohort <- rep(c("Uppsala", "Stockholm"), c(sum(SwedenClinical$cohort == "Uppsala"), ncol(stk.x)))
id_mu$new_id <- make.unique(id_mu$ori_id, sep = "")
id_mu <- subset(id_mu, cohort == "Stockholm")
all(rownames(subset(SwedenClinical, cohort == "Stockholm")) %in% id_mu$new_id)
dic_samples$ext_id_mu <- id_mu$new_id[match(dic_samples$ext_id, id_mu$ori_id)]

SwedenClinical = SwedenClinical[ SwedenClinical$cohort=="Stockholm", ]

phenoData2 = cbind(phenoData,SwedenClinical[dic_samples[match(phenoData$geo_accession,dic_samples[,2]),4],],stk.clinical[dic_samples[match(phenoData$geo_accession,dic_samples[,2]),1],])
phenoData = phenoData2

#phenoData_2 = GSE1456_UPC_2@phenoData@data
#phenoData = rbind(phenoData_1[,1:ncol(phenoData_2)],phenoData_2)

covariate.df= data.frame(title=phenoData$title, stringsAsFactors = F)
rownames(covariate.df)=rownames(phenoData)
covariate.df$adjuvant = phenoData$treat.adjuv
covariate.df$neoadjuvant = F
covariate.df$adjNONE = F
covariate.df$adjNONE[phenoData$treat.adjuv==F]=T
covariate.df$adjNONE[is.na(phenoData$treat.adjuv)]=NA

covariate.df$adjNA = F
covariate.df$adjNA[is.na(phenoData$treat.adjuv)]=T



covariate.df$geoAc = phenoData$geo_accession
covariate.df$cohort = "STK"

covariate.df$age = phenoData$age
ER_1 = phenoData$er
table(ER,ER_1)
table(er_probe,ER_1[idsim])
covariate.df$er = ER_1
covariate.df$pr = PR
covariate.df$pr[idsim] = pr_probe[idsim]


covariate.df$her2 = HER2
covariate.df$her2[idsim] = her2_probe[idsim]

covariate.df$hist = NA
grade = NA

covariate.df$grade = grade


size = NA
lympNodePos = NA





stage = NA

covariate.df$stage = stage
covariate.df$t_stage = NA
covariate.df$m_stage = NA
covariate.df$n_stage= NA
covariate.df$lympNodePos = lympNodePos
covariate.df$lymphNodeNum = NA
covariate.df$size= size
covariate.df$CX = ifelse(covariate.df$adjuvant,"CMF","none")
covariate.df$CX[ is.na(covariate.df$adjuvant)]=NA
covariate.df$RX= NA
covariate.df$HX = ifelse(covariate.df$adjuvant,"Tam","none")
covariate.df$HX[ is.na(covariate.df$adjuvant)]=NA

covariate.df$anthra = F

covariate.df$chemo = covariate.df$CX

covariate.df$ht = covariate.df$HX=="Tam"
covariate.df$ht[ is.na(covariate.df$adjuvant)]=NA

covariate.df$NoTr =  covariate.df$CX=="none"
covariate.df$NoTr[ is.na(covariate.df$adjuvant)]=NA


covariate.df$rfs.e = phenoData$characteristics_ch1.1=="RELAPSE: 1"
covariate.df$rfs.t =  as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.2))*365

covariate.df$rdfs.e =NA
covariate.df$rdfs.t =NA

covariate.df$os.e = phenoData$characteristics_ch1.3=="DEATH: 1"
covariate.df$os.t = as.numeric(gsub("[^0-9,.]", "",phenoData$characteristics_ch1.5))*365
covariate.df$dss.e = phenoData$characteristics_ch1.4=="DEATH_BC: 1"
covariate.df$pCR = NA

covariate.df$Herc = NA
covariate.df$Tam =covariate.df$ht
covariate.df$Tax= NA

# process cells 

dir.create("../data/clinical/cells/GSE1456/GPL96/", recursive = T, showWarnings = F)
flist <- sapply(colnames(rawsetL[[1]]@assayData$exprs),function(x) list.files("../data/clinical/cells/GSE1456/",  paste(x,"*",sep = ""), full.names = TRUE))
file.copy(unlist(flist), "../data/clinical/cells/GSE1456/GPL96/")

dir.create("../data/clinical/cells/GSE1456/GPL97/", recursive = T, showWarnings = F)
flist <- sapply(colnames(rawsetL[[2]]@assayData$exprs),function(x) list.files("../data/clinical/cells/GSE1456/",  paste(x,"*",sep = ""), full.names = TRUE))
file.copy(unlist(flist), "../data/clinical/cells/GSE1456/GPL97/")


library(affy)

fns = list.celfiles("../data/clinical/cells/GSE1456/GPL96/",full.names=T)
celF = ReadAffy(filenames = fns)
eset <- rma(celF)
colnames(eset)=substr(colnames(eset),1,10)
pData(eset)=covariate.df

featureData(eset)=featureData(rawsetL[[1]])
eset = eset[which(rownames(eset) %in% eset@featureData@data$ID),]
dir.create("../data/clinical/rmas/GSE1456",recursive = T)
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/rmas/GSE1456/rmaData_96.RData")

fns = list.celfiles("../data/clinical/cells/GSE1456/GPL97/",full.names=T)
celF = ReadAffy(filenames = fns)
eset <- rma(celF)
colnames(eset)=substr(colnames(eset),1,8)
pData(eset)=cov

featureData(eset)=featureData(rawsetL[[2]])
#eset = eset[which(rownames(eset) %in% eset@featureData@data$ID),]
dir.create("../data/clinical/rmas/GSE1456",recursive = T)
eset_rma = eset
save(eset_rma,file="../data/clinical/rmas/GSE1456/rmaData_97.RData")
