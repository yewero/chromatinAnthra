# generate Chromatin Regulator Genes list

# complexes file (v1.2) was generated on Feb 2018. Changes in GO modify the final list
# of CRGs. For reproductibility purposes, please use file 1.2

library(biomaRt)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


complexes = NULL

# 1. methyl transferases
HMT <- getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0018024', mart = ensembl)
HMT = HMT[HMT$go_id=='GO:0018024',]
HMT = HMT[HMT$chromosome_name %in% c(1:22,"X","Y"),]

HMT$name="METHTRANSF"
HMT$color="#D9D9D9"

complexes = rbind(complexes,HMT)

# 2. histone demethylases
HDMT = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0032452', mart = ensembl)
HDMT = HDMT[HDMT$go_id=='GO:0032452',]
HDMT = HDMT[HDMT$chromosome_name %in% c(1:22,"X","Y"),]

HDMT$name="DEMETH"
HDMT$color="#FFED6F"

complexes = rbind(complexes,HDMT)

# 3. histone deacetylase
HDA =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
              filters = 'go', values = 'GO:0004407', mart = ensembl)
HDA = HDA[HDA$go_id=='GO:0004407',]
HDA = HDA[HDA$chromosome_name %in% c(1:22,"X","Y"),]

HDA$name="DEACETH"
HDA$color="darkorchid1"

complexes = rbind(complexes,HDA)

# 4. histone acetyltransferase
HA =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0004402', mart = ensembl)
HA = HA[HA$go_id=='GO:0004402',]
HA = HA[HA$chromosome_name %in% c(1:22,"X","Y"),]

HA$name="ACETH"
HA$color="cornsilk"

complexes = rbind(complexes,HA)

# 5. histone phosporilation
HP =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0016572', mart = ensembl)
HP = HP[HP$go_id=='GO:0016572',]
HP = HP[HP$chromosome_name %in% c(1:22,"X","Y"),]

HP$name="PHOSPH"
HP$color="brown1"

complexes = rbind(complexes,HP)

# 6. PRC1 comples
PRC1 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0035102', mart = ensembl)
PRC1 = PRC1[PRC1$go_id=='GO:0035102',]
PRC1 = PRC1[PRC1$chromosome_name %in% c(1:22,"X","Y"),]
PRC1$name="PRC1"
PRC1$color="#FFFFB3"

complexes = rbind(complexes,PRC1)

# 7. PRC2 complex
PRC2 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0035098', mart = ensembl)
PRC2 = PRC2[PRC2$go_id=='GO:0035098',]
PRC2 = PRC2[PRC2$chromosome_name %in% c(1:22,"X","Y"),]

PRC2$name="PRC2"
PRC2$color="#BEBADA"

complexes = rbind(complexes,PRC2)

# REMODELERS FAMILY according to http://www.nature.com/nrm/journal/v7/n6/full/nrm1945.html
#SWI/SNF, ISWI, NURD/Mi-2/CHD, INO80 and SWR1.

# 8. BAF (+)
# Tumour suppressor, Differentiation, Development, Elongation, Signalling, Splicing

BAF =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
              filters = 'go', values = 'GO:0016514', mart = ensembl)
BAF = BAF[BAF$go_id=='GO:0016514',]
otherBAF = c("PBRM1","SS18","BCL11A","BCL11B","BCL7A","BCL7B","BCL7C","DPF1","DPF2","DPF3","PHF10","BRD9")
otherBAF <- setdiff(otherBAF, BAF$hgnc_symbol)
BAF2 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = otherBAF, mart = ensembl)
if(length(setdiff(otherBAF, BAF2$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
BAF2$go_id = BAF$go_id[1]
BAF2$name_1006 = BAF$name_1006[1]
BAF=rbind(BAF,BAF2)
BAF = BAF[BAF$chromosome_name %in% c(1:22,"X","Y"),]

BAF$name = "SWI_SNF"
BAF$color = "#8DD3C7"

complexes = rbind(complexes,BAF)

## 9. ISWI  (NURF / ACF / CHRAC  /WICH / NORC /RSF / CERF) (+)
# NURF: SMARCA1 + BPTF
# WICH: SMARCA5 + WSTF
# NORC: SMARCA5 + TIP5 (BAZ2A) + BAZ2B
# RSF: SMARCA5 + RSF1
# ACF: ACF1 + SMARCA5
# CHRAC: SMARCA5 + ACF1 + CHRAC15 + CHRAC17 
# CERF: SMARCA1 + CECR2

iswi_genes = c("SMARCA5","SMARCA1","BPTF","BAZ1B","BAZ2A","BAZ2B","RSF1","BAZ1A","CHRAC1","POLE3","CECR2","C17orf49","RBBP4","HMGXB4","USF1")
ISWI = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = iswi_genes, mart = ensembl)
if(length(setdiff(iswi_genes, ISWI$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
ISWI$go_id = ""
ISWI$name_1006 = "ISWI"
ISWI = ISWI[ISWI$chromosome_name %in% c(1:22,"X","Y"),]

ISWI$name="ISWI"
ISWI$color = "#FDB462"

complexes = rbind(complexes,ISWI)

# 10. CHDs NURD/Mi-2 (+)
# Transcriptional repression and silencing, Development

chdsGenes = c("CHD1","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9","HDAC1","HDAC2","MBD2","MBD3","MTA1", "MTA2","MTA3","GATAD2A","GATAD2B", "RBBP4","RBBP7")
CHD = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
            filters = 'hgnc_symbol', values = chdsGenes, mart = ensembl)
if(length(setdiff(chdsGenes, CHD$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
CHD$go_id = ""
CHD$name_1006 = "CHD_NURD_mi2"
CHD = CHD[CHD$chromosome_name %in% c(1:22,"X","Y"),]

CHD$name="CHD_NURD"
CHD$color = "#FB8072"

complexes = rbind(complexes,CHD)

# 11. INO80 (+)
#
# http://www.nature.com/nrm/journal/v10/n6/fig_tab/nrm2693_T1.html
INO80 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
                filters = 'go', values = 'GO:0031011', mart = ensembl)
INO80 = INO80[INO80$go_id=='GO:0031011',]

ino80_genes = c("ACTB")
INO80_2 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = ino80_genes, mart = ensembl)
if(length(setdiff(ino80_genes, INO80_2$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
INO80_2$go_id = INO80$go_id[1]
INO80_2$name_1006 = INO80$name_1006[1]
INO80=rbind(INO80,INO80_2)
INO80 = INO80[INO80$chromosome_name %in% c(1:22,"X","Y"),]

INO80$name="INO80"
INO80$color="cadetblue2"

complexes = rbind(complexes,INO80)

# 12. SWR1 (+)
# DNA repair
# http://www.nature.com/nrm/journal/v10/n6/fig_tab/nrm2693_T1.html
swr_genes = c("SRCAP","RUVBL1","RUVBL2","ACTB","ACTL6A","ACTR6","YEATS4","DMAP1","ERCC5","VPS72")
SWR1 = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
             filters = 'hgnc_symbol', values = swr_genes, mart = ensembl)
if(length(setdiff(swr_genes, SWR1$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
SWR1 = SWR1[SWR1$chromosome_name %in% c(1:22,"X","Y"),]

SWR1$go_id = ""
SWR1$name_1006 = "SWR1 complex"

SWR1$name="SWR1"
SWR1$color="goldenrod4"

complexes = rbind(complexes,SWR1)

# 13. BAP1/PR-DUB (+)
bap1_bpdub_complex_genes = c("ASXL1","ASXL2","BAP1","BARD1","BRCA1","FOXK1","FOXK2","HCFC1","KDM1B","OGT","YY1")
BAP1_PRDUB = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                   filters = 'hgnc_symbol', values = bap1_bpdub_complex_genes, mart = ensembl)
if(length(setdiff(bap1_bpdub_complex_genes, BAP1_PRDUB$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
BAP1_PRDUB$go_id=""
BAP1_PRDUB$name_1006="BAP1_PR-DUB_Complex"
BAP1_PRDUB = BAP1_PRDUB[BAP1_PRDUB$chromosome_name %in% c(1:22,"X","Y"),]

BAP1_PRDUB$name="BAP1_PRDUB"
BAP1_PRDUB$color="#B3DE69"

complexes = rbind(complexes,BAP1_PRDUB)


# 14. CAF-1
CAF1 =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
               filters = 'go', values = 'GO:0033186', mart = ensembl)

CAF1 = CAF1[CAF1$go_id=='GO:0033186',]

CAF1$name="CAF1"
CAF1$color="#80B1D3"

complexes = rbind(complexes,CAF1)

# 15. cohesin (+)
# http://www.nature.com/nrc/journal/v14/n6/full/nrc3743.html
#http://www.nature.com/nrg/journal/v15/n4/full/nrg3663.html
# check role CTCF in  methylation
cohe_genes = c("STAG1","STAG2","STAG3","SMC1A","SMC1B","SMC3","RAD21","REC8","RAD21L1","WAPL","PDS5A","PDS5B","CDCA5","MTA3")
COHE =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
               filters = 'hgnc_symbol', values = cohe_genes, mart = ensembl)
if(length(setdiff(cohe_genes, COHE$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
COHE$go_id=""
COHE$name_1006="Cohesin_Complex"
COHE <- COHE[COHE$chromosome_name %in% c(1:22,"X","Y"),]

COHE$name="COHE"
COHE$color="#BC80BD"

complexes = rbind(complexes,COHE)

# 16. condensins (+)
#http://genesdev.cshlp.org/content/26/15/1659.full
conde_genes = c("SMC4","NCAPH","NCAPG","NCAPD2","SMC2","NCAPD3","NCAPG2","NCAPH2")
CONDE =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = conde_genes, mart = ensembl)
if(length(setdiff(conde_genes, CONDE$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
CONDE$go_id=""
CONDE$name_1006="Condensin_Complex"

CONDE$name="CONDE"
CONDE$color="#CCEBC5"

complexes = rbind(complexes,CONDE)

# 17. TOPo activity
TOPO = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
             filters = 'go', values = 'GO:0003916', mart = ensembl)
TOPO = TOPO[TOPO$go_id=='GO:0003916',]
TOPO = TOPO[TOPO$chromosome_name %in% c(1:22,"X","Y"),]

TOPO$name="TOPO"
TOPO$color="#FCCDE5"

complexes = rbind(complexes,TOPO)


# 18-1. DNA methyltransferases

DNA_meth = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
                   filters = 'go', values = 'GO:0006306', mart = ensembl)
DNA_meth = DNA_meth[DNA_meth$go_id=='GO:0006306',]
DNA_meth = DNA_meth[DNA_meth$chromosome_name %in% c(1:22,"X","Y"),]

DNA_meth$name="DNA_METH"
DNA_meth$color="aquamarine"

complexes = rbind(complexes,DNA_meth)


# 18-2. DNA demethylases
DNA_demeth = getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position", 'go_id',"name_1006"),
      filters = 'go', values = 'GO:0080111', mart = ensembl)
DNA_demeth = DNA_demeth[DNA_demeth$go_id=='GO:0080111',]
DNA_demeth = DNA_demeth[DNA_demeth$chromosome_name %in% c(1:22,"X","Y"),]

DNA_demeth$name="DNA_DMETH"
DNA_demeth$color="slateblue1"

complexes = rbind(complexes,DNA_demeth)

# 18-3. histone proteins (+)
hist_genes = read.table("../data/histone_genenames.txt",header = T,sep = "\t",stringsAsFactors = F)

# HIST =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
#                filters = 'hgnc_symbol', values = hist_genes$Approved.Symbol, mart = ensembl)
###
HIST =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position","hgnc_id"),
               filters = 'hgnc_id', values = paste0("HGNC:", hist_genes$HGNC.ID), mart = ensembl)
# subset(hist_genes, ! HGNC.ID %in% gsub("^HGNC:", "", HIST$hgnc_id))
# tmp <- data.frame("H2BW4P", 767811, "X", 103975924, 103978357, "HGNC:25757", stringsAsFactors = F)
# colnames(tmp) <- colnames(HIST)
# HIST <- rbind(HIST, tmp)
# rm(tmp)
HIST <- HIST[HIST$chromosome_name %in% c(1:22,"X","Y"),]
HIST <- HIST[match(paste0("HGNC:", hist_genes$HGNC.ID), HIST$hgnc_id), ]
HIST$hgnc_symbol_new <- HIST$hgnc_symbol
HIST$hgnc_symbol <- hist_genes$Approved.Symbol
HIST <- HIST[! is.na(HIST$chromosome_name), ]
HIST_ext <- HIST
HIST <- HIST[, 1:5]
###
if(length(setdiff(hist_genes$Approved.Symbol, HIST$hgnc_symbol)) > 0) {
  warning("Missed genes: ", setdiff(hist_genes$Approved.Symbol, HIST$hgnc_symbol))
}
HIST$go_id=""
HIST$name_1006="HISTONES"

HIST$name="HIST"
HIST$color="tan3"
complexes = rbind(complexes,HIST)

# 18-4. chromatin pioneer factors (+)
pio_genes = c("FOXA1")

pio =   getBM(attributes=c('hgnc_symbol',"entrezgene_id", "chromosome_name","start_position","end_position"),
                filters = 'hgnc_symbol', values = pio_genes, mart = ensembl)
if(length(setdiff(pio_genes, pio$hgnc_symbol)) > 0) {
  warning("Missed genes")
}
pio$go_id=""
pio$name_1006="pio_GENES"

pio$name="pio"
pio$color="violetred1"

complexes = rbind(complexes,pio)

# X. all
if(! exists("complexes_full")) {
  cat("[info] create full copy of complexes\n")
  complexes_full <- complexes
}
complexes <- complexes_full

# remove duplicate entrez
complexes2 = unique(complexes[order(complexes$entrezgene_id),c(1,2)])
idDup = duplicated(complexes2$hgnc_symbol)
falseEntrez = complexes2$entrezgene_id[idDup]
rm(complexes2)

complexes = complexes[!complexes$entrezgene_id %in% falseEntrez,]
complexes = na.omit(complexes) # 471 unique genes (! it will remove the genes with NA in entrez gene ID)
#subset(complexes_full, ! hgnc_symbol %in% complexes_rep$hgnc_symbol)
complexes_rep <- complexes


#save(complexes, file="../data/complexes_v1.3_BRCA.RData")
load("../data/complexes_v1.2_BRCA.RData")
complexes_pub <- complexes
rm(complexes)

### comparison
setdiff(complexes_pub$hgnc_symbol, complexes_rep$hgnc_symbol)
setdiff(complexes_rep$hgnc_symbol, complexes_pub$hgnc_symbol)

group_stat <- merge(as.data.frame(table(complexes_pub$name)), as.data.frame(table(complexes_rep$name)), by = "Var1", all = T)
colnames(group_stat) <- c("name", "pub", "rep")
subset(group_stat, pub != rep)

chr_stat <- merge(as.data.frame(table(complexes_pub$chromosome_name)), as.data.frame(table(complexes_rep$chromosome_name)), by = "Var1", all = T)
colnames(chr_stat) <- c("chr", "pub", "rep")
subset(chr_stat, pub != rep)
###
