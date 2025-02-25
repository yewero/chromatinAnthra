library(org.Hs.eg.db)
library(RTN)

# load gene ID
gene_ID_v22 <- read.table("../../Genomes/human_v22/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_ID_v22) <- c("id", "name")
gene_ID_v22$ensembl_nov <- gsub("\\..*", "", gene_ID_v22$id)
gene_ID_v38 <- read.table("../../Genomes/human_v38/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(gene_ID_v38) <- c("id", "name")
gene_ID_v38$ensembl_nov <- gsub("\\..*", "", gene_ID_v38$id)

load("../data/complexes_v1.2_BRCA.RData")
# check gene names
complexes$hgnc_symbol[complexes$hgnc_symbol == "HIST2H2AA4"] <- "H2AC19"
complexes$hgnc_symbol[complexes$hgnc_symbol == "HIST2H2AA3"] <- "H2AC18"
complexes <- merge(complexes, gene_ID_v38, by.x = "hgnc_symbol", by.y = "name", all.x = T, sort = F)
complexes_1 <- subset(complexes, ! is.na(id))
complexes_2 <- merge(subset(complexes, is.na(id))[, 1:9], gene_ID_v22, by.x = "hgnc_symbol", by.y = "name", all.x = T, sort = F)
complexes <- rbind(complexes_1, complexes_2)
complexes <- merge(complexes, gene_ID_v22, by = "ensembl_nov", sort = F)
colnames(complexes)[colnames(complexes) == "hgnc_symbol"] <- "hgnc_symbol_ori"
colnames(complexes)[colnames(complexes) == "name.y"] <- "hgnc_symbol"
colnames(complexes)[colnames(complexes) == "name.x"] <- "name"
complexes <- complexes[, c(13, 3:10, 1:2)]
all(complexes$hgnc_symbol %in% gene_ID_v22$name)
rm(complexes_1, complexes_2)
cat("Unique selected gene number:", length(unique(complexes$hgnc_symbol)), "\n")

# load expression object
load("../data/vst.RData")

exp.T.vst.e = exp.T.vst
exp.N.vst.e = exp.N.vst
rownames(exp.T.vst.e)=Des$external_gene_name
rownames(exp.N.vst.e)=Des$external_gene_name

# des_raw = select(org.Hs.eg.db, rownames(exp.T.vst), c("ENTREZID","GENENAME","ENSEMBL"), "ENSEMBL")
# Des$entrez = des_raw$ENTREZID[match(Des$ensembl_gene_id,des_raw$ENSEMBL)]
# idna = which(is.na(Des$entrez))
# idna = which(rownames(exp.T.vst.e)=="")
Des.e = Des
# exp.T.vst.e=exp.T.vst.e[-idna,]
# exp.N.vst.e=exp.N.vst.e[-idna,]
# Des.e = Des.e[-idna,]
iddup = duplicated(rownames(exp.T.vst.e))
exp.T.vst.e= exp.T.vst.e[which(!iddup),]
exp.N.vst.e= exp.N.vst.e[which(!iddup),]
Des.e = Des.e[which(!iddup),]
geneNames = intersect(Des.e$external_gene_name, complexes$hgnc_symbol)
names(geneNames)=geneNames
tumor.TNI <- new("TNI", gexp=exp.T.vst.e, regulatoryElements=geneNames)
tumor.rtni <- tni.preprocess(tumor.TNI)


tumor.rtni<-tni.permutation(tumor.rtni)
tumor.rtni<-tni.bootstrap(tumor.rtni)
tumor.rtni<-tni.dpi.filter(tumor.rtni)

save(tumor.rtni,file="../data/allTumor_rtni.RData")

#aracne
load("../data/allTumor_rtni.RData")


tpc2regulon = function(net,filterTFs=NULL,pvalueCutoff= 0.05){
    #net = normal.rtni
    
    if (is.null(filterTFs)) {TFs <- colnames(net@results$tn.dpi)}
    if (!is.null(filterTFs)) {TFs <- filterTFs; TFs <- intersect(TFs, colnames(net@results$tn.dpi))}
    
    
    mode1 = sign(net@results$tn.dpi[,TFs])
    
    scores <- net@results$tn.dpi[,TFs]
    qvalues = net@results$miadjpv[,TFs]
    
    
    aracne <- list()
    for (tf in TFs) {
      reg <- qvalues[,tf]
      which.reg <- reg <= pvalueCutoff
      which.mode = mode1[,tf]!=0
      likelihood <- scores[which.reg&which.mode,match(tf, colnames(scores))]
      tfmode <- mode1[which.reg&which.mode,match(tf, colnames(mode1))]
      aracne[[tf]] <- list("tfmode"=tfmode, "likelihood"=likelihood)
    }
    
    # removing missing data from the aracne regulon
    aracne <- aracne[names(aracne) != "NA"]
    aracne <- lapply(aracne, function(x) {
      filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
      x$tfmode <- x$tfmode[filtro]
      x$likelihood <- x$likelihood[filtro]
      return(x)
    })
    
    regul <- aracne[sapply(aracne, function(x) length(names(x$tfmode)))>0]
    class(regul) <- "regulon"
    
    return(regul)
  }


regul2.T = tpc2regulon(tumor.rtni)
chromatin_regulon = regul2.T
save(chromatin_regulon,file = "../data/chromatin_regulon.RData")





#Supfig2
library("ggplot2")
dir.create("../plots", showWarnings = F)

chromatin_degrees = unlist(lapply(chromatin_regulon,function(x) length(x$tfmode)))
# get histogram of regulon
toPlot = data.frame(gene=names(chromatin_degrees),degree=chromatin_degrees,stringsAsFactors = F)

ggpubr::gghistogram(toPlot,x="degree",y="..count..",xlab = "Degree",rug = TRUE ,add="median",color=viridis::viridis_pal()(3)[2],fill=viridis::viridis_pal()(3)[2] )
ggsave("../plots/SupFig2_histogram_CRG_regulon.pdf",width = 7,height = 4)
# get barplot regulon

ggpubr::ggbarplot(toPlot[toPlot$degree>500,],x="gene",y="degree",sort.val = "asc",
                  color=viridis::viridis_pal()(3)[1],fill=viridis::viridis_pal()(3)[1])+theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("../plots/SupFig2_barplot_CRG_regulon.pdf",width = 7,height = 4)
