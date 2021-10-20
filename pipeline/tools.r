
fixCompName <- function(complexes) {
  # load gene ID
  gene_ID_v22 <- read.table("../../Genomes/human_v22/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
  colnames(gene_ID_v22) <- c("id", "name")
  gene_ID_v22$ensembl_nov <- gsub("\\..*", "", gene_ID_v22$id)
  gene_ID_v38 <- read.table("../../Genomes/human_v38/gene_ID2Name_fixed.txt", header = F, sep = "\t", stringsAsFactors = F)
  colnames(gene_ID_v38) <- c("id", "name")
  gene_ID_v38$ensembl_nov <- gsub("\\..*", "", gene_ID_v38$id)
  
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
  print(all(complexes$hgnc_symbol %in% gene_ID_v22$name))
  return(complexes)
}
