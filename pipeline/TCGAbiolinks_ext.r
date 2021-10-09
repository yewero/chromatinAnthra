
readTranscriptomeProfiling_new <- function (files, data.type, workflow.type, cases, summarizedExperiment, nCPU = 10) 
{
  if (grepl("Gene Expression Quantification", data.type, ignore.case = TRUE)) {
    if (grepl("HTSeq", workflow.type)) {
      # pb <- txtProgressBar(min = 0, max = length(files), 
      #                      style = 3)
      # for (i in seq_along(files)) {
      #   data <- read_tsv(file = files[i], col_names = FALSE, 
      #                    col_types = cols(X1 = col_character(), X2 = col_double()))
      #   if (!missing(cases)) 
      #     colnames(data)[2] <- cases[i]
      #   if (i == 1) 
      #     df <- data
      #   if (i != 1) 
      #     df <- merge(df, data, by = colnames(df)[1], 
      #                 all = TRUE)
      #   setTxtProgressBar(pb, i)
      # }
      # close(pb)
      ### {XXX
      cat("[info] read files using", nCPU, "CPUs.\n")
      cl <- parallel::makeCluster(nCPU, type = "FORK")
      df_LS <- parallel::parLapply(cl, files, function(x) {
        data <- read_tsv(file = x, col_names = F, col_types = cols(X1 = col_character(), X2 = col_double()))
        data <- as.data.frame(data)
        rownames(data) <- data[, 1]
        data <- data[, -1, drop = F]
        return(data)
      })
      parallel::stopCluster(cl); rm(cl)
      cat("[info] all files have been read.\n")
      res_checkrn <- unique(sapply(df_LS, function(x) { (nrow(x) == nrow(df_LS[[1]])) & all(rownames(x) == rownames(df_LS[[1]])) }))
      if(length((res_checkrn) == 1) & res_checkrn) {
        df <- do.call("cbind", df_LS)
        colnames(df) <- cases
        df <- data.frame(X1 = rownames(df), df, stringsAsFactors = F)
        rownames(df) <- NULL
      } else {
        stop("Inconsistent row names.")
      }
      ### XXX}
      if (summarizedExperiment) 
        df <- makeSEfromTranscriptomeProfiling(df, cases, 
                                               workflow.type)
    }
  }
  else if (grepl("miRNA", workflow.type, ignore.case = TRUE) & 
           grepl("miRNA", data.type, ignore.case = TRUE)) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
      data <- read_tsv(file = files[i], col_names = TRUE, 
                       col_types = "cidc")
      if (!missing(cases)) 
        colnames(data)[2:ncol(data)] <- paste0(colnames(data)[2:ncol(data)], 
                                               "_", cases[i])
      if (i == 1) 
        df <- data
      if (i != 1) 
        df <- merge(df, data, by = colnames(df)[1], all = TRUE)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  else if (grepl("Isoform Expression Quantification", data.type, 
                 ignore.case = TRUE)) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
      data <- read_tsv(file = files[i], col_names = TRUE, 
                       col_types = c("ccidcc"))
      if (!missing(cases)) 
        data$barcode <- cases[i]
      else data$file <- i
      if (i == 1) 
        df <- data
      if (i != 1) 
        df <- rbind(df, data)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  return(df)
}
environment(readTranscriptomeProfiling_new) <- asNamespace("TCGAbiolinks")
assignInNamespace("readTranscriptomeProfiling", readTranscriptomeProfiling_new, "TCGAbiolinks")
rm(readTranscriptomeProfiling_new)

colDataPrepare_new <- function (barcode) 
{
  message("Starting to add information to samples")
  if (all(grepl("TARGET", barcode))) 
    ret <- colDataPrepareTARGET(barcode)
  if (all(grepl("TCGA", barcode))) 
    ret <- colDataPrepareTCGA(barcode)
  message(" => Add clinical information to samples")
  patient.info <- NULL
  patient.info <- tryCatch({
    step <- 100
    for (i in 0:(ceiling(length(ret$patient)/step) - 1)) {
      start <- 1 + step * i
      end <- ifelse(((i + 1) * step) > length(ret$patient), 
                    length(ret$patient), ((i + 1) * step))
      if (is.null(patient.info)) {
        patient.info <- getBarcodeInfo(ret$patient[start:end])
      }
      else {
        patient.info <- rbind(patient.info, getBarcodeInfo(ret$patient[start:end]))
      }
    }
    patient.info
  }, error = function(e) {
    step <- 2
    for (i in 0:(ceiling(length(ret$patient)/step) - 1)) {
      start <- 1 + step * i
      end <- ifelse(((i + 1) * step) > length(ret$patient), 
                    length(ret$patient), ((i + 1) * step))
      if (is.null(patient.info)) {
        patient.info <- getBarcodeInfo(ret$patient[start:end])
      }
      else {
        patient.info <- rbind(patient.info, getBarcodeInfo(ret$patient[start:end]))
      }
    }
    patient.info
  })
  ret <- merge(ret, patient.info, by.x = "patient", by.y = "submitter_id", 
               all.x = TRUE)
  ret <- addFFPE(ret)
  if (!"project_id" %in% colnames(ret)) {
    aux <- getGDCprojects()[, 5:6]
    aux <- aux[aux$disease_type == unique(ret$disease_type), 
               2]
    ret$project_id <- as.character(aux)
  }
  if (all(grepl("TARGET", barcode))) {
    ret <- ret[match(barcode, ret$barcode), ]
    rownames(ret) <- ret$barcode
    return(ret)
  }
  ret$sample.aux <- substr(ret$sample, 1, 15)
  out <- NULL
  for (proj in na.omit(unique(ret$project_id))) {
    if (grepl("TCGA", proj, ignore.case = TRUE)) {
      message(" => Adding subtype information to samples")
      tumor <- gsub("TCGA-", "", proj)
      available <- c("ACC", "BRCA", "BLCA", "CESC", "CHOL", 
                     "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", 
                     "KIRP", "LGG", "LUAD", "LUSC", "PAAD", "PCPG", 
                     "PRAD", "READ", "SKCM", "SARC", "STAD", "THCA", 
                     "UCEC", "UCS", "UVM")
      if (grepl(paste(c(available, "all"), collapse = "|"), 
                tumor, ignore.case = TRUE)) {
        subtype <- TCGAquery_subtype(tumor)
        colnames(subtype) <- paste0("subtype_", colnames(subtype))
        if (all(str_length(subtype$subtype_patient) == 
                12)) {
          subtype$sample.aux <- paste0(subtype$subtype_patient, 
                                       "-01")
        }
        ret.aux <- ret[ret$sample.aux %in% subtype$sample.aux, 
        ]
        ret.aux <- merge(ret.aux, subtype, by = "sample.aux", 
                         all.x = TRUE)
        out <- rbind.fill(as.data.frame(out), as.data.frame(ret.aux))
      }
    }
  }
  ret.aux <- ret[!ret$sample %in% out$sample, ]
  ret <- rbind.fill(as.data.frame(out), as.data.frame(ret.aux))
  ret$sample.aux <- NULL
  ret <- ret[match(barcode, ret$barcode), ]
  rownames(ret) <- ret$barcode
  return(ret)
}
environment(colDataPrepare_new) <- asNamespace("TCGAbiolinks")
assignInNamespace("colDataPrepare", colDataPrepare_new, "TCGAbiolinks")
rm(colDataPrepare_new)

get.GRCh.bioMart_new <- function (genome = "hg19", use.gtf = T, as.granges = FALSE) 
{
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    gene.location <- tryCatch({
      host <- ifelse(genome == "hg19", "grch37.ensembl.org", 
                     "www.ensembl.org")
      ensembl <- tryCatch({
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                   host = host)
      }, error = function(e) {
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest", host = host)
      })
      attributes <- c("chromosome_name", "start_position", 
                      "end_position", "strand", "ensembl_gene_id", 
                      "entrezgene_id", "external_gene_name")
      db.datasets <- listDatasets(ensembl)
      description <- db.datasets[db.datasets$dataset == 
                                   "hsapiens_gene_ensembl", ]$description
      message(paste0("Downloading genome information (try:", 
                     tries, ") Using: ", description))
      filename <- paste0(gsub("[[:punct:]]| ", "_", description), 
                         ".rda")
      if (!file.exists(filename)) {
        chrom <- c(1:22, "X", "Y")
        gene.location <- getBM(attributes = attributes, 
                               filters = c("chromosome_name"), values = list(chrom), 
                               mart = ensembl)
        ###
        if(use.gtf) {
          gene.location <- build_gene_loc(gene.location = gene.location)
        }
        ###
        save(gene.location, file = filename)
      }
      else {
        message("Loading from disk")
        gene.location <- get(load(filename))
      }
      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
    })
    if (!is.null(gene.location)) 
      break
  }
  if (tries == 3L) 
    stop("failed to get URL after 3 tries:", "\n  error: ", 
         msg)
  gene.location <- gene.location[! duplicated(gene.location$ensembl_gene_id), ]  # XXX
  if (as.granges) {
    gene.location$strand[gene.location$strand == 1] <- "+"
    gene.location$strand[gene.location$strand == -1] <- "-"
    gene.location$chromosome_name <- paste0("chr", gene.location$chromosome_name)
    gene.location <- makeGRangesFromDataFrame(gene.location, 
                                              seqnames.field = "chromosome_name", start.field = "start_position", 
                                              end.field = "end_position", keep.extra.columns = TRUE)
  }
  return(gene.location)
}
environment(get.GRCh.bioMart_new) <- asNamespace("TCGAbiolinks")
assignInNamespace("get.GRCh.bioMart", get.GRCh.bioMart_new, "TCGAbiolinks")
rm(get.GRCh.bioMart_new)

build_gene_loc <- function(genome_dir = "../../Genomes/human_v22", gene.location) {
  gtf_file <- list.files(path = genome_dir, pattern = "gencode.v[0-9]+.primary_assembly.annotation.gtf", full.names = T)
  cat("[info] use the GTF file:", gtf_file, "\n")
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file, format = "gtf")
  all_genes <- GenomicFeatures::genes(txdb)
  all_genes <- as.data.frame(all_genes)
  rownames(all_genes) <- NULL
  all_genes$ensembl_gene_id <- gsub("\\..*", "", all_genes$gene_id)
  # add gene name
  id_file <- list.files(path = genome_dir, pattern = "gene_ID2Name_fixed.txt", full.names = T)
  gene_name <- read.table(id_file, header = F, sep = "\t", stringsAsFactors = F)
  colnames(gene_name) <- c("id", "name")
  all_genes <- merge(all_genes, gene_name, by.x = "gene_id", by.y = "id", sort = F)
  # add Entrez ID
  gene_entrez <- gene.location[, c("ensembl_gene_id", "entrezgene_id")]
  gene_entrez <- gene_entrez[order(gene_entrez$ensembl_gene_id, gene_entrez$entrezgene_id), ]
  gene_entrez <- gene_entrez[! duplicated(gene_entrez$ensembl_gene_id), ]
  all_genes <- merge(all_genes, gene_entrez, by = "ensembl_gene_id", all.x = T, sort = F)
  # format res
  gene_DF <- all_genes[, c("seqnames", "start", "end", "strand", "ensembl_gene_id", "entrezgene_id", "name")]
  colnames(gene_DF) <- colnames(gene.location)
  gene_DF$chromosome_name <- gsub("^chr", "", gene_DF$chromosome_name)
  gene_DF$strand <- as.integer(ifelse(gene_DF$strand == "+", 1, -1))
  gene_DF <- gene_DF[order(gene_DF$chromosome_name, gene_DF$start_position, gene_DF$end_position), ]
  rownames(gene_DF) <- NULL
  return(gene_DF)
}
environment(build_gene_loc) <- asNamespace("TCGAbiolinks")
