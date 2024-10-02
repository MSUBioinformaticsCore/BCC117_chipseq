# findEnrichedGOBP --------------------------------------------------------
#' @name findEnrichedGOBP
#' @description uses topGO to find enriched GOBP terms for a list of genes.
#' ------------ Pvalues from weight01 fisher test
#' @param goi character vector of genes of interest
#' @param background character vector of all background genes
#' @param org.db the org.db for the organism of interest from bioconductor. 
#' ------------- Must already be installed.
#' @param id_type gene identifier to use. 
#' ------------- c("entrez","genbank","alias","ensembl","symbol","genename")
#' @param min_size smallest number of genes for a GO term to be used
#' @param max_size largest number of genes for a GO term to be used
#' @return a data frame with the pvalue and overlapping genes 
#' ------- between the genes of interest and GOBPs. 

findEnrichedGOBP = function(goi, 
                            background, 
                            org.db, 
                            id_type,
                            min_size,
                            max_size){
  
  require(paste(org.db), character.only = TRUE)
  require(tidyverse)
  require(topGO)
  
  submit = factor(ifelse(background %in% goi, 1, 0))
  names(submit) = background
  
  GO = new("topGOdata",
           ontology = "BP",
           allGenes = submit,
           nodeSize = min_size,
           annot = annFUN.org,
           mapping = org.db,
           ID = id_type)
  
  nTerms = length(usedGO(GO))
  
  Fisher <- runTest(GO,
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  results = GenTable(GO, 
                     pval = Fisher, 
                     topNodes = nTerms,
                     numChar=1000)
  
  results = results %>%
    filter(Annotated <=  max_size)
  
  results$FDR = p.adjust(results$pval, method = "BH")
  
  ann.genes = genesInTerm(GO, results$GO.ID)
  symbol = names(submit)[submit == 1]
  
  results$OverlappingGenes = sapply(results$GO.ID, 
                                    function(x){
                                      sig = ann.genes[[x]][ann.genes[[x]] %in% symbol]
                                      sig = paste(sig, collapse = ", ")
                                      return(sig)}
  )
  return(results)
}

# findEnrichedGOMF --------------------------------------------------------
#' @name findEnrichedGOMF
#' @description uses topGO to find enriched GOMF terms for a list of genes.
#' ------------ Pvalues from weight01 fisher test
#' @param goi character vector of genes of interest
#' @param background character vector of all background genes
#' @param org.db the org.db for the organism of interest from bioconductor. 
#' ------------- Must already be installed.
#' @param id_type gene identifier to use. 
#' ------------- c("entrez","genbank","alias","ensembl","symbol","genename")
#' @param min_size smallest number of genes for a GO term to be used
#' @param max_size largest number of genes for a GO term to be used
#' @return a data frame with the pvalue and overlapping genes 
#' ------- between the genes of interest and GOBPs. 

findEnrichedGOMF = function(goi, 
                            background, 
                            org.db, 
                            id_type,
                            min_size,
                            max_size){
  
  require(paste(org.db), character.only = TRUE)
  require(tidyverse)
  require(topGO)
  
  submit = factor(ifelse(background %in% goi, 1, 0))
  names(submit) = background
  
  GO = new("topGOdata",
           ontology = "MF",
           allGenes = submit,
           nodeSize = min_size,
           annot = annFUN.org,
           mapping = org.db,
           ID = id_type)
  
  nTerms = length(usedGO(GO))
  
  Fisher <- runTest(GO,
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  results = GenTable(GO, 
                     pval = Fisher, 
                     topNodes = nTerms,
                     numChar=1000)
  
  results = results %>%
    filter(Annotated <=  max_size)
  
  results$FDR = p.adjust(results$pval, method = "BH")
  
  ann.genes = genesInTerm(GO, results$GO.ID)
  symbol = names(submit)[submit == 1]
  
  results$OverlappingGenes = sapply(results$GO.ID, 
                                    function(x){
                                      sig = ann.genes[[x]][ann.genes[[x]] %in% symbol]
                                      sig = paste(sig, collapse = ", ")
                                      return(sig)}
  )
  return(results)
}

# findEnrichedGOCC --------------------------------------------------------
#' @name findEnrichedGOMF
#' @description uses topGO to find enriched GOMF terms for a list of genes.
#' ------------ Pvalues from weight01 fisher test
#' @param goi character vector of genes of interest
#' @param background character vector of all background genes
#' @param org.db the org.db for the organism of interest from bioconductor. 
#' ------------- Must already be installed.
#' @param id_type gene identifier to use. 
#' ------------- c("entrez","genbank","alias","ensembl","symbol","genename")
#' @param min_size smallest number of genes for a GO term to be used
#' @param max_size largest number of genes for a GO term to be used
#' @return a data frame with the pvalue and overlapping genes 
#' ------- between the genes of interest and GOBPs. 

findEnrichedGOCC = function(goi, 
                            background, 
                            org.db, 
                            id_type,
                            min_size,
                            max_size){
  
  require(paste(org.db), character.only = TRUE)
  require(tidyverse)
  require(topGO)
  
  submit = factor(ifelse(background %in% goi, 1, 0))
  names(submit) = background
  
  GO = new("topGOdata",
           ontology = "CC",
           allGenes = submit,
           nodeSize = min_size,
           annot = annFUN.org,
           mapping = org.db,
           ID = id_type)
  
  nTerms = length(usedGO(GO))
  
  Fisher <- runTest(GO,
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  results = GenTable(GO, 
                     pval = Fisher, 
                     topNodes = nTerms,
                     numChar=1000)
  
  results = results %>%
    filter(Annotated <=  max_size)
  
  results$FDR = p.adjust(results$pval, method = "BH")
  
  ann.genes = genesInTerm(GO, results$GO.ID)
  symbol = names(submit)[submit == 1]
  
  results$OverlappingGenes = sapply(results$GO.ID, 
                                    function(x){
                                      sig = ann.genes[[x]][ann.genes[[x]] %in% symbol]
                                      sig = paste(sig, collapse = ", ")
                                      return(sig)}
  )
}

# overlapSets ------------------------------------------------------------
#' @name overlap_sets
#' @description performs pairwise hypergeometric tests between multiple sets 
#' of genes i.e. cluster marker genes from two different single cell
#' different single cell data sets
#' @param table1 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 1      
#' @param table2 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Cluster
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param background the number of genes in the gene universe 
#' ie. all expressed genes
#' @return a list with two elements: 
#' [1]intersecting.genes: a tidy table where each row shows two groups 
#' (one from table1 and one from table2) and an intersecting gene or "none"
#' [2] pval.scores: a tidy table where each row is a pair of groups 
#' (one from table1 and one from table2) with the number of genes in each group,
#' the nunber of intersecting genes, the uncorrected pval, and the enrichment score

overlapSets <- function(table1, table2, background){
  require(dplyr)
  
  # make table1 into a list divided by group
  colnames(table1) = c("Gene","Group")
  tab1_groups = unique(table1$Group)
  tab1_list = list()
  
  for(group in tab1_groups){
    tab1_list[[group]] = table1 %>%
      filter(Group == group)
  }
  
  # make table2 into a list divided by group
  colnames(table2) = c("Gene","Group")
  tab2_groups = unique(table2$Group)
  tab2_list = list()
  
  for(group in tab2_groups){
    tab2_list[[group]] = table2 %>%
      filter(Group == group)
  }
  
  # make a list of dfs with the names of genes shared between
  # every group in table1 and every group in table2
  INT = list()
  # make a list of dfs with the number of genes shared between
  # every group in table1 and every group in table2
  nINT = list()
  
  for(group in tab1_groups){
    genes = tab1_list[[group]]$Gene
    INT_list = lapply(tab2_list, 
                      function(x){
                        y = intersect(x$Gene, genes);
                        if(length(y) > 0){
                          z = data.frame(SharedGene = y)}
                        else{z = data.frame(SharedGene = "none")};
                        z
                      })
    
    for(i in 1:length(INT_list)){
      INT_list[[i]]$Group2 = names(INT_list)[i]
    }
    
    INT[[group]] = do.call(rbind,INT_list)
    INT[[group]]$Group1 = group
    
    INT[[group]] = INT[[group]] %>%
      select(Group1, Group2, SharedGene)
    
    nINT_list = lapply(tab2_list, 
                       function(x){
                         y = intersect(x$Gene, genes);
                         z = length(y);
                         z
                       })
    
    nINT_vec = unlist(nINT_list)
    nGroup2 = unlist(lapply(tab2_list, nrow))
    
    nINT[[group]] = data.frame(Group1 = group,
                               nGroup1 = length(genes),
                               Group2 = names(nINT_vec),
                               nGroup2  = nGroup2,
                               nSharedGenes = nINT_vec)
    
    
  }
  
  INTdf = do.call(rbind, INT)
  nINTdf = do.call(rbind, nINT)
  
  nINTdf = nINTdf %>%
    mutate(Pval = phyper(nSharedGenes-1,
                         nGroup2,
                         background-nGroup2,
                         nGroup1,
                         lower.tail= FALSE)
    )
  
  nINTdf = nINTdf %>%
    mutate(Enrichment = 
             log2(
               (nSharedGenes/nGroup2)/
                 (nGroup1/background))
    )
  
  return(list(intersecting.genes = INTdf, 
              pval.scores = nINTdf))
}  

# fgseaEnrichment ---------------------------------------------------------
#' @name fgseaEnrichment
#' @description performs GSEA enrichments pairwise between a ranked lists of  
#' of genes i.e. cluster marker genes with some stat like -log10FDR and other gene sets
#' @param stat_table a dataframe with three columns from dataset 1, 
#' ie Gene, Set, Stat. Should include all genes with no filter.
#' format example
#' --------------  
#' Gene  Set      Stat
#' Gene1 cluster1 10
#' Gene2 cluster1 100
#' Gene3 cluster2 5
#' Gene4 cluster3 20
#' @param set_table a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Set
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param scoreType "pos" if all stats are positive not signed, "std" otherwise
#' ie. all expressed genes
#' @return  a tidy table where each row is a pair of groups 
#' (one from stat_table and one from set_table) with the number of genes in each group,
#' the number of intersecting genes, the uncorrected pval, and the enrichment score

fgseaEnrichment = function(stat_table, set_table, scoreType) {
  
  require(fgsea)
  require(tidyverse)
  
  set_list = list()
  for(set in unique(set_table$Set)){
    
    set_list[[set]] = 
      set_table %>%
      filter(Set == set) %>%
      pull(Gene)
  } 
  
  stat_list = list()
  for(set in unique(as.character(stat_table$Set))){
    
    stat_list[[set]] = 
      stat_table %>%
      filter(Set == set) %>%
      pull(Stat)
    
    stat_list[[set]] =  
      as.numeric(stat_list[[set]])
    
    max_stat = max(stat_list[[set]][is.finite(stat_list[[set]])])
    stat_list[[set]][is.infinite(stat_list[[set]])] = max_stat + 1
    
    names(stat_list[[set]]) =
      stat_table %>%
      filter(Set == set) %>%
      pull(Gene)
  }
  
  apply_fgseaMultilevel = function(stats, pathways, minSize, scoreType){
    
    res = 
      fgseaMultilevel(
        pathways = pathways,
        stats = stats,
        minSize = minSize,
        scoreType = scoreType
      )
    
    return(res)
  }
  
  res_list = lapply(stat_list, 
                    apply_fgseaMultilevel,
                    pathways = set_list,
                    minSize = 1,
                    scoreType = scoreType)
  
  names(res_list) = names(stat_list)
  
  for(i in 1:length(res_list)){
    
    res_list[[i]]$StatDataset = names(res_list)[i]
    
    res_list[[i]] = 
      res_list[[i]] %>%
      select(-leadingEdge) %>%
      rename(SetDataset = pathway) %>%
      select(StatDataset, everything())
  }
  
  res = do.call(rbind, res_list)
  return(res)
}

# entrez2symbol -----------------------------------------------------------
entrez2symbol = function(entrez_ids){
  require("org.Hs.eg.db")
  x <- org.Hs.egSYMBOL
  included_entrez <- names(as.list(x))
  included_oi = entrez_ids[entrez_ids %in% included_entrez]
  xx <- as.list(x[included_oi])
  gene_symbols = unlist(xx)
  unloadNamespace("org.Hs.eg.db")
  detach("package:AnnotationDbi")
  return(gene_symbols)
}
# symbol2entrez -----------------------------------------------------------
symbol2entrez = function(gene_symbols){
  require("org.Hs.eg.db")
  x <- org.Hs.egSYMBOL2EG
  included_symbols <- names(as.list(x))
  included_oi = gene_symbols[gene_symbols %in% included_symbols]
  xx <- as.list(x[included_oi])
  entrez_ids = unlist(xx)
  unloadNamespace("org.Hs.eg.db")
  detach("package:AnnotationDbi")
  return(entrez_ids)
}