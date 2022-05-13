

#read your csv file
diff_exp <- read.csv ("analysis5.6.csv", col.names = c('PROBEID', 't', 'p', 'p_adjust'))

#Using the select() function of the bioconductor package hgu133plus2.db, map the probeset IDs to gene symbols by specifying the appropriate key and column arguments. Some probeset IDs map to the same gene symbol, so reason about and pick a rationale for choosing which probeset ID to use as representative. Add an additional column to the differential expression results that contains one symbol for each probeset ID.
gene_diff <- AnnotationDbi::select(hgu133plus2.db, diff_exp$PROBEID, c('SYMBOL'))


#collapser is a function that combines genesymbols that map to 2 different probeids using '|' 
collapser <- function(a) {a %>% unique %>% sort %>% paste(collapse = '|') }
#using collapser, we combine multiple genesymbols with same probeid
gene_diff <- gene_diff %>% group_by(PROBEID) %>% summarise_each(funs(collapser)) %>% ungroup

#storing and merging biological results with genesymbol 
merge_symbol <- merge(biological_results, genesymbol, on = 'PROBEID')
View(merge_symbol)

#remove duplicated rows
merge_symbol <- merge_symbol[!duplicated(merge_symbol$SYMBOL),]

#removes blank rows
merge_symbol <- merge_symbol[!(is.na(merge_symbol$SYMBOL) | merge_symbol$SYMBOL == ""), ]

#load gene sets 
hallmarks <- getGmt('h.all.v7.5.1.symbols.gmt')
GO <- getGmt('c5.go.v7.5.1.symbols.gmt')
KEGG <- getGmt('KEGG.gmt')

#get geneset length

print(paste0("There are ", length(hallmarks)," gene sets in the Hallmark collection"))
print(paste0("There are ", length(KEGG)," gene sets in the Kegg collection"))
print(paste0("There are ", length(GO)," gene sets in the GO collection"))      

#order the values in decreasing values of t-stats 
merge_symbol <- merge_symbol[order(merge_symbol$t, decreasing = TRUE), ]


#top 1000 up regulated and down regulated
up_1000 <- head(merge_symbol, n = 1000)
down_1000 <- tail(merge_symbol, n = 1000)

#selecting top 10 for report
up_10 <- head(up_1000, n = 10)
down_10 <- tail(down_1000, n = 10)

#store the top10 up and down regulated genes
write.csv(up_10, "10_upregulated_genes.csv")
write.csv(down_10, "10_downregulated_genes.csv")

#extract the genes that were not expressed
not_diffexp_up <- subset(merge_symbol, !merge_symbol$SYMBOL %in% up_1000$SYMBOL)
not_diffexp_down <- subset(merge_symbol, !merge_symbol$SYMBOL %in% down_1000$SYMBOL)

#define function to create contigency table
fishertest <- function(gl, gs, nde)           
{ diffexp_ings <- length(intersect(gl,gs))    
diffexp_notgs <- length(gl) - diffexp_ings    
notde_ings <- length(intersect(nde,gs))       
notde_notgs <- length(nde) - notde_ings      
return(c(diffexp_ings,diffexp_notgs,notde_ings,notde_notgs))}   #returns the fishertest values

#stores results of fisher test for hallmark geneset 
hallmarks_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

#stores the results for hallmark geneset comparison in separate data frame using for loop
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}

hallmarks_results <- hallmarks_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))


#stores results of fisher test for kegg geneset
kegg_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

##stores the results for kegg geneset comparison in separate data frame using for loop
for (i in 1:length(KEGG))
{
  geneid <- geneIds(KEGG[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}

kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))


#stores results of fisher test for kegg geneset
go_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

##stores the results for go geneset comparison in separate data frame using for loop
for (i in 1:length(GO))
{
  geneid <- geneIds(GO[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  go_results[nrow(go_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  go_results[nrow(go_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}

go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))


#adjusting the pvalue usinf benjamini hochberg method and storing the data in seperate file
go_results$BH <- p.adjust(go_results$pvalue, method = "BH", n = length(go_results$pvalue))
write.csv(go_results, "final_go.csv")     

kegg_results$BH <- p.adjust(kegg_results$pvalue, method = "BH", n = length(kegg_results$pvalue))
write.csv(kegg_results, "final_kegg.csv")

hallmarks_results$BH <- p.adjust(hallmarks_results$pvalue, method = "BH", n = length(hallmarks_results$pvalue))
write.csv(hallmarks_results, "final_hallmarks.csv")


top3_kegg <- slice_min(kegg_results, order_by=pvalue, n=3)
top3_go <- slice_min(go_results, order_by=pvalue, n=3)
top3_hm <- slice_min(hallmarks_results, order_by=pvalue, n=3)
top3s <- rbind(top3_kegg, top3_go, top3_hm)
print(top3s)

write.csv(top3s, file="geneset_top_results.csv")

fisher_pvals <- function(gmt){
  pvalues <- c()
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(merge_symbol$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% merge_symbol$SYMBOL)
    in_set_not_differential <- in_set - in_set_differential
    not_in_set_differential <- differentially_expressed - in_set_differential
    not_in_set_not_differential <- 0
    fishervals <- fisher.test(matrix(c(in_set_differential, in_set_not_differential, not_in_set_differential, not_in_set_not_differential), nrow = 2))
    pval <- fishervals$p.value
    pvalues[setname] <- pval
  }
  return(pvalues)
}

# Creating function to generate dataframe of fisher test values
fisher_df <- function(gmt){
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(merge_symbol$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% merge_symbol$SYMBOL)
    in_set_not_differential <- in_set - in_set_differential
    not_in_set_differential <- differentially_expressed - in_set_differential
    not_in_set_not_differential <- 0
    fishervals <- fisher.test(matrix(c(in_set_differential, in_set_not_differential, not_in_set_differential, not_in_set_not_differential), nrow = 2))
    pval <- fishervals$p.value
    est <- fishervals$estimate
    padj <- p.adjust(pval, method="fdr")
    df[[setname]] <- data.frame(geneset = setname, statistic = est, pvalue = pval, p.adj = padj)
  }
  return(df)
}

#KEGG pathway
pvalues_kegg <- fisher_pvals(KEGG)
df_kegg <- fisher_df(KEGG)

# GO Pathways
pvalues_go <- fisher_pvals(GO)
df_go <- fisher_df(GO)

# Hallmark Pathways
pvalues_hallmarks <- fisher_pvals(hallmarks)
df_hallmarks <- fisher_df(hallmarks)
