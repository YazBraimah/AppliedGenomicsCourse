### Load required packages
library(edgeR) 
library(tidyverse)
library(pheatmap)
library(Glimma)

library(clusterProfiler)
library(org.Dm.eg.db)

###------------------------###
###------------------------###

setwd("~/Dropbox/Teaching/Syracuse_University/Applied_Genomics/R_data/")

### Load count data from the FlyAtlas2 project
fa_counts = read.csv("fa_counts.txt", header = T, sep = "\t", row.names = 1)

###------------------------###
###------------------------###

### Select columns from the count table for analysis, e.g. Hindgut tissue:

colnames(fa_counts) %>% str_subset(pattern = "Hindgut|Midgut") -> newCOls
fa_counts_guts <- dplyr::select(fa_counts, all_of(newCOls))

###------------------------###
###------------------------###

### Before starting the differential expression analysis, we need to define replicate groups (https://stringr.tidyverse.org/articles/regular-expressions.html)
groups <- colnames(fa_counts_guts) %>% str_replace("_[:digit:]", "")

### Now start the DGEList object
y <- DGEList(counts = fa_counts_guts, group = groups)

### Filter out genes with low expression
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

### Next, we need to normalize by first calculating normalization factors
y <- calcNormFactors(y)

### We have to set up a design matrix to define groups that will be compared:
design <- model.matrix(~ 0 + groups)
colnames(design) %>%  str_replace("groups", "") -> colnames(design)

### Next, estimate the dispersion parameters:
y <- estimateDisp(y, design)

###------------------------###
###------------------------###

### Let's examine features of the data, such as replicate groupings and dispersion estimates:

## Plot sample correlation
data = log2(y$counts+1)
data = as.matrix(data) # required for base R functions
sample_cor = cor(data)
pheatmap(sample_cor)

### Multi dimensional scaling plot
glMDSPlot(y, groups = y$samples$group)

### Biological coefficient of variation
plotBCV(y)

###------------------------###
###------------------------###

### Now fit a quasi-likelihood GLM to the data:
fit <- glmQLFit(y, design)

### To pair-wise differential expression contrasts, we can set up contrasts:
fH.v.mH <- makeContrasts(Female_Hindgut-Male_Hindgut, levels=design)
fH.v.fM <- makeContrasts(Female_Hindgut-Female_Midgut, levels=design)

### Now perform QL F-test:
qlf.fH.v.mH <- glmQLFTest(fit, contrast = fH.v.mH)
qlf.fH.v.mH.tTags <- topTags(qlf.fH.v.mH, n = NULL) # Need to speicy n as NULL tp get all genes
qlf.fH.v.mH.tTags.table <- qlf.fH.v.mH.tTags$table

qlf.fH.v.fM <- glmQLFTest(fit, contrast = fH.v.fM)
qlf.fH.v.fM.tTags <- topTags(qlf.fH.v.fM, n = NULL)
qlf.fH.v.fM.tTags.table <- qlf.fH.v.fM.tTags$table

### Check the number of DE genes:
summary(decideTests(qlf.fH.v.mH))
summary(decideTests(qlf.fH.v.fM))

### get rid of row names:
qlf.fH.v.mH.tTags.table <- qlf.fH.v.mH.tTags.table %>% rownames_to_column(var = "gene") %>% as_tibble()

qlf.fH.v.fM.tTags.table <- qlf.fH.v.fM.tTags.table %>% rownames_to_column(var = "gene") %>% as_tibble()

### CONTINUE ####

### Add a column indicating whether a given gene is significant:
qlf.fH.v.mH.tTags.table %>% 
  mutate(significant = if_else((logFC > 1 & FDR < 0.05) | (logFC < -1 & FDR < 0.05), "yes", "no")) -> qlf.fH.v.mH.tTags.table

qlf.fH.v.fM.tTags.table %>% 
  mutate(significant = if_else((logFC > 1 & FDR < 0.05) | (logFC < -1 & FDR < 0.05), "yes", "no")) -> qlf.fH.v.fM.tTags.table

### Summarize the DE data with a volcano plot:
qlf.fH.v.mH.tTags.table %>% 
  ggplot(aes(x = logFC, y = -log10(PValue), colour = significant)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = filter(qlf.fH.v.mH.tTags.table, -log10(PValue) > 9), aes(label = gene), size = 8) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "orange"))

###------------------------###
###------------------------###

### We can now explore the data in more detail by looking at individual genes that are significant

### Load TPM (transcripts per million) table
fa_tpm = read.csv("fa_tpm.txt", header = T, sep = "\t", row.names = 1)
fa_tpm_guts <- dplyr::select(fa_tpm, all_of(newCOls))

fa_tpm_guts <- fa_tpm_guts %>% rownames_to_column(var = "gene") %>% as_tibble()

fa_tpm_guts %>% gather(key = sample, value = TPM, -gene) %>% separate(sample, c("sample", "gut", "replicate")) -> fa_tpm_guts_gathered

### For heatmaps, we often have to perform some transformations of the TPM values

# make a vector of significant genes
qlf.fH.v.mH.tTags.table %>%  filter(significant == "yes") %>% pull(gene) -> sigGenes_fH.v.mH
# select significant genes only
data = subset(fa_tpm, rownames(fa_tpm) %in% sigGenes_fH.v.mH)
# log transform the TPM values
data = log2(data + 1)
# Calculate the median centered TPM value
data = as.data.frame(t(scale(t(data), scale = F)))
# limit the maximum/minimum log2 value displayed to 4
data[data < -2] = -2
data[data > 2] = 2

## Plot it
pheatmap(mat = data, 
         border_color = NA, 
         show_colnames = TRUE, 
         show_rownames = F, 
         drop_levels = TRUE, 
         annotation_names_row = F, 
         fontsize = 8)

### Single gene plot:
fa_tpm_guts_gathered %>% 
  filter(gene == "FBgn0032171") %>% 
  group_by(gene, sample, gut) %>% 
  summarize(mean = mean(TPM), n = n(), se = sd(TPM)/sqrt(n)) %>% 
  ggplot(aes(fill = gut)) +
    geom_col(aes(sample, mean), position = "dodge") +
    geom_errorbar(aes(x = sample, ymin = mean - se, ymax = mean + se), position = position_dodge(0.9), width = 0.3) +
    geom_point(data = filter(fa_tpm_guts_gathered, gene == "FBgn0032171"), aes(sample, TPM), position = position_dodge2(width = 0.5)) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


### ... If we want to explore expression of many genes, the above code would be too cumbersome. We should create a function:
plotGene <- function(id) {
  fa_tpm_guts_gathered %>% 
    filter(gene == id) %>% 
    group_by(gene, sample, gut) %>% 
    summarize(mean = mean(TPM), n = n(), se = sd(TPM)/sqrt(n)) %>% 
    ggplot(aes(fill = gut)) +
    geom_col(aes(sample, mean), position = "dodge") +
    geom_errorbar(aes(x = sample, ymin = mean - se, ymax = mean + se), position = position_dodge(0.9), width = 0.3) +
    geom_point(data = filter(fa_tpm_guts_gathered, gene == id), aes(sample, TPM), position = position_dodge2(width = 0.5)) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> myPlot ### must make a plot object for a function to work, then use the "return command"
  return(myPlot)
}

## now simply plot with:
plotGene("FBgn0032171")


### For functional enrichment analysis, we often require universal forms of the gene IDs. 

aux_IDs = read.csv("fbgn_auxiliary_ids.txt", header = T, sep ="\t", na.string = "")
annotation = read.csv("fbgn_annotation_ID_fb_2018_04_mod.tsv", header = T, sep ="\t", na.string = "")

### Create a vector of gene names for enrichment analysis:
as_tibble(qlf.fH.v.fM.tTags.table) %>% filter(significant == "yes") %>% pull(gene) -> sigGenes_fH.v.fM

### For GO enrichment analysis, we need to extract the test set and the background set
as_tibble(aux_IDs) %>% filter(primary_FBgn %in% sigGenes_fH.v.fM & na_based_protein_accession != "NA") %>% pull(na_based_protein_accession) -> sigGenes_ncbi

as_tibble(aux_IDs) %>% filter(na_based_protein_accession != "NA") %>% pull(na_based_protein_accession) -> all_ncbi


## Now run GO enrichment:
ego <- enrichGO(gene          = sigGenes_ncbi,
                universe      = all_ncbi,
                # keyType       = 'ENTREZID',
                OrgDb         = org.Dm.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, 
                readable      = TRUE)

ego %>% as_tibble()

dotplot(ego)

### For KEGG enrichment analysis, we can also use the clusterProfile package, but this time use the FLyBase CG ID's. We can check the "organism" code first:

CG_IDs = filter(annotation, primary_FBgn %in% sigGenes_fH.v.fM)
CG_IDs = paste(CG_IDs$organism_abbreviation, CG_IDs$annotation_ID, sep = "_")
my_kegg <- enrichKEGG(gene = CG_IDs, organism = 'dme', qvalueCutoff = 0.05)

### Make a custom ggplot of the KEGG data
filter(my_kegg@result, qvalue < 0.05) %>% as_tibble() %>% 
  separate(GeneRatio, c("geneDE", "background")) %>% 
  mutate(GeneRatio = as.numeric(geneDE)/as.numeric(background)) %>% 
  ggplot(aes(reorder(Description, -log10(p.adjust)), -log10(p.adjust), size = Count, colour = GeneRatio)) +
  geom_point() +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank())


