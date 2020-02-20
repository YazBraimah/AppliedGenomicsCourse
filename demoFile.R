
### Load required packages
library(edgeR) 
library(tidyverse)
library(patchwork)

setwd("~/Dropbox/Teaching/Syracuse_University/Applied_Genomics/RNA-seq demo/fastqc")


### Load count data
counts = read.csv("hcc1395_counts.matrix", header = T, sep = "\t", row.names = 1)


### Examine table dimensions:
dim(counts)

## Check correlation between replicates:
n1 <- ggplot(data = counts, aes(x = log2(normal_rep1), y = log2(normal_rep2))) + geom_point(colour = "red", alpha = 0.75)
n2 <- ggplot(data = counts, aes(x = log2(normal_rep1), y = log2(normal_rep3))) + geom_point(colour = "green", alpha = 0.75)
n3 <- ggplot(data = counts, aes(x = log2(normal_rep2), y = log2(normal_rep3))) + geom_point(colour = "blue", alpha = 0.75)
n1 + n2 + n3

### Create Differential Gene Expression List Object (DGEList) object
d0 <- DGEList(counts)


### Derive experiment metadata from the sample names
snames <- colnames(counts) # Sample names
snames

### Create a new variable “group” that combines the replicates into samples
group <- substr(snames, 1, 3)
group

### Calculate normalization factors
d0 <- calcNormFactors(d0)
d0$samples


### Filter genes ith low expression
cutoff <- 3
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left


### Examine the relationship among samples
plotMDS(d, col = as.numeric(group))

### make log CPM object for plotting
logcpm <- cpm(d, prior.count=2, log=TRUE)

### define a matrix for differential expression
mm <- model.matrix(~0 + group)
mm

### Look at mean-variance relationships
y <- voom(d, mm, plot = T)


### lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))

### Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
contr <- makeContrasts(groupnor - grouptum, levels = colnames(coef(fit)))
contr


### Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 30)


### How many genes are differentially expressed?
length(which(top.table$adj.P.Val < 0.05))


### represent the differentiall expressed genes in a volcano plot:
top.table$sig = ifelse(top.table$adj.P.Val < 0.05 & (top.table$logFC > 1 | top.table$logFC < -1), "yes", "no")
ggplot(top.table, aes(logFC, -log10(P.Value), colour = sig)) + geom_point() + geom_vline(xintercept = 0, linetype = "dashed") +
scale_colour_manual(values = c("gray", "purple")) 
