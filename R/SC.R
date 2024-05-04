# (0) Libraries #####
library(edgeR)
library(tidyverse)
library(GenomicFeatures)
library(ggbiplot)

cutoff <- c(-1,1)

# Session > Set working directory > Source file location

# (1) Counts ####

sample_info <- read.delim("Counts_C.txt")
targets <- data.frame(files = sample_info$files, group = paste0(sample_info$group),
                      description = paste0(sample_info$replicate))
row.names(targets) <- targets$description
d <- readDGE(targets, header = FALSE) 

keep <- filterByExpr(d)
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d)
plotMDS(d, col=as.numeric(d$samples$group)) # plot leading fold changes - to get an idea if data have batch effects

# Uncomment the next 3 lines for exporting the MDS plot
#png(filename = "MDSplot.png", res = 300, width = 1920, height = 1440)
#plotMDS(d, col=rep(1:length(targets[,2]), each=(length(targets[,2])/length(unique(targets[,2])))))
#dev.off()


# (2) Model fitting ####

## WITHOUT BATCH EFFECT CORRECTION

# Create design matrix
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
my.contrasts <- makeContrasts(
  VC_VCP = VC-VCP,
  V_VCP = V-VCP,
  # include more comparisons here if necessary
  levels=design)

# Estimating the dispersions and plot them
d <- estimateDisp(d, design, robust=TRUE) ## does both dispersion estimations in one step - suggested in edgeR manual;
plotBCV(d) # plots the results of the dispersion estimation

# Uncomment the next 3 lines to export the dispersion plot
#pdf("disp.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotBCV(d)
#dev.off()# closes and saves the pdf file with the plot

fit <- glmQLFit(d, design, robust = TRUE) # fits a model onto the data
plotQLDisp(fit) # plots the result of the model fitting

qlf <- glmQLFTest(fit, contrast=my.contrasts[,c("VC_VCP")]) # differential expression analysis with the comparison (=contrast) specified in quotes (multiple comparisons can be done separated by a comma)
topTags(qlf) # shows the top differentially expressed genes
summary(decideTests(qlf)) # gives a summary of the differential expression analysis

plotMD(qlf) # plots the results
abline(h=c(-1, 1), col="blue")

# Uncomment the next 4 lines to export the results plot
#pdf("diff_expr_results.pdf")# change the name of the .pdf file; this opens a pdf file to store the next plot
#plotMD(qlf)
#abline(h=cutoff, col="blue")
#dev.off()

diff_results <- as.data.frame(topTags(qlf, n=Inf)) # saves the results in the variable diff_results

# Export results
write_tsv(diff_results, "diff_results.tsv")

gff_file <- "sc.gff"
GeneData <- rtracklayer::import.gff(gff_file)
SC_Entrez <- as.data.frame(elementMetadata(GeneData[GeneData$type == "gene", c(5,6,7)]))
SC_Entrez$Dbxref <- gsub(pattern = "GeneID:", replacement = "", x = SC_Entrez$Dbxref)
colnames(SC_Entrez) <- c("locus_tag", "EntrezID", "GeneName")

diff_results2 <- merge(diff_results, SC_Entrez, by.x=0, by.y=1, all.x=TRUE)

# Export results with gene names
write_tsv(diff_results2, "diff_results_wNames.tsv")



##  WITH BATCH EFFECT CORRECTION
# TBA
Ratio <- sample_info$group
Ratio <- as.factor(Ratio)
Ratio <- relevel(Ratio, ref = "V")
Time <- sample_info$time
Repl <- sample_info$replicate

# (3) PCA ####
cpmill <- cpm(d$counts, normalized.lib.size=TRUE)
cpmill_transp <- data.frame(t(cpmill))
colnames(cpmill_transp) <- row.names(cpmill)
cpmill_transp$group <- as.factor(targets$group)
cpmill_transp$time <- as.factor(targets$time)

vars <- rownames(cpmill_transp)
cpmill_transp.pca <- prcomp(dplyr::select(cpmill_transp,-group,-time), scale. = TRUE)
plot <- ggbiplot(cpmill_transp.pca, choices = c(1,2), var.axes = FALSE, groups = cpmill_transp$group:cpmill_transp$time, colour = cpmill_transp$group:cpmill_transp$time)
plot + 
  geom_point(aes(color = cpmill_transp$group:cpmill_transp$time),size = 2)+
  scale_color_manual(values = c("V:6" = "#44AA99", "V:24" = "#88CCEE","VC:6" = "#882255","VC:24" = "#CC6677","VCP:6" = "#699729","VCP:24" = "#67E117","VP:6" = "#B03833","VP:24" = "#EF9EAC"))


# Uncomment these 3 lines to export plot
#pdf(file = "PCA.pdf")
#plot
#dev.off()
# (4) Gene Ontology ####
up <- diff_results2[diff_results2$FDR < 0.05 & diff_results2$logFC > 0,] # change logFC if necessary
names <- up[,8]
writeLines(upnames, "GO_CA.txt")

# Use gene names from text file on ShinyGO:
## http://bioinformatics.sdstate.edu/go/
# (5) Scatterplot ####
### Enter axis labels ###
xlb = "Average counts per million reads (log2)"
ylb = "Fold change (log2)"

### Plotting
plot_data <- diff_results2
plot_data$Group <- ifelse((plot_data$FDR < 0.05), 'diff_reg', 'not_reg')

myplot <- ggplot(data=plot_data, aes(x=logCPM, y=logFC)) +
  geom_point() +
  geom_point(data = plot_data[plot_data$Group == 'diff_reg',], color = 'red') +
  geom_hline(yintercept = c(-1,1), color='blue', linetype="dashed") +
  #xlim(c(-5,15)) +
  #ylim(c(-10,15)) +
  xlab(xlb) +
  ylab(ylb) +
  cleanup
myplot

myplot_lab <- myplot + geom_text(aes(label=ifelse(((logFC >= 5 & FDR < 0.05) | 
                                                     (logFC <= -2.5 & FDR < 0.05)),
                                                  as.character(GeneName),'')),
                                 size = 3  ,hjust=-0,vjust=-0.5, colour="black")
myplot_lab

# Uncomment next 3 lines to export plot
#pdf('MyPlot.pdf')
#myplot_lab
#dev.off()

# (6) Volcano plot ####

xlb_volc = "log2 fold change hir1 vs wt"
ylb_volc = "-log10(adjusted p-value)"

volcplot <- ggplot(data=plot_data, aes(x=logFC, y=-log10(FDR))) +
  geom_point()+
  geom_point(data = plot_data[plot_data$Group == 'diff_reg' & plot_data$logFC >= 0,], color = 'darkred') +
  geom_point(data = plot_data[plot_data$Group == 'diff_reg' & plot_data$logFC < 0,], color = 'darkblue') +
  geom_vline(xintercept = c(-1,1), color='black', linetype="dashed",size=0.4) +
  geom_vline(xintercept = 0, color='black', linetype="dashed",size=0.6) +
  geom_hline(yintercept = -log10(0.05), color='black', linetype="dashed",size=0.4) +
  #xlim(c(-5,15)) +
  #ylim(c(0,130)) +
  xlab(xlb_volc) +
  ylab(ylb_volc)
volcplot

#volcplot_lab <- volcplot + geom_text(aes(label=ifelse(((logFC >= 2 & FDR < 0.05) | (logFC <= -2 & FDR < 0.05)),
               #                                       as.character(GeneName),'')),
               #                      size = 3  ,hjust=-0,vjust=-0.5, colour="black")
#volcplot_lab

# Uncomment next 3 lines to export plot
#pdf('MyPlot.pdf')
#volcplot_lab
#dev.off()