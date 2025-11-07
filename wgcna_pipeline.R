library(WGCNA)
allowWGCNAThreads()     # enable multithreading
options(stringsAsFactors = FALSE)

setwd("C:/Users/Shreemaran V/microarray_normalization")

# Read expression data
newdata <- read.csv("normalized_expression_matrix.csv", check.names = FALSE, stringsAsFactors = FALSE)
cat("Read expression file: ", dim(newdata), " (rows,cols)\n")

# Identify probe ID column (assume first column is probe IDs)
probe.col <- colnames(newdata)[1]
cat("Using probe ID column:", probe.col, "\n")

# Create datExpr0: rows = samples, cols = probes 
# Remove probe ID column and transpose once
expr_matrix <- as.data.frame(newdata[ , -1, drop = FALSE])  # samples are columns here
# After read.csv, columns 2..n are sample columns; we want rows = samples -> transpose
datExpr0 <- as.data.frame(t(expr_matrix))
# set names (probe IDs) and sample rownames
colnames(datExpr0) <- newdata[[probe.col]]
rownames(datExpr0) <- colnames(newdata)[-1]

# Make sure values are numeric
datExpr0[] <- lapply(datExpr0, function(x) as.numeric(as.character(x)))
if (anyNA(datExpr0)) cat("Note: some NAs present in datExpr0 after coercion.\n")

cat("datExpr0 dimension (samples x genes):", dim(datExpr0), "\n")

# Quality check 
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
cat("goodSamplesGenes allOK:", gsg$allOK, "\n")
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    cat("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "), "\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "), "\n")
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering (outlier detection) 
sampleTree <- hclust(dist(datExpr0), method = "average")
dev.off(ifelse(dev.cur()>1, dev.cur(), 1))  # close non-default device if any
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

cat("datExpr0 dims:", dim(datExpr0), "original newdata dims:", dim(newdata), "\n")

# Soft-threshold power selection 
options(stringsAsFactors = FALSE)
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# plotting
png("softThreshold_plot.png", width = 900, height = 450)
par(mfrow = c(1,2)); cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold(power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()
cat("Soft-threshold diagnostic saved to softThreshold_plot.png\n")

softPower <- 14
cat("Using softPower =", softPower, "\n")

# Build network 
net = blockwiseModules(datExpr0, power = softPower, TOMType = "unsigned",
                       minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "newdata-resistTom", verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

# Plot dendrogram and colors 
png("dendro_colors.png", width = 1200, height = 800)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
cat("Dendrogram+colors saved to dendro_colors.png\n")

moduleColors = labels2colors(net$colors)
moduleLabels = net$colors
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
cat("nSamples=", nSamples, "nGenes=", nGenes, "\n")

# Module eigengenes 
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Read trait data 
datTraits <- read_excel("trait_data.xlsx")
cat("Trait data:", dim(datTraits), "\n")

# Ensure SampleID column exists
if (!"SampleID" %in% colnames(datTraits)) {
  stop("trait_data.xlsx must have a column named 'SampleID'")
}

# Set rownames and remove SampleID column
rownames(datTraits) <- datTraits$SampleID
datTraits <- datTraits[ , !(colnames(datTraits) %in% "SampleID"), drop = FALSE]

# Align sample order with expression data
commonSamples <- intersect(rownames(datExpr0), rownames(datTraits))
if (length(commonSamples) == 0) stop("No matching sample names between expression and trait data")

datExpr0 <- datExpr0[commonSamples, , drop = FALSE]
datTraits <- datTraits[commonSamples, , drop = FALSE]
cat("Aligned samples:", nrow(datExpr0), "\n")

# Module–trait correlation 
nSamples <- nrow(datExpr0)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png("moduleTraitHeatmap.png", width = 1200, height = 800)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = "Module–Trait Relationships")
dev.off()
cat("✅ Module–trait heatmap saved to moduleTraitHeatmap.png\n")

# Gene significance and module membership 
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Correlate each gene with each trait
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

# Combine all results
probes <- colnames(datExpr0)
geneInfo <- data.frame(Probe = probes,
                       moduleColor = moduleColors[match(probes, probes)],
                       geneModuleMembership,
                       MMPvalue,
                       geneTraitSignificance,
                       GSPvalue,
                       stringsAsFactors = FALSE)

write.csv(geneInfo, "geneInfo.csv", row.names = FALSE)
cat("✅ geneInfo.csv saved successfully\n")

# Save session
save(MEs, moduleLabels, moduleColors, geneTree, datExpr0, datTraits, geneInfo, file = "final_wgcna_results.RData")
cat("✅ final_wgcna_results.RData saved successfully\n")
