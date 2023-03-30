# Install necessary libraries
#source("http://bioconductor.org/biocLite.R")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("fontconfig")
#BiocManager::install("freetype2")
#BiocManager::install("tidyverse")
#BiocManager::install("dplyr")
#BiocManager::install("RColorBrewer")
#BiocManager::install("ggrepel")

# Load necessary libraries
#library(tidyverse)
library(dplyr)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Set up the correct working directory
setwd("~/Documents/Master/ml/project/TCGA_dataset/TCGA_dataset/downstream")

# Create empty directory "output"
# The directory is for all output of the R script
dir.create("output")


# Load count data into data frame
seqdata <- read.csv("seqdata_tpm.csv", stringsAsFactors = FALSE)

# Load top features
top_features <- read.table("top_features.txt", stringsAsFactors = FALSE)

# Load meta data into data frame
metadata <- read.csv("prediction_file_crc.csv", stringsAsFactors = FALSE)

filtered_meta <- filter(metadata, msi_status == "MSI-H" | msi_status == "MSS" | msi_status == "MSI-L")

msi_h.names <- filter(filtered_meta, msi_status == "MSI-H")$X
mss_msi_l.names <- filter(filtered_meta, msi_status == "MSS" | msi_status == "MSI-L")$X

# Set gene names as row names
countdata <- seqdata[,(-1)]
rownames(countdata) <- seqdata[,(1)]

# Subset of TPM data with only the most important genes (features)
countdata <- countdata[which(rownames(countdata) %in% top_features$V1),]

# Show 6 first lines of raw count data
#head(countdata)

# Show number of genes (rows) and samples (columns) of count data
dim(countdata)

# Change "-" in sample names to "."
msi_h.names <- gsub("\\-", ".", msi_h.names)
mss_msi_l.names <- gsub("\\-", ".", mss_msi_l.names)

# Calculate the mean and standard deviation of each gene for the samples classified as MSI-L and MSS as control
#   while the samples classified as MSI-H as knockdown
control.mean = apply(countdata[mss_msi_l.names],1,mean)
control.sd =  apply(countdata[mss_msi_l.names],1,sd)
knockdown.mean = apply(countdata[msi_h.names],1,mean)
knockdown.sd =  apply(countdata[msi_h.names],1,sd)


# Plot results for control samples in blue
plot(control.mean,control.sd,log = 'xy',pch = 20,cex=.5,col='blue',xlab = 'Mean Rate',ylab = 'Stand dev')
# and the knockdown samples in red
points(knockdown.mean,knockdown.sd,pch=20,cex=0.5,col='red')
title('Standard deviation vs mean rate')

# Calculates the slope to correct for
all.fit = lm(sd~mean, data.frame(mean=log(c(control.mean,knockdown.mean)),sd=log(c(control.sd,knockdown.sd))))
coef(all.fit)

# Transforms the TPM counts with the slope
transformed = as.data.frame(countdata^(1-coef(all.fit)['mean']))
#head(transformed)

#mean and sd for transformed values
control.t.mean = apply(transformed[mss_msi_l.names],1,mean)
control.t.sd = apply(transformed[mss_msi_l.names],1,sd)
knockdown.t.mean = apply(transformed[msi_h.names],1,mean)
knockdown.t.sd = apply(transformed[msi_h.names],1,sd)

# Plot sd vs mean of transformed data
plot(x=control.t.mean, y=control.t.sd, pch=20, cex=0.5, col='blue', ylim = c(-0.03,0.3),
     xlab='Mean transformed value', ylab='Sd of transformed value')
points(x=knockdown.t.mean, y=knockdown.t.sd, pch=20, cex=0.5, col='red')
title('SD versus mean rate of transformed data')

# T-test on all rows, and extracting the p-values from the test result
pvalues.t <- apply(transformed, 1, function(row) 
{t.test(x=row[mss_msi_l.names], y=row[msi_h.names], var.equal=TRUE)$p.value})

# Get the total counts of each sample
totalcounts = colSums(countdata)
head(totalcounts)

# Calculate log ratios of the average counts of kd samples / ct samples
log2ratio <- apply(
  countdata, 1, function(x) {log(sum(x[msi_h.names])/sum(totalcounts[msi_h.names]),2) -
      log(sum(x[mss_msi_l.names])/sum(totalcounts[mss_msi_l.names]),2)})

head(log2ratio)

# combine log2ratios and p-values to create a volcano plot
par(mfrow=c(1,1))
p05 = -log10(0.05)
plot(y=-log10(pvalues.t),x=log2ratio)
lines(x=c(-1,-1),y=c(p05,50),col='gray')
lines(x=c(1,1),y=c(p05,50),col='gray')
lines(x=c(-6,-1),y=c(p05,p05),col='gray')
lines(x=c(1,4),y=c(p05,p05),col='gray')
text(x=1.75, y=p05+2, label="P = 0.05")
text(x=1.15, y=40, label="FC = 2", srt=90)
title('volcano plot')


countdata$pval <- pvalues.t
countdata$log2fc <- log2ratio


# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
countdata$delabel <- ifelse(rownames(countdata) %in% noquote(names(head(pvalues.t[order(pvalues.t)], 30))), noquote(rownames(countdata)), NA)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
countdata$diffexpressed <- "NO"
  # if log2Foldchange > 1 and pvalue < 0.05, set as "UP"<br /><br /><br />
countdata$diffexpressed[countdata$log2fc > 1 & countdata$pval < 0.05] <- "UP"
  # if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"<br /><br /><br />
countdata$diffexpressed[countdata$log2fc < -1 & countdata$pval < 0.05] <- "DOWN"
  
#head(countdata[order(pvalues.t) & countdata$diffexpressed == 'DOWN', ])
  

myvolcanoplot <- ggplot(data = countdata, mapping = aes(x = log2ratio, y = -log10(pvalues.t), col = diffexpressed, label = delabel)) +

geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +

geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +

geom_point(size = 2) +

# to set the colours of our variable
# to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
# since some genes can have minuslog10padj of inf, we set these limits
scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 50), xlim = c(-3, 3)) +
  labs(color = 'Severe', #legend_title
       x = expression("log"[2]*" fold change"), y = expression("-log"[10]*"p-value")) +

  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('DE genes in TCGA samples with TPM counts (MSI-H vs MSS & MSI-L)') + # Plot title
  geom_text_repel(max.overlaps = Inf) # To show all labels 

# Open the file that will contain your plot (the name is up to you)
png(file = "output/top_500_volcano_plot.png", width = 500, height = 800) # you can change the size of the output file
# Execute the plot
myvolcanoplot
# Close the file that will contain the plot
dev.off()
