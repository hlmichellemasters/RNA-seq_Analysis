# RNA-seq_Analysis
library(gplots)
library(RColorBrewer)
library(ape)

#########################################################
### B) Read in your Rockhopper data and transform it into matrix format
#########################################################

data <- read.csv("NC_999999_transcripts.txt", sep="\t", stringsAsFactors=FALSE)	# This reads the data, which is in a tab-delimited format. Change the filename as necessary.
count = 0

rnames <- data[,6]                            # assign labels in column 6 to "rnames" - these are the gene IDs

for(i in data[,6]) 
  count = count + 1
  rnames[[count]] <- data[[count,6]]
  if(i == "-") {
    rnames[[count]] <- data[[count,8]]			# If no gene ID is available, this replaces it with the gene product description - you can change this if you would prefer
  }
  if(rnames[[count]] == "-") {
    rnames[[count]] <- "Unknown Function"		# This is a hack to deal with all the unannotated genes in our genomes
  }
  
  ratio <- data[,9]/data[,10]					# This is just to establish the size/structure of "ratio" to be filled in later
  
  qvalue <- data[,11]
  
  for(i in rnames) {
    if(data[[count,9]] == 0) {
      ratio[[count]] <- log(data[[count,10]],2)			# If expression.1 condition is zero counts, make the ratio equal to the negative
    }
    if(data[[count,10]] == 0) {
      ratio[[count]] <- log(data[[count,9]],2)
    }
    if(data[[count,10]] != 0) {
      if(data[[count,9]] != 0) {
        ratio[[count]] <- log(data[[count,9]]/data[[count,10]], 2)		# These calculate the count-ratio, while being careful not to divide by zero
      }
    }
  }
  
  mat_data <- data.matrix(data[,9:ncol(data)])  # transform column 9 and on into a matrix
  rownames(mat_data) <- rnames                  # assign row names 
  
  total_genes <- length(rnames)		# Useful if you want to calculate proportion of differentially expressed genes
  
  sig_values <- which(qvalue<0.01)	# Include only the genes that have significant qvalues
  
  # sig_values <- which(qvalue<2)		# Use this as an alternative if you would rather include all samples, not just the significant ones
  
  total_diff <- length(sig_values)
  
  sig_mat_data <- mat_data[sig_values,]	# Here we are isolating just the values that correspond with the significant q values
  
  sig_ratios <- ratio[sig_values]
  
  sig_rnames <- rnames[sig_values]
  
  sig_mat_plot <- sig_mat_data[,1:2]
  
  rownames(sig_mat_plot) <- sig_rnames
  
  known_sig_values <- which(sig_rnames!="-")
  
  known_sig_mat_data <- sig_mat_data[known_sig_values,]
  
  known_sig_rnames <- sig_rnames[known_sig_values]
  
  known_sig_ratios <- sig_ratios[known_sig_values]
  
  known_sig_mat_plot <- known_sig_mat_data[,1:2]
  
  #########################################################
  ### C) Customizing and plotting the heat map
  #########################################################
  
  # creates a own color palette from red to green
  # my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  
  # (optional) defines the color breaks manually for a "skewed" color transition
  # col_breaks = c(seq(-1,0,length=100),  # for red
  #  seq(0,0.8,length=100),              # for yellow
  #  seq(0.8,1,length=100))              # for green
  
  # creates a 5 x 5 inch image
  pdf("heatmaps_in_r.pdf"    # create PNG for the heat map        
  )     
  
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  lhei = c(1.5,4,1)
  lwid = c(1.5,4)
  
  dend <- heatmap.2(known_sig_mat_plot, 
                    #  cellnote = known_sig_mat_plot,  # Adds the text for the values in each heat map cell
                    #  notecex = 0.5,		# Sets the font size for the cell entries
                    lmat = lmat,
                    lwid = lwid,
                    lhei = lhei,
                    notecol="black",      # change font color of cell labels to black
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    margins =c(11,20),     # widens margins around plot
                    col=colorRampPalette(c("blue","purple","pink"))(400),       # Defines color palette - there are many options for this, and it really comes down to preference. cm.colors, topo.colors, terrain.colors, rainbow, are a few other options. 
                    cexCol=.7,				# Sets font for sample labels
                    cexRow=.55,				# Sets font for gene labels
                    #  breaks=col_breaks,    # enable color transition at specified limits
                    #  dendrogram="row",     # only draw a row or column dendrogram - will show how multiple samples (if column) or genes (if row) cluster in their expression patterns
                    Colv="NA",		 # turn off column dendrogram
                    Rowv="NA"
  )           
  
  title("Expression Change", line= 1, cex.main=1.5)
  
  dev.off()               # close the PNG device
  
  # dend <- heatmap.2(known_sig_mat_plot)
  
  # row.dend <- dend$rowDendrogram
  # row.hclust <- as.hclust(row.dend)
  # row.tree <- as.phylo(row.hclust)
  # write.tree(row.tree, file="Row Cluster Tree.txt")
  
  pdf("Known Significant Genes.pdf")	# Creates a pdf file with a barplot of all the expression ratios for the significant genes.
  
  par(mar=c(4, max(nchar(known_sig_rnames)/1), 2, 2))	# This is a hack to get the margins big enough to fit the labels
  
  # par(mar=c(4, 5, 2, 2)
  
  barplot(known_sig_ratios, main="Known significantly differentially expressed genes", names=known_sig_rnames, cex.names=0.8, xlab="Expression Ratio log2(1/2)", las=2, col="blue", horiz=TRUE)
  
  dev.off()
  
  # pdf("Total_Diff_Genes.pdf")
  
  # barplot(total_diff, col="green", main="Total number of differentially expressed genes")
  
  # dev.off()
  
  sig_ratio_matrix <- rbind(sig_rnames, sig_ratios)
  
  write("Total Differences", file="total_different.txt")
  write(total_diff, file="total_different.txt", append=TRUE)	# This outputs a text file with the total number of differentially expressed genes, and the proportion of all genes tested that were significant
  write(sig_ratio_matrix, file="total_different.txt", append=TRUE, sep="\t", ncolumns=2)
  
  known_sig_ratio_matrix <- rbind(known_sig_rnames, known_sig_ratios)
  write("Known Differences", file="Known_differences.txt")
  write(known_sig_ratio_matrix, file="Known_differences.txt", append=TRUE, sep="\t", ncolumns=2)
  
