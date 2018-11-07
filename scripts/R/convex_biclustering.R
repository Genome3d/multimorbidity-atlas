#!/usr/bin/env Rscript

#Convex Biclustering of phenotypes according to shared eGenes or eQTLs
suppressPackageStartupMessages(
  for (package in c('gplots', 'cvxbiclustr', 'stringi', 'RColorBrewer')) {
    if (!require(package, character.only=T, quietly=T, warn.conflicts=F)) {
      install.packages(package, repos="http://cran.us.r-project.org")
      library(package, character.only=T)
    }
  }
)


#Make matrix
cat('Reading eGene matrix...\n')
pb <- txtProgressBar(min=0, max=100, style=3)
mdata <- read.delim("../../results/egene_matrix.txt", header = T)
mdata = subset(mdata, mdata$Genes1 > 3 & mdata$Genes2 > 3 & mdata$Common_eGenes > 3)
mdata$Trait1 <- stri_replace_all_fixed(mdata$Trait1,
                                       c("_", "  "),
                                       c(" ", ", "), 
                                       vectorize_all = FALSE)
mdata$Trait2 <- stri_replace_all_fixed(mdata$Trait2,
                                       c("_", "  "),
                                       c(" ", ", "), 
                                       vectorize_all = FALSE)
xtraits <- sort(unique(mdata$Trait1))
ytraits <- sort(unique(mdata$Trait2))
setTxtProgressBar(pb,100)
close(pb)
trait.matrix <- matrix(, nrow = length(ytraits), ncol = length(ytraits), dimnames = list(ytraits, ytraits))

cat('Creating matrix...\n')
pb <- txtProgressBar(min=0, max=100, style=3)
z = 0
for (i in rownames(trait.matrix)){
  z <- z + 1
  y <- (z * 100) / length(ytraits)
  setTxtProgressBar(pb,y)
  for (j in colnames(trait.matrix)){
    i.index <- which(mdata$Trait1 == i)
    for (x in i.index){
      if (mdata[x,2] == j){
        trait.matrix[i,j] <- mdata[x,6]
      }
      
    }
  }
}
close(pb)
trait.matrix <- replace(trait.matrix, is.na(trait.matrix), 0)


if (!require("RColorBrewer")) {
install.packages("RColorBrewer")
library(RColorBrewer)
}


#Annotation for heatmap
types <- colnames(trait.matrix)
ty <- as.numeric(factor(types))
cols <- colorRampPalette(brewer.pal(12, "Paired"))(625)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(100)

# Construct weights and edge-incidence matrices
cat('Constructing weights...\n')
wts <- gkn_weights(trait.matrix) # combines Gaussian kernel weights with k-nearest neighbor weights
w_row <- wts$w_row  # Vector of weights for row graph
w_col <- wts$w_col  # Vector of weights for column graph
E_row <- wts$E_row  # Edge-incidence matrix for row graph
E_col <- wts$E_col  # Edge-incidence matrix for column graph

# Initialize solution path gamma parameters 
gammaSeq <- 10^seq(0, 3, length.out = 100)


#Perform validation to select regularization parameter
cat('Performing validation...\n')
solution <- cobra_validate(trait.matrix,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)


# Plot validation error
verr <- solution$validation_error
plot(verr)


## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- solution$groups_row[[ix]]
groups_col <- solution$groups_col[[ix]]

M <- biclust_smooth(trait.matrix,groups_row,groups_col)

cat('Ploting graph...')
pdf("../../results/egenes_bicluster.pdf", width=100, height=100)
heatmap.2(M,col=hmcols,labRow = rownames(trait.matrix), labCol = colnames(trait.matrix),ColSideCol=cols[ty],margins=c(40,30),
          main = "Smooth Convex Biclustering of Phenotypes according to shared eGenes", trace="none",
          keysize = 0.2, density.info = "none", key.xlab = "Proportion of shared eGenes",# key.par = list(mar=c(3,3,3,3)),
          lmat=rbind(c(0,4,0), c(0, 1,0), c(3,2,0), c(0,0,5)), lwid=c(0.1,1,0.1), lhei=c(0.1,0.01,1,0.05))
dev.off()
cat('\nDone.')

#End of eGene Convex Biclustering

#End of eQTL Convex Biclustering
mdata <- read.delim("../../results/eqtl_matrix.txt", header = T, row.names=NULL)
mdata$Trait2 <- mdata$Trait1
mdata$Trait1 <- mdata$row.names
mdata = subset(mdata, mdata$eQTLs1 > 3 & mdata$eQTLs2 > 3)


mdata$Trait1 <- stri_replace_all_fixed(mdata$Trait1,
                                       c("_", "  "),
                                       c(" ", ", "), 
                                       vectorize_all = FALSE)
mdata$Trait2 <- stri_replace_all_fixed(mdata$Trait2,
                                       c("_", "  "),
                                       c(" ", ", "), 
                                       vectorize_all = FALSE)

xtraits <- sort(unique(mdata$Trait1))
ytraits <- sort(unique(mdata$Trait2))
setTxtProgressBar(pb,100)
close(pb)
trait.matrix <- matrix(, nrow = length(ytraits), ncol = length(ytraits), dimnames = list(ytraits, ytraits))

cat('Creating matrix...\n')
pb <- txtProgressBar(min=0, max=100, style=3)
z = 0
for (i in rownames(trait.matrix)){
  z <- z + 1
  y <- (z * 100) / length(ytraits)
  setTxtProgressBar(pb,y)
  for (j in colnames(trait.matrix)){
    i.index <- which(mdata$Trait1 == i)
    for (x in i.index){
      if (mdata[x,3] == j){
        trait.matrix[i,j] <- mdata[x,6]
      }
      
    }
  }
}
close(pb)
#trait.matrix[!is.numeric(trait.matrix)] <- NA
trait.matrix <- replace(trait.matrix, is.na(trait.matrix), 0)




#Annotation for heatmap
types <- colnames(trait.matrix)
ty <- as.numeric(factor(types))
cols <- colorRampPalette(brewer.pal(12, "Paired"))(625)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(100)

# Construct weights and edge-incidence matrices
cat('Constructing weights...\n')
wts <- gkn_weights(trait.matrix) # combines Gaussian kernel weights with k-nearest neighbor weights
w_row <- wts$w_row  # Vector of weights for row graph
w_col <- wts$w_col  # Vector of weights for column graph
E_row <- wts$E_row  # Edge-incidence matrix for row graph
E_col <- wts$E_col  # Edge-incidence matrix for column graph

# Initialize solution path gamma parameters 
gammaSeq <- 10^seq(0, 3, length.out = 100)


#Perform validation to select regularization parameter
cat('Performing validation...\n')
solution <- cobra_validate(trait.matrix,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)


# Plot validation error
verr <- solution$validation_error
plot(verr)


## Heatmap of data smoothed at the model selected to minimize validation error
ix <- which.min(verr)
groups_row <- solution$groups_row[[ix]]
groups_col <- solution$groups_col[[ix]]

M <- biclust_smooth(trait.matrix,groups_row,groups_col)

cat('Ploting graph...')
pdf("../../results/eqtls_bicluster.pdf", width=100, height=100)
heatmap.2(M,col=hmcols,labRow = rownames(trait.matrix), labCol = colnames(trait.matrix),ColSideCol=cols[ty],margins=c(40,30),
          main = "Smooth Convex Biclustering of Phenotypes according to shared eQTLs", trace="none",
          keysize = 0.2, density.info = "none", key.xlab = "Proportion of shared eQTLs",# key.par = list(mar=c(3,3,3,3)),
          lmat=rbind(c(0,4,0), c(0, 1,0), c(3,2,0), c(0,0,5)), lwid=c(0.1,1,0.1), lhei=c(0.1,0.01,1,0.05))
dev.off()
cat('\nDone.')

#End of eQTL Convex Biclustering

