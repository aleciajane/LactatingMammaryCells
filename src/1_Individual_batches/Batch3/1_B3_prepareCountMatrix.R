## This script extracts a cleaned count matrix from the 10X output and is largely based on:
## https://github.com/MarioniLab/DropletUtils

#For each data file (3x NMC1B + 1x LMC9), I have repeated the below code

# First remove reads derived from index swapping
library(DropletUtils)
#####################NMC1B_1################################
samples <-  c("../../data/Batch3/RB_M39/1/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
		     barcode.length=bc.length,
		     get.swapped=get.swapped,
		     min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("RB1")

finalCounts <- list()

for (i in seq_along(cleanedCounts)) {
    m <- cleanedCounts[[i]]

    # test for empty droplets
    out <- emptyDrops(m)#,BPPARAM=bpparam)

    # set FDR threshold
    nonempty <- out$FDR < 0.01

    #remove NAs from NA FDR
    nonempty[is.na(nonempty)] <- FALSE

    # cant use logic vector to subset dgCMatrix?
    nonempty <- rownames(out[nonempty,])

    # Clean matrix
    m.clean <- m[,nonempty]
    finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../../data/Batch3/RB_M39/1/CountMatrix_RB1.rds")

#####################NMC1B_2################################
samples <-  c("../../data/Batch3/RB_M39/2/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
                              barcode.length=bc.length,
                              get.swapped=get.swapped,
                              min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("RB2")

finalCounts <- list()

for (i in seq_along(cleanedCounts)) {
  m <- cleanedCounts[[i]]
  
  # test for empty droplets
  out <- emptyDrops(m)#,BPPARAM=bpparam)
  
  # set FDR threshold
  nonempty <- out$FDR < 0.01
  
  #remove NAs from NA FDR
  nonempty[is.na(nonempty)] <- FALSE
  
  # cant use logic vector to subset dgCMatrix?
  nonempty <- rownames(out[nonempty,])
  
  # Clean matrix
  m.clean <- m[,nonempty]
  finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../../data/Batch3/RB_M39/2/CountMatrix_RB2.rds")

#####################NMC1B_3################################
samples <-  c("../../data/Batch3/RB_M39/3/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
                              barcode.length=bc.length,
                              get.swapped=get.swapped,
                              min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("RB3")

finalCounts <- list()

for (i in seq_along(cleanedCounts)) {
  m <- cleanedCounts[[i]]
  
  # test for empty droplets
  out <- emptyDrops(m)#,BPPARAM=bpparam)
  
  # set FDR threshold
  nonempty <- out$FDR < 0.01
  
  #remove NAs from NA FDR
  nonempty[is.na(nonempty)] <- FALSE
  
  # cant use logic vector to subset dgCMatrix?
  nonempty <- rownames(out[nonempty,])
  
  # Clean matrix
  m.clean <- m[,nonempty]
  finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../../data/Batch3/RB_M39/3/CountMatrix_RB3.rds")

#####################LMC9_17025################################
samples <-  c("../../data/Batch3/HMC_17025/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
                              barcode.length=bc.length,
                              get.swapped=get.swapped,
                              min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("HMC")

finalCounts <- list()

for (i in seq_along(cleanedCounts)) {
  m <- cleanedCounts[[i]]
  
  # test for empty droplets
  out <- emptyDrops(m)#,BPPARAM=bpparam)
  
  # set FDR threshold
  nonempty <- out$FDR < 0.01
  
  #remove NAs from NA FDR
  nonempty[is.na(nonempty)] <- FALSE
  
  # cant use logic vector to subset dgCMatrix?
  nonempty <- rownames(out[nonempty,])
  
  # Clean matrix
  m.clean <- m[,nonempty]
  finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../../data/Batch3/HMC_17025/CountMatrix_HMC.rds")