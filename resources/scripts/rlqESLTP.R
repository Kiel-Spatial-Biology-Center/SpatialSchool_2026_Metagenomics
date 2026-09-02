# Extended RLQ (rlqESLTP)
#
# Reproduced verbatim from the supplementary code of:
#   Pavoine S., Rodriguez-Girones M.A. & Cabana J. (2011) Extended RLQ:
#   linking environment, space, traits and phylogeny. Journal of Ecology.
#   (JEC_1743_sm_apps5.R)
#
# Requires: ade4
#
# The five ordinations it combines:
#   dudiE  environmental variables of the sites   (PCA)
#   dudiS  spatial structure of the sites         (PCA of neighbour scores)
#   dudiL  the community table itself             (correspondence analysis)
#   dudiT  functional traits of the species       (PCoA)
#   dudiP  phylogenetic distances between species (PCoA)

rlqESLTP <- function(dudiE, dudiS, dudiL, dudiT, dudiP, ...){
	tabE <- dudiE$li/sqrt(dudiE$eig[1])
	tabS <- dudiS$li/sqrt(dudiS$eig[1])
	tabP <- dudiP$li/sqrt(dudiP$eig[1])
	tabT <- dudiT$li/sqrt(dudiT$eig[1])
	tabES <- cbind.data.frame(tabE, tabS)
	names(tabES) <- c(paste("E", 1:ncol(tabE), sep = ""),
		paste("S", 1:ncol(tabS), sep = ""))
	tabTP <- cbind.data.frame(tabT, tabP)
	names(tabTP) <- c(paste("T", 1:ncol(tabT), sep = ""),
		paste("P", 1:ncol(tabP), sep = ""))
	pcaES <- dudi.pca(tabES, scale = F, row.w = dudiL$lw, scan = FALSE,
		nf = (length(dudiE$eig) + length(dudiS$eig)))
	pcaTP <- dudi.pca(tabTP, scale = F, row.w = dudiL$cw, scan = FALSE,
		nf = (length(dudiT$eig) + length(dudiP$eig)))

	X <- rlq(pcaES, dudiL, pcaTP, ...)

	U <- as.matrix(X$l1) * unlist(X$lw)
	U <- data.frame(as.matrix(pcaES$tab[, 1:ncol(tabE)]) %*% U[1:ncol(tabE), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
	names(U) <- names(X$lR)
	X$lR_givenE <- U
	

	U <- as.matrix(X$l1) * unlist(X$lw)
	U <- data.frame(as.matrix(pcaES$tab[, -(1:ncol(tabE))]) %*% U[-(1:ncol(tabE)), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
	names(U) <- names(X$lR)
	X$lR_givenS <- U


	U <- as.matrix(X$c1) * unlist(X$cw)
	U <- data.frame(as.matrix(pcaTP$tab[, 1:ncol(tabT)]) %*% U[1:ncol(tabT), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
	names(U) <- names(X$lQ)
	X$lQ_givenT <- U

	U <- as.matrix(X$c1) * unlist(X$cw)
	U <- data.frame(as.matrix(pcaTP$tab[, -(1:ncol(tabT))]) %*% U[-(1:ncol(tabT)), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
	names(U) <- names(X$lQ)
	X$lQ_givenP <- U

	X$row.w <- dudiL$lw

	X$col.w <- dudiL$cw
	
	class(X) <- c("rlqESLTP", "rlq", "dudi")

	return(X)
}

