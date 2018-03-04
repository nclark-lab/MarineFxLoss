#Generate a set of {ngenes} matched pseudogenes for a single gene's relative rate of loss
#Input: gene name
#Output: table of pseudogenes for input to BayesTraits
args <- commandArgs(TRUE)
ngenes = 100000 #request 01/06/18: Can we also have this be an input?
wgene <- args[1] #gene to match rates for in the simulations
pseudotable <- 'pseudogenes_bt_out.tsv' #table containing processed BayesTraits output for full dataset
require(dplyr)
require(data.table)
if (!'geiger' %in% installed.packages()) {
    install.packages("geiger",repos="http://cran.r-project.org")
    require(geiger)
} else {
    require(geiger) #added 05/08/17
}
psdUse <- fread(pseudotable, header = T)
ratetree <- read.nexus(paste0('./bt_trees/',wgene,'.nex')) #requires a bt_trees directory containing pruned trees for each gene
rlr <- filter(psdUse, gene == wgene)$genelossrate #loss rate for the gene to simulate
print(paste("Simulating",ngenes,"genes matched for",wgene,"loss rate (",rlr,")"))

Q = matrix(c(-1,0,1,0),2,2)

simonegene_nofilter <- function(rtree, tomult, tmat, ntosim){
        #generate matched data set based on a single gene's loss rate profile
	#rtree: rate tree for gene
	#tomult: rate by which to multiply branch lengths (gene's independent loss rate)
	#tmat: transition matrix (Q)
        #ntosim: # of simulations to run in sim.char
        grtree <- rtree
        grtree$edge.length <- rtree$edge.length*tomult
        stre <- sim.char(grtree,tmat,model='discrete',root=1,nsim=ntosim)
        sprof = aperm(stre,c(1,3,2))
        tswz <- t(sprof[,,1])
        tswz[tswz == 2] <- 0
	return(tswz)
}

simonegene <- function(rtree, tomult, tmat, ntosim, ntorep, minzeros, maxzeros){
        #generate matched data set based on a single gene's loss rate profile
        #ntosim: # of simulations to run in sim.char
        #ntorep: # of rows (P/A per species) to output
        #minzeros and maxzeros: limits on how many zeros to include in a row (2 and 29)
        grtree <- rtree
        grtree$edge.length <- rtree$edge.length*tomult
        stre <- sim.char(grtree,tmat,model='discrete',root=1,nsim=ntosim) #make nsim lower? just need 1
        sprof = aperm(stre,c(1,3,2))
        tswz <- t(sprof[,,1])
        tswz[tswz == 2] <- 0
        rst <- rowSums(tswz)
	minspec <- max(length(grtree$tip.label)-maxzeros,0)
	maxspec <- length(grtree$tip.label)-minzeros
        spui <- intersect(which(rst >= minspec),which(rst <= maxspec))
        if (length(spui) >= ntorep) {
                matchedsim <- tswz[sample(spui,ntorep),]
                return(matchedsim)
        } else {
                print("Not enough genes matched filtering criteria")
                return(c(rep(NA,ncol(tswz))))
        }
}

simonegene_old <- function(rtree, tomult, tmat, ntosim, ntorep, minspec, maxspec){
	#generate matched data set based on a single gene's loss rate profile
	#ntosim: # of simulations to run in sim.char
	#ntorep: # of rows (P/A per species) to output
	#minspec and maxspec: limits on how many non-missing species to include a row (31 and 56)
	grtree <- rtree
	grtree$edge.length <- rtree$edge.length*tomult
	stre <- sim.char(grtree,tmat,model='discrete',root=1,nsim=ntosim) #make nsim lower? just need 1
	sprof = aperm(stre,c(1,3,2))
	tswz <- t(sprof[,,1])
	tswz[tswz == 2] <- 0
	rst <- rowSums(tswz)
	spui <- intersect(which(rst >= minspec),which(rst <= maxspec))
	if (length(spui) >= ntorep) {
		matchedsim <- tswz[sample(spui,ntorep),]
		return(matchedsim)
	} else {
		print("Not enough genes matched filtering criteria")
		return(c(rep(NA,ncol(tswz))))
	}
}

#Run all simulations, with or without filtering
#tsimwithzero <- simonegene(ratetree,rlr,Q,ngenes*10,ngenes,2,29) #modified 01/06/18 to provide input to new function with min and max zeros, rather than species
tsimwithzero <- simonegene_nofilter(ratetree,rlr,Q,ngenes) #keep all simulations, no filtering
outmat = cbind(paste0('sim',1:ngenes),tsimwithzero)
#write.table(outmat, file = paste('outfiles/matchedsimulatedgenes_n',ngenes,'_g',wgene,'.tsv',sep=''), quote = F, row.names = F, col.names = T, sep = '\t')
write.table(outmat, file = paste('nfoutfiles/matchedsimulatedgenes_n',ngenes,'_g',wgene,'_nofilter.tsv',sep=''), quote = F, row.names = F, col.names = T, sep = '\t')

print(paste("Finished simulating",ngenes,"genes matched for",wgene,"loss rate (",rlr,")"))
