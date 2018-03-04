#Input files ####
nucmisstab <- "NucMissingExcludeExon25.tsv" #missingness using nucleotide sequence
aamisstab <- "AAMissingExcludeExon25.tsv" #missingness using amino acid sequence
specieslist <- "SpeciesToInclude58Eutherian.txt" #list of species to include
stoptable <- 'LatestStopsExcludeExon25.tsv' #called stop codons using same exclusion criteria
shifttable <- 'LatestFrameshiftExcludeExon25.tsv' #called frameshifts using same exclusion criteria

#Output files ####
maxmisstab <- "MaxMissingExcludeExon25.tsv" #maximum missingness across nuc and aa for relevant species
excludedduplist <- "DuplicateGenesWithMultipleIsoformsSameMissingData.txt" #duplicate isoforms to exclude
allduptable <- "DuplicateGeneStatus.tsv" #table of duplicate genes, ucids, and inclusion/exclusion
maxmisstabnodups <- "MaxMissingExcludeExon25_nodups.tsv"
lesiontabnodups <- "LesionStatusExcludeExon25_nodups.tsv"

#Import subfunctions ####
source("subfunctions_missing_pseudo.R")

#Get maximum proportion of missing data from nuc and aa data ####
pmissnuc <- read.table(nucmisstab,header=T,as.is=T)
pmissaa <- read.table(aamisstab,header=T,as.is=T)
inclsp <- scan(specieslist,what="character")
colstoincl <- c(1,2,which(names(pmissnuc) %in% inclsp))
pmissnuc <- pmissnuc[,colstoincl]
pmissaa <- pmissaa[,colstoincl]
checknucaa <- check_data_frames(pmissnuc,pmissaa,"ucid")
stopifnot(checknucaa) 
pmissmax <- pmissnuc
for (sp in c(3:ncol(pmissmax))) {
  pmissmax[,sp] <- pmax(pmissaa[,sp],pmissnuc[,sp]) #get maximum for each cell
}
write.table(pmissmax,file=maxmisstab,row.names=FALSE,quote=FALSE)

#Determine which orthologs to filter out for duplicate genes ####
dupgenes <- unique(pmissmax$gene[duplicated(pmissmax$gene)])
duptab <- pmissmax[pmissmax$gene %in% dupgenes,]
duptabtowrite <- duptab[,1:2]
duptabtowrite$include <- c(rep(TRUE,nrow(duptabtowrite)))
duptab$totmiss <- rowSums(duptab[,3:ncol(duptab)])
rowstorem <- list(length(dupgenes))
nondiffgenes <- list(length(dupgenes))
for (w in 1:length(dupgenes)) {
  mindup <- duptab[which(duptab$gene==dupgenes[w]),]
  ucidrem <- mindup$ucid[which(mindup$totmiss!=min(mindup$totmiss))]
  if (length(ucidrem)==(nrow(mindup)-1)) {
    rowstorem[[w]] <- which(pmissmax$ucid %in% ucidrem)
    nondiffgenes[[w]] <- NA
  } else {
    print(paste("Multiple isoforms for",dupgenes[w],"have minimum missingness (",min(mindup$totmiss),")"))
    ucidrem <- as.vector(mindup$ucid)
    rowstorem[[w]] <- which(pmissmax$ucid %in% ucidrem)
    nondiffgenes[[w]] <- dupgenes[w]
  }
  duptabtowrite$include[which(duptabtowrite$ucid %in% ucidrem)] <- FALSE
}
nondiffgenes <- unlist(nondiffgenes)
nondiffgenes <- nondiffgenes[is.na(nondiffgenes)==F]
write(nondiffgenes,file=excludedduplist,ncolumns=1)
write.table(duptabtowrite,file=allduptable,row.names=FALSE,quote=FALSE)
duprowstorem <- unlist(rowstorem)
duprowstorem <- duprowstorem[is.na(duprowstorem)==F] #save these values and eliminate from all working tables

#Read in and combine stops and frameshifts ####
pstop <- read.table(stoptable,header=T,as.is=TRUE)
pshift <- read.table(shifttable,header=T,as.is=TRUE)
pstop <- pstop[,colstoincl]
pshift <- pshift[,colstoincl]
checkstopshift <- check_data_frames(pstop,pshift,"ucid")
stopifnot(checkstopshift)
pboth <- pstop
pboth[,3:ncol(pboth)] <- pstop[,3:ncol(pstop)]+pshift[,3:ncol(pshift)]
pboth[,3:ncol(pboth)] <- pboth[,3:ncol(pboth)] > 0

#Remove duplicate genes from all tables and write to files ####
pmissmax <- pmissmax[-duprowstorem,]
pboth <- pboth[-duprowstorem,]
write.table(pmissmax,file=maxmisstabnodups,row.names=FALSE,quote=FALSE)
write.table(pboth,file=lesiontabnodups,row.names=FALSE,quote=FALSE)

#Determine error rates for a given missingness threshold ####
checkmisslesion <- check_data_frames(pmissmax,pboth,"ucid")
stopifnot(checkmisslesion)
#Find the threshold where there will be no more than 1 error (in any species) for every 10 genes
accepterrorrate <- 1/(10*(ncol(pboth)-2))
props <- c(1:20)/20
err <- list(length(props))
bestval <- NA
for (p in 1:c(length(props))) {
  err[[p]] <- calcerrorrates(props[p],pmissmax[,3:ncol(pmissmax)],pboth[,3:ncol(pboth)])
  if (err[[p]][1] < accepterrorrate) {
    bestval <- p
  } else {
    break
  }
}
print(paste("The largest missingness cutoff that will result in no more that one",
            "mis-called functional gene across all species for every 10 genes in",
            "the dataset (functional false positive rate =",accepterrorrate,") is",
            props[bestval],"."))
print("The error rates associated with this missingness threshold are:")
print(err[[bestval]])
if (bestval < length(props)) {
  newprops <- props[bestval]+.01*c(1:4)
  newerr <- list(length(newprops))
  newbestval <- NA
  for (p in 1:c(length(newprops))) {
    newerr[[p]] <- calcerrorrates(newprops[p],pmissmax[,3:ncol(pmissmax)],pboth[,3:ncol(pboth)])
    if (newerr[[p]][1] < accepterrorrate) {
      newbestval <- p
    } else {
      break
    }
  }
  print(paste("The largest missingness cutoff that will result in no more that one",
              "mis-called functional gene across all species for every 10 genes in",
              "the dataset (functional false positive rate =",accepterrorrate,") is",
              newprops[newbestval],"."))
  print("The error rates associated with this missingness threshold are:")
  print(newerr[[newbestval]])
}

#Plot some estimates of the distribution of missing data among pseudogenes and non-pseudogenes ####
tlpmiss <- unlist(pmissmax)
tlpboth <- unlist(pboth)
pmisspseudo <- tlpmiss[which(tlpboth==TRUE)]
pmissnonpseudo <- tlpmiss[which(tlpboth==FALSE)]

#Plot histogram of missing data
hist(as.numeric(pmisspseudo),breaks=100,xlab="Proportion missing data",
     main="Histogram of proportion missing for callable pseudogenes\n(stops and frameshifts)")
qcuts <- quantile(as.numeric(pmisspseudo),c(.9,.95,.99))
abline(v=qcuts,col=c("navy","blue","lightblue"))
legend(0.8,50000,legend=c(paste("q90:",qcuts[1]),paste("q95",qcuts[2]),paste("q99",qcuts[3]))
       ,col=c("navy","blue","light blue"),lty=1,cex=0.7,bty="n")
#Viewed as cdf
epmiss <- ecdf(as.numeric(pmisspseudo))
enpmiss <- ecdf(as.numeric(pmissnonpseudo))
epmissnoext <- ecdf(as.numeric(pmisspseudo[pmisspseudo %in% c("0","1")==F]))
enpmissnoext <- ecdf(as.numeric(pmissnonpseudo[pmissnonpseudo %in% c("0","1")==F]))
plot(epmiss,xlab="Proportion missing",ylab="Cumulative fraction",main="CDF of proportion missing")
plot(enpmiss,add=T,col="blue")
plot(epmissnoext,add=T,col="gray")
plot(enpmissnoext,add=T,col="light blue")
legend("right",legend=c("Called pseudogenes\n(stops & frameshifts)","All other",
                        "Called pseudogenes, no 0","All other, no 0 or 1"),lty=1,
       col=c("black","blue","gray","light blue"))
