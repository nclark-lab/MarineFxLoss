#Combine two pseudogene lists ####
stoptable <- 'pseudo_stops.tsv'
shifttable <- 'recall_frameshift_miss.tsv' 
maxmisstabnodups <- "MaxMissingExcludeExon25_nodups.tsv" #created from AnalyzeMissingData
lesiontabnodups <- "LesionStatusExcludeExon25_nodups.tsv" #created from AnalyzeMissingData
missthresh <- 0.16 #proportion missing threshold above which to exclude a gene ("x")
minabsence <- 2 #minimum number of species with gene absent
maxabsence <- 29 #maximum number of species with gene absent
mincallable <- 29 #minimum number of species with less than the maximum allowable missing data
mincallablestrict <- 39
outsuff1 <- 'miss16' #suffix for output filename
outsuff2 <- 'miss16_abs2-29_call29' #suffix for output filename
outsuff3 <- 'miss16_abs2-29_call39' #suffix for output filename with more stringent missingness filter
#pstop <- read.table(stoptable,header=T,as.is=TRUE)
#pshift <- read.table(shifttable,header=T,as.is=T)

#Much of this is now in the 'AnalyzeMissingDataForPseudogenes.R' script
source("subfunctions_missing_pseudo.R")

#Generate lesion table, if it doesn't already exist ####
if (!file.exists(lesiontabnodups)) {
  pstop <- read.table(stoptable,header=T,as.is=TRUE)
  pshift <- read.table(shifttable,header=T,as.is=TRUE)
  checkstopshift <- check_data_frames(pstop,pshift,"ucid")
  stopifnot(checkstopshift)
  pboth <- pstop
  pboth[,3:ncol(pboth)] <- pstop[,3:ncol(pstop)]+pshift[,3:ncol(pshift)]
  pboth[,3:ncol(pboth)] <- pboth[,3:ncol(pboth)] > 0
} else {
  pboth <- read.table(lesiontabnodups,header=T,as.is=T)
}

#Use missing data threshold to exclude gene x species pairs with missing data ####
#pmiss <- read.table("proportion_gaps_out_nospacenames.tsv",header=T,as.is=T)
pmissmax <- read.table(maxmisstabnodups,header=T,as.is=T)
checkmisslesion <- check_data_frames(pmissmax,pboth,"ucid")
stopifnot(checkmisslesion)

#replace with x if greater than missing data threshold
pwiththresh <- pboth
pwiththresh[,3:ncol(pwiththresh)] <- 1-pwiththresh[,3:ncol(pwiththresh)]
for (c in 3:ncol(pwiththresh)) {
  pwiththresh[which(pmissmax[,c]>=missthresh),c] <- "x"
}
write.table(pwiththresh,paste("pseudogenes_",outsuff1,".tsv",sep=""),row.names=F,quote=F)

#Determine how many genes can be included in the analysis using inclusion criteria:
#1) at least 2 species with "0"
whichincl1 <- which(rowSums(data.matrix(pwiththresh)==0,na.rm=T)>=minabsence)
#2) at most 29 species with "0"
whichincl2 <- which(rowSums(data.matrix(pwiththresh)==0,na.rm=T)<=maxabsence)
#3) at least 29 species with "0" or "1"
whichincl3 <- which((58-rowSums(pwiththresh=="x"))>=mincallable)
whichincl3strict <- which((58-rowSums(pwiththresh=="x"))>=mincallablestrict)
whichint <- intersect(intersect(whichincl1,whichincl2),whichincl3)
incltab <- pwiththresh[whichint,]
write.table(incltab,paste("pseudogenes_",outsuff2,"_filtered.tsv",sep=""),row.names=F,quote=F)
incltabnoucid <- incltab[,-which(names(incltab)=="ucid")]
write.table(incltabnoucid,paste("pseudogenes_",outsuff2,"_filtered_noucid.tsv",sep=""),row.names=F,quote=F)
#Use a more stringent callable filter
whichintstrict <- intersect(intersect(whichincl1,whichincl2),whichincl3strict)
incltabstrict <- pwiththresh[whichintstrict,]
write.table(incltabstrict,paste("pseudogenes_",outsuff3,"_filtered.tsv",sep=""),row.names=F,quote=F)
incltabstrictnoucid <- incltabstrict[,-which(names(incltabstrict)=="ucid")]
write.table(incltabstrictnoucid,paste("pseudogenes_",outsuff3,"_filtered_noucid.tsv",sep=""),row.names=F,quote=F)
            
#Assess distribution of number of species with missing data ####
numMissing = rowSums(pwiththresh=="x",na.rm=T)
hist(numMissing,breaks=58,xlab = "# species with missing data",
main = paste("Histogram of number of species with >",missthresh,"missing"))
quantile(numMissing,c(0.75,0.9,0.95,0.99))

#Compare with previous table of filtered pseudogenes ####
st1 <- read.csv("ST1_NoDates.csv",header=T,as.is=T)
names(st1)[which(names(st1)=="Gene")] <- "gene"
statcol <- which(names(st1) %in% c("OldGene_WithDates","N.species.","LL_Independent",          
                                  "GeneLossRate_Independent","blank",                   
                                  "LL_Dependent","MarineLossRate",          
                                  "TerrestrialLossRate", "LRT",                    
                                  "P.value","Description"))
st1comp <- st1[,-statcol]
convtable <- read.csv("SpeciesNameConversion.csv",header=T,as.is=T)
for (c in 2:ncol(st1comp)) {
  if (names(st1comp)[c] %in% convtable$nameinst1) {
    names(st1comp)[c] <- convtable$abbr[which(convtable$nameinst1==names(st1comp)[c])]
  } else {
    print(paste(names(st1comp)[c],"not found in conversion table"))
  }
}
if (ncol(st1comp) != ncol(pboth)) {
  if (ncol(st1comp) > ncol(pboth)) {
    missp <- names(st1comp)[which(names(st1comp) %in% names(pboth) == F)]
    print("Species in original Supp Table 1 missing from new pseudogene calls:")
    print(missp)
    st1comp <- st1comp[,-which(names(st1comp) %in% missp)]
  } else {
    missp <- names(pboth)[which(names(pboth) %in% names(st1comp) == F)]
    print("Species in new pseudogene calls missing from Supp Table 1:")
    print(missp)
    pboth <- pboth[,-which(names(pboth) %in% missp)]
  }
}

#Determine overlap between sets
samegenes <- intersect(st1comp$gene,incltab$gene)
print(paste("Of the",length(unique(st1comp$gene)),"original and",length(unique(incltab$gene)),
            "new genes for analysis,",length(unique(samegenes)),"genes are overlapping \n(",
            length(unique(samegenes))/length(unique(st1comp$gene)),",",
            length(unique(samegenes))/length(unique(incltab$gene)),")."))
#st1comp <- st1comp[st1comp$gene %in% samegenes,]
#st1comp <- st1comp[order(st1comp$gene),]
#pboth <- pboth[pboth$gene %in% samegenes,]
#pboth <- unique(pboth)
#pboth <- pboth[order(pboth$gene),]
#if (sum(st1comp$gene==pboth$gene) != nrow(st1comp)) {
#  print("Issue with keeping same genes in same order")
#}
#st1comp <- st1comp[,order(names(st1comp))]
#pboth <- pboth[,order(names(pboth))]
#if (sum(names(st1comp)==names(pboth)) != ncol(st1comp)) {
#  print("Issue with keeping species in same order")
#}
