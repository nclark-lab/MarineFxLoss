#!/bin/bash
# FULL LIST
#panTro4 gorGor3 ponAbe2 nomLeu3 rheMac3 macFas5 papHam1 chlSab1 calJac3 saiBol1 otoGar3 tupChi1 speTri2 jacJac1 micOch1 criGri1 mesAur1 mm10 rn5 hetGla2 cavPor3 chiLan1 octDeg1 oryCun2 ochPri3 susScr3 vicPac2 camFer1 turTru2 orcOrc1 panHod1 bosTau7 oviAri3 capHir1 equCab2 cerSim1 felCat5 canFam3 musFur1 ailMel1 odoRosDiv1 lepWed1 pteAle1 pteVam1 myoDav1 myoLuc2 eptFus1 eriEur2 sorAra2 conCri1 loxAfr3 eleEdw1 triMan1 chrAsi1 echTel2 oryAfe1 dasNov3


# download some genomes in 'Chains' file from UCSC
./download_genomes.pl


# download more genomes from UCSC in .fa.gz format
for i in panTro4 ponAbe2 nomLeu3 rheMac3 macFas5 tupChi1 micOch1 criGri1 rn5 ochPri3 vicPac2 bosTau7 musFur1 lepWed1 eriEur2
do
	wget http://hgdownload-test.cse.ucsc.edu/goldenPath/$i/bigZips/$i.fa.gz
done


# download more genomes from UCSC in .2bit format
for i in ponAbe2 mm10 equCab2 odoRosDiv1 conCri1
do
	wget http://hgdownload-test.cse.ucsc.edu/goldenPath/$i/bigZips/$i.2bit
done

# Change 2bit to fa format and zip
for i in ponAbe2 mm10 equCab2 odoRosDiv1 conCri1
do
	./twoBitToFa $i.2bit $i.fa
	gzip $i.fa
done


# download remaining genomes from NCBI
 # chlSab1 jacJac1 mesAur1 chiLan1 octDeg1 camFer1 panHod1 capHir1 pteAle1 myoDav1 eptFus1 eleEdw1 chrAsi1 oryAfe1
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/260/255/GCF_000260255.1_OctDeg1.0/GCF_000260255.1_OctDeg1.0_genomic.fna.gz -O octDeg1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/311/805/GCF_000311805.1_CB1/GCF_000311805.1_CB1_genomic.fna.gz -O camFer1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/400/835/GCF_000400835.1_PHO1.0/GCF_000400835.1_PHO1.0_genomic.fna.gz -O panHod1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.1_ARS1/GCF_001704415.1_ARS1_genomic.fna.gz -O capHir1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/325/575/GCF_000325575.1_ASM32557v1/GCF_000325575.1_ASM32557v1_genomic.fna.gz -O pteAle1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/327/345/GCF_000327345.1_ASM32734v1/GCF_000327345.1_ASM32734v1_genomic.fna.gz -O myoDav1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/308/155/GCF_000308155.1_EptFus1.0/GCF_000308155.1_EptFus1.0_genomic.fna.gz -O eptFus1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/299/155/GCF_000299155.1_EleEdw1.0/GCF_000299155.1_EleEdw1.0_genomic.fna.gz -O eleEdw1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/296/735/GCF_000296735.1_ChrAsi1.0/GCF_000296735.1_ChrAsi1.0_genomic.fna.gz -O chrAsi1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/298/275/GCF_000298275.1_OryAfe1.0/GCF_000298275.1_OryAfe1.0_genomic.fna.gz -O oryAfe1.fa.gz

exit
