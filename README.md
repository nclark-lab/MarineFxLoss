# MarineFxLoss
This repo contains scripts associated with a genome-wide screen for convergent loss of gene function in marine mammals.

These scripts and associated data files (downloadable [here](https://pitt.box.com/s/mb9z57sxu3o2cgeru5s7u0bek2nbw1lr)) were used to perform analyses reported in the following manuscript:
```
Meyer WK, Jamison J, Richter R, Woods SE, Partha R, Kronk C, Chikina M, Bonde RK, Crocker DE, Gaspard J, Lanyon JM, Marsillach J, Furlong CE, and Clark NL. Ancient convergent losses of Paraoxonase 1 yield potential risks for modern marine mammals. In revision, Science.

```

## Contents

The sub-directories contain scripts for the following sets of analyses:
1) PseudogeneIdentification: scripts to parse 100-way alignment files and generate tables with pseudogene/functional/not assigned categorization for each species and each gene, and to estimate appropriate missingness filters
2) BayesTraits: script to run BayesTraits for association between gene functional loss and terrestrial-marine transition, with associated control files
3) Simulations: scripts to generate matched simulated genes for real genes based on gene trees and estimated functional loss rates, and to analyze the resulting simulated genes to estimate empirical p-values and FDR
4) FunctionalEnrichment: scripts to perform hypergeometric test of functional enrichment for the top N genes from a ranked gene list against functional databases 

## Authors

* **Maria Chikina** - [mchikina](https://github.com/mchikina)
* **Nathan Clark** - [nclark-lab](https://github.com/nclark-lab)
* **Wynn Meyer** - [sorrywm](https://github.com/sorrywm)
* **Raghavendran Partha** - [raghavendranpartha](https://github.com/raghavendranpartha)

## Acknowledgments

* Analysis of association between gene functional loss and terrestrial-marine transition was performed using the BayesTraits software (v3):
```
M. Pagel, A. Meade, Bayesian Analysis of Correlated Evolution of Discrete Characters by Reversible‐Jump Markov Chain Monte Carlo. Am. Nat. 167, 808–825 (2006).

```

* Scripts for functional enrichment and estimation of FDR were developed for analyses in the following paper:
```
Chikina M, Robinson JD, Clark NL. Hundreds of Genes Experienced Convergent Shifts in Selective Pressure in Marine Mammals. Mol Biol Evol. 2016;33: 2182–92. doi:10.1093/molbev/msw112

```

* Methods for estimating FDR rely on ideas from the following paper:
```
V. G. Tusher, R. Tibshirani, G. Chu, Significance analysis of microarrays applied to the ionizing radiation response. Proc. Natl. Acad. Sci. 98, 5116–5121 (2001).

```

* Data from the following resources provide the gene sets analyzed for functional enrichment:
```
A. Subramanian et al., Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A. 102, 15545–50 (2005).

J. A. Blake et al., Mouse Genome Database (MGD)-2017: community knowledge resource for the laboratory mouse. Nucleic Acids Res. 45, D723–D729 (2017).

```

* Gene alignments are derived from the following source:
```
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/multiz100way/

For information about all genomes included in this alignment, see the following:
http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_alignment
```
