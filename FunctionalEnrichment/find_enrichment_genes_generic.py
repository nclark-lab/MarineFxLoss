#!/usr/bin/python
import sys, glob
from scipy.special import gammaln 
from math import exp,log
import math

import operator
import random


def get_pi0(pv, lambdas):
    """
    Compute Storey's C{pi0} from p-values C{pv} and C{lambda}.

    this function is equivalent to::
    
        m = len(pv)
        return [sum(p >= l for p in pv)/((1.0-l) * m) for l in lambdas]
        
    but the above is C{O(m*n)}, while it needs only be C{O(m+n)
    (n = len(lambdas))}

    @type pv: list
    @param pv: B{SORTED} p-values vector
    @type lambdas: list
    @param lambdas: B{SORTED} lambda values vector

    @rtype: list
    @return: estimated proportion of null hypotheses C{pi0} for each lambda
    """
    m = len(pv)
    i = m - 1
    pi0 = []
    for l in reversed(lambdas):
        while i >= 0 and pv[i] >= l:
            i -= 1
        pi0.append((m-i-1)/((1.0-l)*m))
    pi0.reverse()
    return pi0

def storey_qvalues(pv, l=None, pi0=None):
    """
    Return Storey FDR q-values corresponding to p-values C{pv}.

    The main difference between B-H's and Storey's q-values is that the latter
    are weighted by the estimated proportion C{pi0} of true null hypotheses.

    @type pv: list
    @param pv: p-values from a multiple statistical test
    @type l: float
    @param l: lambda value for C{pi0} (fraction of null p-values) estimation

    @rtype: list
    @return: storey q-values corresponding to C{pv}
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))

    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    if pi0 is None:
      if l is None:
          # optimal lambda/pi0 estimation
          lambdas = [i/100.0 for i in xrange(0, 95, 5)]
          n = len(lambdas)
          pi0 = get_pi0(pv, lambdas)
          min_pi0 = min(pi0)
          mse = [0] * n        
          for i in xrange(1, 101):
              # compute bootstrap sample with replacement
              pv_boot = [pv[int(random.random()*m)] for j in xrange(m)]
              pi0_boot = get_pi0(sorted(pv_boot), lambdas)
              for j in xrange(n):
                  mse[j] += (pi0_boot[j] - min_pi0) * (pi0_boot[j] - min_pi0)
          min_mse = min(mse)
          argmin_mse = [i for i, mse_i in enumerate(mse) if mse_i == min_mse]
          pi0 = min(pi0[i] for i in argmin_mse)
          
          pi0 = min(pi0, 1)
      

      else:
          try:
              l = float(l)
          except ValueError:
              raise TypeError("lambda must be a number")
          if l < 0 or l >= 1:
              raise ValueError("lambda must be within [0,1)")
          pi0 = get_pi0(pv, [l])
          pi0 = min(pi0[0], 1)
    sys.stderr.write("pi0="+str(pi0)+"\n")
    qvalues = m * [0]
    mincoeff = pi0 * pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = pi0*m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues



def lnchoose(n, m):
  nf = gammaln(n + 1)
  mf = gammaln(m + 1)
  nmmnf = gammaln(n - m + 1)
  return nf - (mf + nmmnf)

def hypergeometric_gamma(k, n1, n2, t):
  if k > n1 or k > n2:
    return 0
  else:
    res=0
    for i in range(k, min(n1, n2)+1):
        c1 = lnchoose(n1,i)
        c2 = lnchoose(t-n1, n2 - i)
        c3 = lnchoose(t ,n2)
        res+=exp(c1 + c2 - c3)
    return res


def bincoeff1(n, r):
  if r < n - r:
    r = n - r
  x = 1
  for i in range(n, r, -1):
    x *= i
  for i in range(n - r, 1, -1):
    x /= i
  return x

def hypergeometric(k, n1, n2, t):
  if t > n1 + n2:
    t = n1 + n2
  if k > n1 or k > t:
    return 0
  elif t > n2 and ((k + n2) < t):
    return 0
  else:
    res=0
    for i in range(k, min(n1, n2)):
        c1 = log(raw_bincoeff1(n1,k))
        c2 = log(raw_bincoeff1(n2, t - k))
        c3 = log(raw_bincoeff1(n1 + n2 ,t))
        res+=exp(c1 + c2 - c3)
    return res

class Enrichment:
     def __init__(self, go, found, back, numgenes, pval, adj_pval, genes=[]):
         self.go = go
         self.found=found
         self.back=back
         self.pval=pval
         self.adj_pval=adj_pval
         self.genes=genes
         self.numgenes=numgenes
     def __cmp__(self, other):
         return cmp(self.pval, other.pval)
     def __str__(self):
         return "%s\t%d\t%d\t%4.4e\t%4.4e" % (self.go,self.found,self.back,float(self.pval), float(self.adj_pval))+"\t"+", ".join(self.genes)
     def allStr(self, numGenes, numBack):
         return "%s\t%d\t%4.4f\t%d\t%4.4f\t%4.4f\t%4.4e\t%4.4e" % (self.go,self.found,self.found/float(numGenes), self.back, self.back/float(numBack),self.found/float(numGenes)/(self.back/float(numBack)), float(self.pval), float(self.adj_pval))+"\t"+", ".join(self.genes)
     @staticmethod
     def headerString():
         return "name\tnumber in foreground\tfraction in foreground\tnumber in background\tfraction in background\tfold enrichment\tp-value\tcorrected significance"
     



from optparse import OptionParser
usage="usage: options <file1> <file2>"
parser = OptionParser(usage)
parser.add_option("-a", "--anns", dest="anns",
                      help="annotation file")
parser.add_option("-t", "--trans", dest="trans",
                      help="translation file")
parser.add_option("-b", "--back", dest="back",
                      help="background file")
parser.add_option("-p", "--pval", type="float", dest="pval",
                      help="P-value cutoff", default=0.05)
parser.add_option("-r", "--remove", type="int", dest="remove",
                      help="Remove categories bellow count", default=3)
parser.add_option("-f", "--no-fdr",  action="store_false", dest="fdr", default=True)
parser.add_option("-q", "--no-qval",  action="store_false", dest="qval", default=True)
parser.add_option("-o", "--output", dest="output", help="output string")
#parser.add_option("-o", "--output",  action="store_true", dest="output", default=False, help="Create an output file with extension \'.genes_en\'")
parser.add_option("-A", "--all",  action="store_true", dest="all", default=False, help="Use the whole background list regardless of annotation status")


(options, args) = parser.parse_args()
dSeq2Common=dict()
if (options.pval):
  pval=options.pval
else:
  pval=0.05

if (options.trans):
  fhTrans=open(options.trans, 'rb')
  for line in fhTrans:
    lAnn=line[:-1].split('\t')
    if (lAnn[0] and lAnn[1]):
      dSeq2Common[lAnn[0]]=lAnn[1]
      

if (len(args)<1):
    print args[1]
lAnn=list()
#sGo=set()
fhAnn=open(options.anns, 'rb')
dAnnGene=dict()
for line in fhAnn:
    lAnn=line[:-1].split('\t')
    gene=lAnn.pop(0)
    if (not dAnnGene.has_key(gene)):
        sAnn=set()
        dAnnGene[gene]=sAnn
    else:
        sAnn=dAnnGene[gene]
    for a in lAnn:
        if not a == "":
            sAnn.add(a)
   #         sGo.add(a)

totalGenes=len(dAnnGene);

#compute background
dBack=dict();
sInBack=set();
dTooSmall=dict();
dFound=dict()
if (options.back):
  sys.stderr.write('reading in background\n')
  fhBack=open(options.back, 'rb')
  for line in fhBack:
    lTok=line[:-1].split('\t')
    gene=lTok[0]
    if (dAnnGene.has_key(gene) or options.all):
      sInBack.add(gene)
      if (dAnnGene.has_key(gene)):
        sAnn=dAnnGene[gene]
        for a in sAnn:
          dFound[a]=0;

else:
  for gene in dAnnGene.keys():
    sInBack.add(gene)

totalGenes=len(sInBack);
    
    
for k,tmpS in dAnnGene.items(): 
  if k in sInBack or not options.back:
    for go in tmpS:
      if not dBack.has_key(go):
        dBack[go]=1;
      else:
        dBack[go]+=1
for ann in dBack.keys():
  if (dBack[ann]<options.remove):
    dTooSmall[ann]=1
    #    sGo.remove(ann)

sGeneGO=set()
dGenes=dict()
aEnrich=[]
presentGenes=dict()
numGenes=0;

for genesF in args:
  presentGenes.clear()
  for go,found in dFound.items():
    dFound[go]=0
  dGenes.clear()
  aEnrich=[]
  numGenes=0;
  genesFH=open(genesF, 'rb')

  for line in genesFH:
    lTok=line[:-1].split('\t')
    gene=lTok[0]
    if (dSeq2Common.has_key(gene)):
      gene=dSeq2Common[gene]
    if (not presentGenes.has_key(gene)):
      presentGenes[gene]=1      
      if (dAnnGene.has_key(gene) and gene in sInBack):
        numGenes+=1;
        for go in dAnnGene[gene]:
          if (not dTooSmall.has_key(go)):
            if not dFound.has_key(go) or dFound[go]==0:
              dFound[go]=1;
              genesA=[]
              dGenes[go]=genesA;
              dGenes[go].append(gene)
            else:
              dFound[go]+=1
              dGenes[go].append(gene)
      elif (options.all and gene in sInBack):
        numGenes+=1
  for go,found in dFound.items():
    if (found>1):
      aEnrich.append(Enrichment(go, found, dBack[go],numGenes, hypergeometric_gamma(found, dBack[go], numGenes, totalGenes),-1,dGenes[go]))

  lMulti=len(aEnrich)

  sys.stderr.write( "processed genes = "+str( numGenes)+" ;background= "+str(totalGenes)+"\n")
  genesFH.close()    
  aEnrich.sort()
  sys.stderr.write( "tested "+str(lMulti)+" hypotheses\n")

  pvals=list()

  for i in range(0, len(aEnrich)):
    pvals.append( aEnrich[i].pval)
  if (options.fdr):
    if (options.qval):
      pvals=storey_qvalues(pvals)
    else:
      pvals=storey_qvalues(pvals, None,1)
  
  for i in range(0, len(aEnrich)):
    aEnrich[i].adj_pval=pvals[i]


  if (options.output):
    fhOut=open(genesF+".genesEN"+options.output, 'wb')
  else:
      fhOut=sys.stdout
    
  for i in range(0, len(aEnrich)):
    if (aEnrich[i].adj_pval>pval):
      del aEnrich[i:len(aEnrich)]
      break           
    if (i==0):
        fhOut.write(Enrichment.headerString()+"\n")
    fhOut.write(aEnrich[i].allStr(numGenes, totalGenes)+"\n")


