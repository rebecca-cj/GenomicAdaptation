######################################################
#                                                    #
#  Process and summarise population differentiation  #
#          (XtX) analysis from Bayenv2               #
#                                                    #
#  Copyright Rebecca Jordan 2016                     #
#                                                    #
######################################################

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


# -- INFORMATION -- #
# Input:  Output file from Population Differentiation analysis; 1st column = SNP number, 2nd column = XtX
# Output 1: Text file ('output_prefix'_TopXtX.txt) with top 1, 5 and 10% of SNPs based on ranked XtX
# Output 2: Histogram of XtX values.  Vertical line to indicte XtX cutoff for top 10% (solid), 5% (dashed), and 1% (dotted) of SNPs (based on ranked XtX values) 
# Usage: XtX_summary.py <input_file> <output_prefix> 

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

fileIN=sys.argv[1]
prefixOUT=sys.argv[2]

# import data & create dictionary of SNPs (numbers) and XtX
# n.b. Because Bayenv numbers SNPs from 0, SNP numbers increased by one to match SNP numbers in .numbered.map file 
#      .numbered.map file used to reconcile SNP numbers with actual loci with SNP numbers starting from 1
dataIN=open(fileIN).readlines()
SNPs={}
for lineX in dataIN:
    raw=lineX.strip().split()
    SNPs[int(raw[0][-4:])+1]=float(raw[1])

# calculate XtX cutoff value for top 1% (cut99), 5% (cut95) and 10% (cut90) of SNPs
inSORTED=sorted(SNPs.values())
cut90=inSORTED[int(len(inSORTED)*0.9)]
cut95=inSORTED[int(len(inSORTED)*0.95)]
cut99=inSORTED[int(len(inSORTED)*0.99)]

# create lists of top 1%, 5% and 10% of SNPs (SNP numbers)
topSNP={}
topSNP["top10percent"]=[]
topSNP["top5percent"]=[]
topSNP["top1percent"]=[]
for snps in SNPs.keys():
    if SNPs[snps] >= cut99:
        topSNP['top1percent'].append(snps)
    if SNPs[snps] >= cut95:
        topSNP['top5percent'].append(snps)
    if SNPs[snps] >= cut90:
        topSNP['top10percent'].append(snps)

# output as three separate files containing list of SNP numbers
for groups in topSNP.keys():
    outfile=open(prefixOUT+"_XtX."+groups[:-7],'w')    
    for i in range(0,len(topSNP[groups])):
        outfile.write(str(topSNP[groups][i]) + '\n')
    outfile.close()

# create histogram of XtX values, indicating XtX value cutoff of top 1% (dotted), 5% (dashed) and 10% (solid) of SNPs
plt.hist(SNPs.values(),bins=(range(int(min(SNPs.values()))-1,int(max(SNPs.values()))+2)),color='0.9')
plt.xlabel("XtX - Differentiation")
plt.ylabel("Number SNPs")
plt.axvline(x=int(cut90),color='r',linewidth=1.25)
plt.axvline(x=int(cut95),color='r',linestyle='dashed',linewidth=1.25)
plt.axvline(x=int(cut99),color='r',linestyle='dotted',linewidth=1.25)
plt.savefig(prefixOUT + 'XtXhist.png',bbox_inches='tight')
