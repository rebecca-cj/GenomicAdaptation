###################################################
##                                               ##
##     Calculate genotype error by SNP           ##
##  by comparing replicate sample genotypes      ##
##                                               ##
##   (average % similarity between replicates    ##
##       allowing for half match between         ##
##        homozygotes and heterozygotes)  	     ##
##                                               ##
##  Copyright Sept 2016 R Jordan                 ##
##                                               ##
###################################################

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

##--- INPUT ---###
# input_genotypes = vcftools .012 output file with sample names added as 2nd column (1st = number; 2nd = sample name; 3rd++ = genotypes in 012 format)
# input_samplenames = text file of sample names (excluding replicate prefix or suffix); one sample name per line
# output_name = prefix for output file

# Usage: python RepCheck_similariy_bylocus.py <input_genotypes> <input_samplenames> <output_name>

import sys

genotypes=sys.argv[1]
samplereps=sys.argv[2]
outputname=sys.argv[3]

sampleList = open(samplereps).readlines()
genodata = open(genotypes).readlines()
numloci = len(genodata[0].strip().split()[2:])

lociDict={}
for yesterday in range(0,numloci):
    lociDict[yesterday]=[]

for lineX in sampleList:
    sampleX=lineX.strip()
    datamat={}
    for lineG in genodata:
        currentline=lineG.strip().split()
        if sampleX in currentline[1]:
            datamat[currentline[1]]=currentline[2:]
    for locus in range(0,numloci):
        tmplocus=[]
        for reps in datamat.keys():
            tmplocus.append(datamat[reps][locus])
            scorex = []
        for numx in range(0,(len(tmplocus)-1)):
            for numy in range((numx+1),len(tmplocus)):
                if int(tmplocus[numx]) == -1 or int(tmplocus[numy]) ==-1:
                    scorex.append("NA")
                elif tmplocus[numx] == tmplocus[numy]:
                    scorex.append(1)
                elif abs(int(tmplocus[numx])-int(tmplocus[numy])) ==1:
                    scorex.append(0.5)
                else:
                    scorex.append(0)
        nonmiss = [num for num in scorex if num != "NA"]
        if len(nonmiss) > 0:        
            lociDict[locus].append(float(sum(nonmiss))/float(len(nonmiss))*100)
        else:
            lociDict[locus].append("NA")

outfile=open(outputname + '_locierror.csv','w')
outfile.write("LocusNum,Missing,AvgSame\n")
for loci in lociDict.keys():
    miss = lociDict[loci].count("NA")
    nonmiss = [num for num in lociDict[loci] if num != "NA"]
    outfile.write(str(loci) + "," + str(miss) + ",")
    if len(nonmiss) > 0:
        outfile.write(str(float(sum(nonmiss))/float(len(nonmiss))) + "\n")
    else:
        outfile.write(str("na") + "\n")
outfile.close()
