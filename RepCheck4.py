#### Compare genotypes of reps - to gain an estimate of genotype error (for given alignment and SNP calling pipeline) ##
#### Allows for half score between homozygote and heterozygote matches 
#
#  Copyright 2016 Rebecca Jordan
#
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
#
#####

# NOTE:  raw counts for 'missing', 'fullmatch', 'halfmatch' and 'diff' are not correct for samples with > 2 replicates

##--INPUT--##
# input_genotypes = vcftools .012 output file with sample names added as 2nd column (1st = number; 2nd = sample name; 3rd++ = genotypes in 012 format)
# input_samplenames = text file of sample names (excluding replicate prefix or suffix); one sample name per line
# output_name = prefix for output file

# Usage: python RepCheck.py <input_genotypes> <input_samplenames> <output_name>


import sys
import fnmatch

genotypes=sys.argv[1]
samplereps=sys.argv[2]
outputname=sys.argv[3]

sampleList = open(samplereps).readlines()
genodata = open(genotypes).readlines()
numloci = len(genodata[0].strip().split()[2:])

errorDict={}
for lineX in sampleList:
    sampleX=lineX.strip()
    datamat={}
    for lineG in genodata:
        currentline=lineG.strip().split()
        if sampleX in currentline[1]:
            datamat[currentline[1]]=currentline[2:]
    missing=0
    halfmatch=0
    fullmatch=0
    diff=0
    locusSim=[]
    for locus in range(0,numloci):
        tmplocus=[]
        for reps in datamat.keys():
            tmplocus.append(datamat[reps][locus])
	scorex=[]
	for numx in range(0,(len(tmplocus)-1)):
	    for numy in range((numx+1),len(tmplocus)):
        	if int(tmplocus[numx]) == -1 or int(tmplocus[numy]) == -1:
		    scorex.append("NA")
		    missing = missing +1
		elif tmplocus[numx] == tmplocus[numy]:
                    scorex.append(1)
		    fullmatch = fullmatch +1
		elif abs(int(tmplocus[numx])-int(tmplocus[numy])) ==1:
                    scorex.append(0.5)
		    halfmatch = halfmatch +1
		else:
                    scorex.append(0)
		    diff = diff +1
	nonmiss = [num for num in scorex if num != "NA"]    
   	if len(nonmiss) > 0:
	    locusSim.append(float(sum(nonmiss))/float(len(nonmiss))*100)
    errorDict[sampleX]=str(numloci) + "," + str(missing) + "," + str(fullmatch) + "," + str(halfmatch) + "," + str(diff) + "," + str(float(sum(locusSim))/float(len(locusSim)))

outfile=open(outputname + '_error.csv','w')
outfile.write("NumLoci,Missing,FullMatch,HalfMatch,Diff,%Similarity\n")    
for samples in errorDict.keys():
	outfile.write(samples + "," + errorDict[samples] + "\n")	
outfile.close()

