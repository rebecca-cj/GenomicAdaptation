######################################################
#                                                    #
#  Filter Bayenv2 environmental association results  #
#  for SNPs significantly associated with variables  #
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
# Assumes fewer than 10000 SNPs being tested
# Associations filtered using Kass and Raftery interpretation of Bayes Factors (as per Bayenv2 manual)
# 'Positive' = 3 < BF < 20; 'Strong' = 20 <= BF < 150; 'VeryStrong = BF > 150
# SNP numbers adjusted by 1 to match numbered loci in .numbered.map file, used to convert SNP numbers to SNP loci (Bayenv counts from 0, .numbered.map counts from 1)
# Input: Bayenv2 environmental association OUTFILE (by default = "bf_environ.ENVIRONFILE.txt"); SNPs as rows, Environmental variables as columns
# Input: indicate (Y/N) whether Spearman's p was also calculated during analysis (Y = BF & Spearman's p; N = just BF)
# Usage: bayenv_enviro_summary.py <input_file> <output_prefix> <number_variables> <rho_calculated_Y/N>

import sys
import os

dataIN=sys.argv[1]
prefixOUT=sys.argv[2]
countVar=sys.argv[3]
rhocalc=sys.argv[4]

print "number of variables considered:  " + countVar

if rhocalc=="Y":
    calccolumns=range(1,((int(countVar)*3)+1),3)
elif rhocalc=="N":
    calccolumns=range(1,(int(countVar)+1))

print "data columns used to assess BF (first column = 0):  " + str(calccolumns)

fileIN=open(dataIN).readlines()
EnviroVar={}
for count in range(0,int(countVar)):
    datacolumn=calccolumns[count]
    var=int(count)+1
    EnviroVar[var]={}
    tmp3=[]
    tmp20=[]
    tmp150=[]
    for lineX in fileIN:
        raw=lineX.strip().split()
        if float(raw[datacolumn])>3 and float(raw[datacolumn])<20:
            tmp3.append(int(raw[0][-4:])+1)
        elif float(raw[datacolumn])>=20 and float(raw[datacolumn])<150:
            tmp20.append(int(raw[0][-4:])+1)
        elif float(raw[datacolumn])>=150:
            tmp150.append(int(raw[0][-4:])+1)
    EnviroVar[var]['positive']=tmp3
    EnviroVar[var]['strong']=tmp20
    EnviroVar[var]['verystrong']=tmp150

# create results directory
os.mkdir("./"+prefixOUT+"_snps")
    
# Output count of number of SNPs associated with each variable at each "association strength"
outsummary=open("./"+prefixOUT+"_snps/"+prefixOUT+'_BFsummary.txt','a')
outsummary.write("Variable\tPositive\tStrong\tVeryStrong\n")
for var in EnviroVar.keys():
    outsummary.write(str(var) + "\t" + str(len(EnviroVar[var]['positive'])) + "\t" + str(len(EnviroVar[var]['strong'])) + "\t" + str(len(EnviroVar[var]['verystrong'])) + "\n")
outsummary.close()

# Output SNPs (numbers) associated with each variable at each "association strength"
outSNPs=open("./"+prefixOUT+"_snps/"+prefixOUT+'_SNPassoc.txt','a')
outSNPs.write("Variable\tBFstrength\tSNPs\n")
for var in EnviroVar.keys():
    for BF in EnviroVar[var].keys():
        outSNPs.write(str(var) + "\t" + str(BF) + "\t")
        for i in range(0,len(EnviroVar[var][BF])):
            if int(EnviroVar[var][BF][i])<10:
                outSNPs.write(str(EnviroVar[var][BF][i])[-1:] + "\t")
            elif int(EnviroVar[var][BF][i])<100:
                outSNPs.write(str(EnviroVar[var][BF][i])[-2:] + "\t")
            elif int(EnviroVar[var][BF][i])<1000:
                outSNPs.write(str(EnviroVar[var][BF][i])[-3:] + "\t")
            else:
                outSNPs.write(str(EnviroVar[var][BF][i]) + "\t")
        outSNPs.write("\n")
outSNPs.close()

# Output individual files for each variable and SNP strength
for var in EnviroVar.keys():
    for BF in EnviroVar[var].keys():
        outfile=open("./"+prefixOUT+"_snps/"+str(var)+"."+BF+".SNPs",'a')
        for i in range(0,len(EnviroVar[var][BF])):
            if int(EnviroVar[var][BF][i])<10:
                outfile.write(str(EnviroVar[var][BF][i])[-1:] + "\n")
            elif int(EnviroVar[var][BF][i])<100:
                outfile.write(str(EnviroVar[var][BF][i])[-2:] + "\n")
            elif int(EnviroVar[var][BF][i])<1000:
                outfile.write(str(EnviroVar[var][BF][i])[-3:] + "\n")
            else:
                outfile.write(str(EnviroVar[var][BF][i]) + "\n")
        outfile.close()
