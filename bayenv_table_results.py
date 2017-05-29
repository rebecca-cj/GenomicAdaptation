###################################################################
#                                                                 #
#   Convert Bayenv2 results to table of loci vs variables         #
#   with associated strengths of association (where applicable)   #
#                                                                 #
#   Copyright Rebecca Jordan 2016                                 #
#                                                                 #
###################################################################

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


# Input file:  bayenv results summaried into file with three columns - 1st = locus, 2nd = associated variable, 3rd = strength of association
# Usage bayenv_table_results.py <input_file> <number_of_variables>

import sys

INPUTfile=sys.argv[1]
OUTprefix=INPUTfile[:-9]
NumberVar=int(sys.argv[2])

fileIN=open(INPUTfile).readlines()

# create dictionary of loci, with all variable slots = 0
Loci={}
for lineX in range(1,len(fileIN)):
    raw=fileIN[lineX].strip().split()
    Loci[raw[0]]={}
    for varX in range(1,(NumberVar+1)):
        Loci[raw[0]][varX]=0

# populate dictionary with significant environmental associations (where applicable)
for lineX in range(1,len(fileIN)):
    raw=fileIN[lineX].strip().split()
    Loci[raw[0]][int(raw[1])]=raw[2]

outfile=open(OUTprefix + '.results.table','w')
outfile.write("Locus\t")
for varX in range(1,(NumberVar+1)):
    outfile.write("Var" + str(varX) + "\t")
outfile.write("\n")
for locusX in Loci.keys():
    outfile.write(str(locusX) + "\t")
    for PC in Loci[locusX].keys():
        outfile.write(str(Loci[locusX][PC]) + "\t")
    outfile.write("\n")
outfile.close()
