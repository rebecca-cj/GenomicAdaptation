## Create files of directional outliers and neutral SNPs from Arlequin 'fdist2' output (no loci under balancing selection output) ##

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


# Input: fdist2_ObsOut.txt Arlequin output file with additional 7th column of qvalues (adjusted FST pvalues) added (can have different name)
# Usage: arl_outliers.sh <meanFst> <arlequin_output_file> <output_prefix> <path/'numbered'.map_file>

# p=0.05
	awk -v fst=$1 '{if($4<=0.05 && $3>fst) print}' $2 > $3_p05.outliers
# p=0.01
	awk -v fst=$1 '{if($4<=0.01 && $3>fst) print}' $2 > $3_p01.outliers
# q=0.1
        awk '{if($7<=0.1 && $3>fst) print}' $2 > $3_qval1.outliers
# 'Neutral' SNPs (p > 0.1; theoretically not significant)
	awk '{if($4>0.1) print}' $2 > $3_neutral.outliers

# Get list of actual outlier (and neutral) SNPs for further analysis; using 'numbered' .map file to convert SNP 'numbers' (from Arlequin output) to actual SNP loci
	for f in $3_p05.outliers $3_p01.outliers $3_qval1.outliers $3_neutral.outliers; do awk 'FNR==NR{a[$1]=$0;next}; $1 in a {print a[$1]}' $4 $f | awk '{split($3,a,":"); print a[1] "\t" a[2]}' - > ${f%.*}.loci; done
