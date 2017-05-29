############################################
##                                        ##
##   Get directional selection outliers   ##
##    from Arlequin fdist2 results        ##
##                                        ##
##   Copyright Feb 2016  R Jordan         ##
##     	    	       	       	       	  ##
############################################

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


# Create lists of 'outlier' loci and associated significance values from Arlequin (directional outliers only; no loci under balancing selection in this study)

# Input: fdist2_ObsOut.txt Arlequin output file with additional 7th column of qvalues (adjusted FST pvalues) added (can have different name)
# Input: statistic and cutoff value to use to create the outlier list; p = p-value, q = q-value

# Usage: arl_outliers_plus_significance.sh <meanFst> <arlequin_output_file> <output_prefix> <statistic> <cutoff> <path/'numbered'.map_file>

if [ $4 == "p" ]
then
    echo "Using p-values with cutoff p<=$5"
    awk -v fst=$1 -v cutoff=$5 '{if($4<=cutoff && $3>fst) print $1 "\t" $4}' $2 > tmp.outfile
elif [ $4 == "q" ]
then
    echo "Using q-values with cutoff q<=$5"
    awk -v fst=$1 -v cutoff=$5 '{if($7<=cutoff && $3>fst) print $1 "\t" $7}' $2 > tmp.outfile
fi

# Get list of actual outlier SNPs for further analysis; using 'numbered' .map file to convert SNP 'numbers' (from arlequin output) to actual SNP $

    awk 'FNR==NR{a[$1]=$3;next} $1 in a {print a[$1]"\t" $2}' $6 tmp.outfile > $3.$4$5.significance

# Cleanup

    rm tmp.outfile
