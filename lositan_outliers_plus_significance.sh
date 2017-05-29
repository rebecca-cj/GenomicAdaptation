## Create lists of 'outlier' loci and associated significance values from Lositan (directional selection only)

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


# Input: Prefix of FDR corrected Lostian output files; corrections made using 'correct.pval.dataframe' function, from getPvalue.R by Lotterhos & Whitlock 2014 (Mol. Ecol. 23:2178-92)
#        Input format = output format of 'correct.pval.dataframe'
#        Files should be named 'PREFIX*_loci_plus_q.txt' where * = 1, 2 and 3 for the three separate runs
# Input: cutoff value q-value to use to create the outlier list
# Usage:  lositan_outliers_plus_significance.sh <input_file_prefix> <output_prefix> <cut_off_qvalue>

echo "Using cutoff of qvalue<=$3"
    # create list of loci found in all three runs
        :> los.tmp
        for f in $1*_loci_plus_q.txt; do awk -v cutoff=$3 'NR>1 {if($7<=cutoff && $5=="R") print $1}' $f >> los.tmp; done
        sort los.tmp | uniq -c | awk '{if($1==3) print $2 "\t"}' - > tmp.outfile
        rm los.tmp
    # add significance from each run to locus
        for f in $1*_loci_plus_q.txt; do awk -v cutoff=$3 'NR>1 {if($7<=cutoff && $5=="R") print $1 "\t" $7}' $f | while read -r Locus BF;
        do sed -i "s/^$Locus\t/$Locus\t$BF\t/g" tmp.outfile; done; done

# Calculate average for statistic across three runs

    awk '{print $1, ($2+$3+$4)/3}' tmp.outfile > $2".q"$3".significance"

# Cleanup

    rm tmp.outfile

