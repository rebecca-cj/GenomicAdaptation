## Create lists of 'outlier' loci and associated significance values from BayeScan (directional selection only; no balancing selection found in study this script was originally written for)

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


# Input: Prefix of BayeScan output files (use bayescan_compare.sh to compare multiple runs prior to running this script)
#        Files should be named 'PREFIX-*_fst.txt' where * = 1, 2 and 3 for the three separate runs
# Input: statistic and cutoff value to use to create the outlier list; BF = BayesFactor, Q = q-values
# Usage:  bayescan_outliers_plus_significance.sh <BayeScan_output_file_prefix> <output_prefix> <statistic> <cut_off_value> <path/'numbered'.map_file>

if [ $3 == "BF" ]
then
    echo "Using BF with cutoff $4"
    # create list of loci found in all three runs
        :> BF.tmp
        for f in $1-*_fst.txt; do awk -v cutoff=$4 'NR>1 {if($3>=cutoff && $5>0) print $1}' $f >> BF.tmp; done
        sort BF.tmp | uniq -c | awk '{if($1==3) print $2 "\t"}' - > tmp.outfile
        rm BF.tmp
    # add significance from each run to locus
        for f in $1-*_fst.txt; do awk -v cutoff=$4 'NR>1 {if($3>=cutoff && $5>0) print $1 "\t" $3}' $f | while read -r Locus BF;
        do sed -i "s/^$Locus\t/$Locus\t$BF\t/g" tmp.outfile; done; done
elif [ $3 == "Q" ]
then
    echo "Using q-value with cutoff $4"	
    # create list of loci found in all three runs
	:> q.tmp
        for f in $1-*_fst.txt; do awk -v cutoff=$4 'NR>1 {if($4<=cutoff && $5>0) print $1}' $f >> q.tmp; done
	sort q.tmp | uniq -c | awk '{if($1==3) print $2 "\t"}' - > tmp.outfile
	rm q.tmp
    # add significance from each run to	locus
        for f in $1-*_fst.txt; do awk -v cutoff=$4 'NR>1 {if($4<=cutoff && $5>0) print $1 "\t" $4}' $f | while read -r Locus Qval;       
        do sed -i "s/^$Locus\t/$Locus\t$Qval\t/g" tmp.outfile; done; done
fi

# Calculate average for statistic across three runs

    awk -v stat=$3 'BEGIN{print "Locus\tavg"stat}{print $1, ($2+$3+$4)/3}' tmp.outfile > tmp.outfile2

# Get list of actual outlier SNPs for further analysis; using 'numbered' .map file to convert SNP 'numbers' (from BayeScan output) to actual SNP loci

    awk 'FNR==NR{a[$1]=$3;next} $1 in a {print a[$1]"\t" $2}' $5 tmp.outfile2 > $2.$3$4.significance

# Clean up

#    rm tmp.outfile tmp.outfile2
