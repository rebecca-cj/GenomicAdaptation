##################################################################################
#                                                                                #
#  Create lists of neutral and 'outlier' loci found across all runs of BayeScan  #
#      	               	       	    	       	       	       	      	       	 #
#      (directional selection only; no balancing selection found in study        #
#                this script was originally written for)                         #
#      	               	       	    	       	       	       	      	       	 #
#     Copyright Oct 2016  Rebecca Jordan                                         #
#      	               	       	    	       	       	       	      	       	 #
##################################################################################

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
# Usage:  bayescan_outliers.sh <BayeScan_output_file_prefix> <output_prefix> <path/'numbered'.map_file>

# Create list of SNPs under directional selection
    # logBF>0.5
	:> 05.tmp
        for f in $1-*_fst.txt; do awk 'NR>1 {if($3>=0.5 && $5>0) print $1}' $f >> 05.tmp; done
	sort 05.tmp | uniq -c | awk '{if($1==3) print $2}' - > $2.substantial
	rm 05.tmp
    # logBF>2
	:> 2.tmp
        for f in $1-*_fst.txt; do awk 'NR>1 {if($3>=2 && $5>0) print $1}' $f >> 2.tmp; done
 	sort 2.tmp | uniq -c | awk '{if($1==3) print $2}' - > $2.decisive
	rm 2.tmp
    # qval<0.05
	:> q.tmp
        for f in $1-*_fst.txt; do awk NR>1 {if($4<=0.05 && $5>0) print $1}' $f >> q.tmp; done
	sort q.tmp | uniq -c | awk '{if($1==3) print $2}' - > $2.qval05
	rm q.tmp

# Create list of 'neutral' SNPs based on q>0.2
	:> n.tmp
        for f in $1-*_fst.txt; do awk '{if($4>0.2) print $1}' $f >> n.tmp; done
	sort n.tmp | uniq -c | awk '{if($1==3) print $2}' - > $2.neutral
	rm n.tmp

# Get list of actual outlier (and neutral) SNPs for further analysis; using 'numbered' .map file to convert SNP 'numbers' (from BayeScan output) to actual SNP loci
        for f in $2.substantial $2.decisive $2.qval05 $2.neutral; do awk 'FNR==NR{a[$1]=$0;next}; $1 in a {print a[$1]}' $3 $f | awk '{split($3,a,":"); print a[1] "\t" a[2]}' - > $f.loci; done
