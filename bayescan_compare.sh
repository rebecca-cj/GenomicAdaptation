## Check consistency between BayeScan runs ##

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


# Input: BayeScan output files from three independant runs (same parameters) to be compared
#        Files should be named 'PREFIX-*_fst.txt' where * = 1, 2 and 3 for the three separate runs
# Usage:  bayescan_compare.sh <BayeScan_output_prefix> <output_prefix>  

# Number of 'outlier' loci (directional selection only)
        echo ">> Number of outliers <<" > $2.compare
        for f in $1-*_fst.txt; do printf "$f\t" >> $2.compare; awk '{if($3>=0.5 && $5>0) PO05++; if($3>=2 && $5>0) PO2++; if($4<=0.05 && $5>0) qval++}END{print "log(PO)>=0.5\t" PO05 "\tlog(PO)>=2.0\t" PO2 "\tqval<=0.05\t" qval}' $f >> $2.compare; done

# Compare outputs to check consistency between runs
    # logBF>0.5
        echo ">> Number shared loci logBF>0.5 <<" >> $2.compare
        for f in $1-*_fst.txt; do awk 'NR>1 {if($3>=0.5 && $5>0) print}' $f > ${f%_fst.txt}.05outloci; done
	:> $2_logBF05.loci
        for f in $1-*.05outloci; do cut -d ' ' -f 1 $f >> $2_logBF05.loci; done; sort $2_logBF05.loci | uniq -c | awk '{if($1==3) all++; if($1==2) two++; if($1==1) one++}END{print "all\t" all "\ttwo\t" two "\tone\t" one "\ttotal\t" NR}' -  >> $2.compare
    # logBF>2
        echo ">> Number shared loci logBF>2 <<" >> $2.compare
        for f in $1-*_fst.txt; do awk 'NR>1 {if($3>=2 && $5>0) print}' $f > ${f%_fst.txt}.2outloci; done
	:> $2_logBF2.loci
        for f in $1-*.2outloci; do cut -d ' ' -f 1 $f >> $2_logBF2.loci; done; sort $2_logBF2.loci | uniq -c | awk '{if($1==3) all++; if($1==2) two++; if($1==1) one++}END{print "all\t" all "\ttwo\t" two "\tone\t" one "\ttotal\t" NR}' - >> $2.compare
    # qval<0.05
        echo ">> Number shared loci qval<=0.05 <<" >> $2.compare
        for f in $1-*_fst.txt; do awk 'NR>1 {if($4<=0.05 && $5>0) print}' $f > ${f%_fst.txt}.q05outloci; done
	:> $2_qval05.loci
        for f in $1-*.q05outloci; do cut -d ' ' -f 1 $f >> $2_qval05.loci; done; sort $2_qval05.loci | uniq -c | awk '{if($1==3) all++; if($1==2) two++; if($1==1) one++}END{print "all\t" all "\ttwo\t" two "\tone\t" one "\ttotal\t" NR}' - >> $2.compare
