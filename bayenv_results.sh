#######################################################################
#                                                                     #
#   Process and summaries Bayenv2 environmental association results   #
#                                                                     #
#	Copyright 2016  Rebecca Jordan                                    #
#                                                                     #
#######################################################################

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


### INFORMATION ###

# Assumes fewer than 10000 SNPs being tested
# Associations filtered using Kass and Raftery interpretation of Bayes Factors (as per Bayenv2 manual)
# 'Positive' = 3 < BF < 20; 'Strong' = 20 <= BF < 150; 'VeryStrong = BF > 150
# SNP numbers adjusted by 1 to match numbered loci in .numbered.map file, used to convert SNP numbers to SNP loci (Bayenv counts from 0, .numbered.map counts from 1)

# REQUIRES the following scripts
#    bayenv_enviro_summary.py
#    bayenv_summary_to_loci.sh
#    bayenv_table_results.py

### USAGE ###

# Usage: bayenv_results.sh <script_dir> <input_file> <output_prefix> <num_var> <rho_calc> <path/numbered.map>

# Input 'script_dir': directory containing additional scripts required
# Input 'input_file': Bayenv2 environmental association OUTFILE (by default = "bf_environ.ENVIRONFILE.txt"); SNPs as rows, Environmental variables as columns
# Input 'output_prefix': prefix to append to all files created
# Input 'num_var': number of variables analysed
# Input 'rho_calc': indicate (Y/N) whether Spearman's p was also calculated during analysis (Y = BF & Spearman's p; N = just BF)
# Input 'path/numbered.map':  "Numbered" .map (plink) file (output from adjustPED.sh) where additional column created such that column 1 = SNP 'number'

### OUTPUT ###

# <prefix>_snps/: directory of results
# <prefix>_BFsummary.txt:  count of SNPs associated with each variable at each association 'strength'
# <prefix>_SNPassoc.txt:  list, by number, of SNPs associated with each variable, for each 'strength' (numbers adjusted to start from 1 rather than 0 as in Bayenv2 results)
# 'variable'.'strength'.SNPs:  set of files containing all SNPs, by number, associated with particular variable at particular association strength
# 'variable'.'strength'.loci:  set of files containing all SNPs, by actual locus, associated with particular variable at particular association strength
# <prefix>.combined:  single file of all associations, one row per SNP, 1st column = SNP, 2nd column = variable, 3rd = association strength.  SNPs can appear multiple times for associations with multiple variables
# <prefix>.results.table:  results as a table, SNPs by row, variables by column.  0 = no association (BF<3).  "Strength" as per Kass and Raftery (see above).

### SCRIPT ###

# summarise results including creation of files with list of SNPs per variable per strength
	echo "performing initial summary, incl. creating separate files for SNPs by variable"
	python $1/bayenv_enviro_summary.py $2 $3 $4 $5

# convert SNP numbers to SNP loci
	echo "converting SNP numbers to loci"
	cd $3_snps/
	for i in *.SNPs; do bash $1/bayenv_summary_to_loci.sh $i $6 ${i%.SNPs}; done

# combine results into single file
	echo "combining results into single file"
	printf "Locus\tVariable\tStrength\n" > $3.combined; for v in `seq 1 $4`; do for s in positive strong verystrong; do awk -v num=$v -v strength=$s '{print $1 ":" $2 "\t" num "\t" strength}' $v.$s.loci >> $3.combined; done; done

# convert single file, with multiple entries per SNP, to table format
	echo "converting combined file to table"
	python $1/bayenv_table_results.py $3.combined $4

	echo "completed"
