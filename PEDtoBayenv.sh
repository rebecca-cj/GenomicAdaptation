## Convert PED format to Bayenv2 input file using PGDSpider ##

# Usage:  PEDtoBayenv.sh <plink-prefix> <path/.spid_file> <output_prefix>

input=${1##*/}

# adjust .spid file to indicate correct .map file 
	sed -i 's/PED_PARSER_MAP_FILE_QUESTION=.*$/PED_PARSER_MAP_FILE_QUESTION='$input'.map/' $2

# Convert PED file to Bayenv2 format using PGDSpider
	java -Xmx2g -jar /volume/softwares/pgdspider/PGDSpider_2.0.9.0/PGDSpider2-cli.jar -inputfile $1.ped -inputformat PED -outputfile $3.bayenv.txt -outputformat BAYENV -spid $2
