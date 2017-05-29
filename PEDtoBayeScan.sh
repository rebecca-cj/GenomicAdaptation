## Convert PED format to BayeScan input file using PGDSpider ##

# Usage:  PEDtoBayeScan.sh <plink-prefix> <path/.spid_file> <output_prefix>

# adjust .spid file to indicate correct .map file 
	sed -i 's/PED_PARSER_MAP_FILE_QUESTION=.*$/PED_PARSER_MAP_FILE_QUESTION='$1'.map/' $2

# Convert PED file to BayeScan format using PGDSpider
	java -Xmx2g -jar ~/data/java/PGDSpider_2.0.7.4/PGDSpider2-cli.jar -inputfile $1.ped -inputformat PED -outputfile $3.BayeS.txt -outputformat GESTE_BAYE_SCAN -spid $2
