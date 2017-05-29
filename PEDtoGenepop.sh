## Convert (adjusted as necessary) PED files to GENEPOP format using PGDSpider ##

# Usage:  PEDtoFDist2.sh <plink-prefix> <path/.spid_file> <output_prefix>

# adjust .spid file to indicate correct .map file 
	sed -i 's/PED_PARSER_MAP_FILE_QUESTION=.*$/PED_PARSER_MAP_FILE_QUESTION='$1'.map/' $2

# convert from PED to FDist2 using PGDSpider
       java -Xmx2g -jar ~/data/java/PGDSpider_2.0.7.4/PGDSpider2-cli.jar -inputfile $1.ped -inputformat PED -outputfile $3.txt -outputformat GENEPOP -spid $2
