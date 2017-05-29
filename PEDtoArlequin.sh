## Convert (adjusted as necessary) PED files to Arlequin input file using PGDSpider ##

# Note script must be run from same location as PED files (doesn't accept files in different directory)

# Usage:  PEDtoArlequin.sh <plink-prefix> <path/.spid_file> <output_prefix>

# adjust .spid file to indicate correct .map file
	sed -i 's/PED_PARSER_MAP_FILE_QUESTION=.*$/PED_PARSER_MAP_FILE_QUESTION='$1'.map/' $2

# convert from PED to Arlequin using PGDSpider
       java -Xmx2g -jar ~/data/java/PGDSpider_2.0.7.4/PGDSpider2-cli.jar -inputfile $1.ped -inputformat PED -outputfile $3.arp -outputformat ARLEQUIN -spid $2
