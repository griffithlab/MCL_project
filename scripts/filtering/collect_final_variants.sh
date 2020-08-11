#This script will collect the finalized variants.

#Loop through applicable directories.
for directory in /gscmnt/gc2547/griffithlab/bvli/MCL/filtering_analysis_gnomad_0.01/TWGE*
do
	#STEP 1: Change into directory.
	cd "$directory"

	#STEP 2: Grab sample name.
	currentsample=$(pwd | sed 's/\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\///')	

	#STEP 3: Use gawk on "$currentsample"_manually_reviewed_variants.mgibed.
	gawk -v sample="$currentsample" 'NR == 1 {for(i=0;i<=NF;i++) {if($i == "POS") {pos=i;} else if ($i == "SYMBOL") {symbol=i;} else if ($i == "Consequence") {consequence=i;} else if ($i == "#CHROM") {chromosome=i;} else {continue} } } NR > 1 { print sample "\t" $chromosome "\t" $pos "\t" $consequence "\t" $symbol "\t" $chromosome":"$pos }' <(cat "$currentsample"_manually_reviewed_variants.mgibed | sed -re '/^##/d') > "$currentsample"_manually_reviewed_variants_final.mgibed

	#STEP 4: Grab lines from pasted_"$currentsample".vep to avoid any de-duplication efforts.
	while read line
	do
		#egrep out lines from vep using current line.
		gawk -v start="$line" '{if($2 ~ start) {print}}' pasted_"$currentsample".vep >> "$currentsample"_manually_reviewed_variants_final.vep
	done < <(cat "$currentsample"_manually_reviewed_variants_final.mgibed | cut -f 6 | uniq | rev | sed 's/:/|/' | rev)
done