#!/bin/bash

for directoryname in /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1031 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-009 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-390 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-263 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-233 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-910 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-519 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1070 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-411 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-040 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-881 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-187 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1122 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-362 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1441 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-526 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1497 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-686 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-355 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1174 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1369 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-496 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1027 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-886 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-724 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-115 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-602 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-930 
do
	#Change directory.
	cd $directoryname
	
	#Grab the normal file.
        normalcalcthresholdsfile=$(find . -name "calculated*normal*")

	#Grab the tumor file.
        tumorcalcthresholdsfile=$(find . -name "calculated*tumor*")

	#Grab current sample name.
	samplename=$(pwd | sed 's/\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\///')
	
	#Remove .\ from normalcalcthresholdsfile.
	newnormalcalcthresholdsfile=$(echo $normalcalcthresholdsfile | sed 's/\.\///')

	#Remove .\ from tumorcalcthresholdsfile.
        newtumorcalcthresholdsfile=$(echo $tumorcalcthresholdsfile | sed 's/\.\///')
	
	#Run awk to average coverage depths for normal.
	gawk -v sample="$samplename" '{x += $7; y += $8; z += $9; a += $10; b += $11; c += $12; d += $13; e += $14;} END {print sample"_normal" "\t" x/NR "\t" y/NR "\t" z/NR "\t"  a/NR "\t"  b/NR "\t" c/NR "\t" d/NR "\t" e/NR}' $newnormalcalcthresholdsfile > averaged_"$newnormalcalcthresholdsfile"

	#Run awk to average coverage depths for tumor.
        gawk -v sample="$samplename" '{x += $7; y += $8; z += $9; a += $10; b += $11; c += $12; d += $13; e += $14;} END {print sample"_tumor" "\t" x/NR "\t" y/NR "\t" z/NR "\t"  a/NR "\t"  b/NR "\t" c/NR "\t" d/NR "\t" e/NR}' $newtumorcalcthresholdsfile > averaged_"$newtumorcalcthresholdsfile"

	#Add this to "master" normal file.
	cat averaged_"$newnormalcalcthresholdsfile" >> /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/new_normal_master_averaged.bed

	#Add this to "master" tumor file.
        cat averaged_"$newtumorcalcthresholdsfile" >> /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/new_tumor_master_averaged.bed
done

#Change directory.
cd /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/

#Average across the normal_master_averaged.bed file.
gawk '{x += $2; y += $3; z += $4; a += $5; b += $6; c += $7; d += $8; e += $9;} END {print $1 "\t" x/NR "\t" y/NR "\t" z/NR "\t"  a/NR "\t"  b/NR "\t" c/NR "\t" d/NR "\t" e/NR}' new_normal_master_averaged.bed  > final_new_normal_master_averaged.bed

#Average across the tumor_master_averaged.bed file.
gawk '{x += $2; y += $3; z += $4; a += $5; b += $6; c += $7; d += $8; e += $9;} END {print $1 "\t" x/NR "\t" y/NR "\t" z/NR "\t"  a/NR "\t"  b/NR "\t" c/NR "\t" d/NR "\t" e/NR}' new_tumor_master_averaged.bed  > final_new_tumor_master_averaged.bed
