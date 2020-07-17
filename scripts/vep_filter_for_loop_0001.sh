#Loop through all TWGE* directories.
for directory in /gscmnt/gc2547/griffithlab/bvli/MCL/filtering_analysis_gnomad_0.01/TWGE*
do
	#Change directoiry.
	cd "$directory"	

	#Grab sample name.
	samplename=$(pwd | sed 's/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\///')

	#Grab annotated.filtered.vcf.gz.
	VCF=$(find . -type f -name 'annotated_filtered.vcf.gz' | sed 's/\.\///')	

	#Grab .vep file.
	VEP=$(ls | egrep '^TWGE-08-0075-[0-9][0-9][0-9].vep$|^TWGE-08-0075-[0-9][0-9][0-9][0-9].vep$')

	#Paste the two together.
	paste <(zcat "$VCF" | sed -re '/^#|^##/d') <(cat "$VEP" | sed -re '/^#|^##/d') > pasted_"$samplename".vep
done

#Loop through all TWGE* directories.
for directory in /gscmnt/gc2547/griffithlab/bvli/MCL/filtering_analysis_gnomad_0.01/TWGE* 
do
	#STEP 1: Change into directory.
	cd "$directory"

	#STEP 2: Grab the pasted_ file.
	pastedfile=$(find . -type f -name "pasted_*" | sed 's/\.\///')	

	#STEP 3: Grab sample name.
	samplename=$(pwd | sed 's/\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\/[^/]*\///')

	#STEP 4: Filter pasted_* using gawk.
	gawk '{m=split($9,a,":"); n=split($10,b,":"); o=split($11,c,":"); for(i=0;i<=m;i++) {if(a[i] == "DP") {normaldepth=b[i]; tumordepth=c[i];} if(a[i] == "AD") {split(b[i],d,","); split(c[i],e,","); normalvaf=(d[2])/(d[1]+d[2]); tumorvarsupportreads=e[2];}} if((normaldepth >= 20 && tumordepth >= 20) && (normalvaf < 0.05) && (tumorvarsupportreads >= 5)) {print}}' "$pastedfile" > "$samplename"_basic_filtered.vep

	#STEP 5: Filter vepfile using gawk.
	gawk '!/incomplete_terminal_codon_variant|start_retained_variant|stop_retained_variant|synonymous_variant|mature_miRNA_variant|5_prime_UTR_variant|3_prime_UTR_variant|non_coding_transcript_exon_variant|intron_variant|non_coding_transcript_variant|upstream_gene_variant|downstream_gene_variant|TFBS_ablation|TFBS_amplification|TF_binding_site_variant|regulatory_region_ablation|regulatory_region_amplification|feature_elongation|feature_elongation|regulatory_region_variant|feature_truncation|intergenic_variant|NMD_transcript_variant/' "$samplename"_basic_filtered.vep > "$samplename"_basic_filtered_consequence_filtered.vep

	#STEP 6: Filter "$samplename"_consequence_filtered.vep using gawk on gnomAD values (if they exist).
	gawk 'BEGIN {numericcount=0;count=0;gnomadcount=0;} {if($0 ~ /^#|^##/) {print;count++} else {l=split($14,a,";"); for(x=1;x<=l;x++) {m=split(a[x],b,"="); if(b[1] ~ /AF|exomes/) {gnomadcount++; numericfield=b[2]; decimalnot=sprintf("%.14f",numericfield); if(decimalnot ~ /^-?[0-9]*([.][0-9]+)?$/) {if(decimalnot > numericcount) {numericcount=decimalnot}}}}} if(numericcount < 0.001 && count < 1 && gnomadcount > 0) {print} if(numericcount == 0 && count < 1 && gnomadcount == 0) {print} numericcount=0;count=0;decimalnot=0;gnomadcount=0;}' "$samplename"_basic_filtered_consequence_filtered.vep > "$samplename"_basic_filtered_consequence_filtered_gnomAD_filtered_0001.vep
done
