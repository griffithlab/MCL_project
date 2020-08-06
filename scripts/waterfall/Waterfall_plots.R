#Load GenVisR and data.table libraries.
library(GenVisR)
library(ggplot2)
library(data.table)
library(plyr)
library(grid)
library(gridExtra)

#Load old_tumor_only manual review final data.
old_tumor_only <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/old_tumor_only/final_all_old_tumor_only_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from old_tumor_only.
final_old_tumor_only <- old_tumor_only[,c('sample','gene','mutation')]

#Create Waterfall plot for final_old_tumor_only, and draw it.
final_old_tumor_only_waterfall <- Waterfall(final_old_tumor_only,recurrence = 0.20,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_old_tumor_only_waterfall)


#Load new_tumor_only manual review final data.
new_tumor_only <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/new_tumor_only/final_all_new_tumor_only_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from new_tumor_only.
final_new_tumor_only <- new_tumor_only[,c('sample','gene','mutation')]

#Create Waterfall plot for final_new_tumor_only, and draw it.
final_new_tumor_only_waterfall <- Waterfall(final_new_tumor_only,recurrence = 0.20,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_new_tumor_only_waterfall)


#Load somatic manual review final data.
somatic <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/somatic/final_all_somatic_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from somatic.
final_somatic <- somatic[,c('sample','gene','mutation')]

#Create Waterfall plot for final_somatic, and draw it.
final_somatic_waterfall <- Waterfall(final_somatic,recurrence = 0.20,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_somatic_waterfall)


#Load all (new_tumor_only, old_tumor_only, somatic) manual review final data.
all <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/all_variants/final_all_variants_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from all.
final_all <- all[,c('sample','gene','mutation')]

#Create Waterfall plot for final_all, and draw it.
final_all_waterfall <- Waterfall(final_all,recurrence = 0.20,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_all_waterfall)


#Load all (new_tumor_only, old_tumor_only, somatic) manual review final data.
forpresentation_all <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/all_variants/corrected_forpresentation_final_all_variants_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from all.
final_forpresentation_all <- forpresentation_all[,c('sample','gene','mutation')]

#Create Waterfall plot for final_all, and draw it.
final_forpresentation_all_waterfall <- Waterfall(final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_forpresentation_all_waterfall)


#Load all (new_tumor_only, old_tumor_only, somatic) manual review final data (FINAL).
final_forpresentation_all <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/all_variants/finalcolumns_corrected_final_all_variants_150reddygeneintersected.tsv",header=T,sep='\t')

#Grab just sample, gene, and mutation columns from all.
final_final_forpresentation_all <- final_forpresentation_all[,c('sample','gene','mutation')]

#Create Waterfall plot for final_all, and draw it.
final_final_forpresentation_all_waterfall <- Waterfall(final_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(final_final_forpresentation_all_waterfall)



#Add nanostring subtype and cr/pr or sd/pd information to patients on waterfall plot using bvlen_metadata_final.
#Pull bvlen_pt_clinical_data_nano_somatic_updated with pfs stats.xlsx into R.
bvlen_metadata_final <-readxl::read_xlsx("/Users/matthewmosior/Research/BV-Len/BV-Len clinical data PHI-free 8_23_2019 with pfs stats updated.xlsx")
#Remove last row and last column.
bvlen_metadata_final$`Date Locked 6/30/2019` <- bvlen_metadata_final$...19 <- bvlen_metadata_final$...20 <- NULL
bvlen_metadata_final <- bvlen_metadata_final[-nrow(bvlen_metadata_final),]
#Rename column headers to get rid of spaces.
bvlen_metadata_final <- rename(bvlen_metadata_final,c("ID #"="ID","Nanostring subtype"="nanostring_subtype","c-Myc positive"="c-myc_positive","Bcl-2 positive"="bcl-2_positive","Bcl-6 positive"="bcl-6_positive","No FISH Abnormalities"="no_FISH_abnormalities","Number of prior Tx"="number_of_prior_tx","Cycle to best response"="cycle_to_best_response","Time to Tx failure"="time_to_tx_failure","Best Response"="best_response","Event = 0 Censor = 1"="pfs_stat"))
#Add column to bvlen_metadata by adding TWBE... prefix to ID column.
bvlen_metadata_final$sample <- gsub("^","TWBE-BV_Len_",bvlen_metadata_final$ID)
#Change expired column ("N" -> "0" and "Y" -> "1").
bvlen_metadata_final$Expired <- gsub("N","0",bvlen_metadata_final$Expired)
bvlen_metadata_final$Expired <- gsub("Y","1",bvlen_metadata_final$Expired)
#Change pfs_stat column ("0" -> "1" and "1" -> "0").
bvlen_metadata_final$pfs_stat_recoded <- bvlen_metadata_final$pfs_stat
bvlen_metadata_final$pfs_stat_recoded[bvlen_metadata_final$pfs_stat_recoded==0] <- 3
bvlen_metadata_final$pfs_stat_recoded[bvlen_metadata_final$pfs_stat_recoded==1] <- 2
bvlen_metadata_final$pfs_stat_recoded[bvlen_metadata_final$pfs_stat_recoded==3] <- 1
bvlen_metadata_final$pfs_stat_recoded[bvlen_metadata_final$pfs_stat_recoded==2] <- 0

#Load all (new_tumor_only, old_tumor_only, somatic) manual review final data (FINAL).
final_forpresentation_all <- read.table("/Users/matthewmosior/Research/BV-Len/final_post_bam_readcount_correction/all_variants/finalcolumns_corrected_final_all_variants_150reddygeneintersected.tsv",header=T,sep='\t')

#Add nanostring_subtype and best_response to final_forpresentation_all.
true_final_forpresentation_all <- merge(final_forpresentation_all,bvlen_metadata_final[,c("nanostring_subtype","best_response","sample")],by.x="sample",by.y="sample")

#Simplify the sample name in true_final_forpresentation_all.
true_final_forpresentation_all$sample <- gsub(pattern="TWBE-BV_Len_",replacement="",true_final_forpresentation_all$sample)

#Simplify mutation field to get rid of stuff after comma.
true_final_forpresentation_all$mutation <- gsub(pattern=",.*",replacement="",true_final_forpresentation_all$mutation)

#Simplify best_response field to collapse completeresponse/partialresponse and progressivedisease/stabledisease (create new column for this). 
#true_final_forpresentation_all$best_response[true_final_forpresentation_all$best_response == "Complete Response"] <- "Responders"
#true_final_forpresentation_all$best_response[true_final_forpresentation_all$best_response == "Partial Response"] <- "Responders"
#true_final_forpresentation_all$best_response[true_final_forpresentation_all$best_response == "Progressive Disease"] <- "Non-Responders"
#true_final_forpresentation_all$best_response[true_final_forpresentation_all$best_response == "Stable Disease"] <- "Non-Responders"

true_final_forpresentation_all$twocat_best_response[true_final_forpresentation_all$best_response == "Complete Response"] <- "Responders"
true_final_forpresentation_all$twocat_best_response[true_final_forpresentation_all$best_response == "Partial Response"] <- "Responders"
true_final_forpresentation_all$twocat_best_response[true_final_forpresentation_all$best_response == "Progressive Disease"] <- "Non-Responders"
true_final_forpresentation_all$twocat_best_response[true_final_forpresentation_all$best_response == "Stable Disease"] <- "Non-Responders"

#Change best_response and nanostring_subtype column headers to Best Response and Nanostring Subtype.
colnames(true_final_forpresentation_all)[colnames(true_final_forpresentation_all)=="best_response"] <- "Four Category Best Response"
#colnames(true_final_forpresentation_all)[colnames(true_final_forpresentation_all)=="nanostring_subtype"] <- "Nanostring Subtype"
colnames(true_final_forpresentation_all)[colnames(true_final_forpresentation_all)=="Nanostring Subtype"] <- "DLBCL Subtype (nanostring)"
colnames(true_final_forpresentation_all)[colnames(true_final_forpresentation_all)=="twocat_best_response"] <- "Two Category Best Response"

#Copy true_final_forpresentation_all to true_final_forpresentation_all_test.
#true_final_forpresentation_all_test <- true_final_forpresentation_all[,c("sample","Nanostring Subtype","Four Category Best Response","Two Category Best Response")]
#true_final_forpresentation_all_test <- true_final_forpresentation_all[,c("sample","DLBCL Subtype (nanostring)","Four Category Best Response","Two Category Best Response")]
true_final_forpresentation_all_test <- true_final_forpresentation_all[,c("sample","DLBCL Subtype (nanostring)","Four Category Best Response")]

#Add dummy variables.
#true_final_forpresentation_all_test[nrow(true_final_forpresentation_all_test)+1,] <- c("002-002","unclassified","dummy1","dummy3")
#true_final_forpresentation_all_test[nrow(true_final_forpresentation_all_test)+1,] <- c("002-002","unclassified","dummy2","dummy4")

#Add clinical data to clinical_true_final_forpresentation_all.
#clinical_true_final_forpresentation_all <- Clinical(inputData=true_final_forpresentation_all[,c("sample","Nanostring Subtype","Best Response")])
clinical_true_final_forpresentation_all <- Clinical(inputData=true_final_forpresentation_all_test,clinicalLayer=list(guides(fill=guide_legend(ncol=2)),scale_fill_manual(values=c("GCB"="red","ABC"="brown","unclassified"="deeppink4","NA"="gray7","Complete Response"="darkgoldenrod1","Partial Response"="darkgoldenrod4","Stable Disease"="dodgerblue4","Progressive Disease"="deepskyblue"))))

#Grab just the sample, gene, and mutation columns for all.
true_final_final_forpresentation_all <- true_final_forpresentation_all[,c('sample','gene','mutation')]

#Set waterfall object.
#Final_Waterfall <- Waterfall(true_final_final_forpresentation_all,clinical=clinical_true_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"=c("missense_variant","inframe_deletion","frameshift_variant","stop_gained","start_lost","inframe_insertion","splice_donor_variant","splice_acceptor_variant"),"color"=c("blue","green","pink","yellow","turquoise","gold","orange","red")))
#Final_Waterfall <- Waterfall(true_final_final_forpresentation_all,clinical=clinical_true_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"=c("missense_variant","stop_gained","start_lost","inframe_insertion","inframe_deletion","frameshift_variant","splice_donor_variant","splice_acceptor_variant"),"color"=c("blue","cadetblue4","cadetblue","darkgreen","chartreuse3","chartreuse","orange","red")))
Final_Waterfall <- Waterfall(true_final_final_forpresentation_all,clinical=clinical_true_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"=c("frameshift_variant","splice_donor_variant","splice_acceptor_variant","stop_gained","start_lost","inframe_deletion","inframe_insertion","missense_variant"),"color"=c("firebrick1","darkorchid1","darkorchid4","darkgreen","darkolivegreen","cyan","cyan4","chartreuse")))

#Adjust the legend label colors (make dummy1 and dummy2 white to hide them).
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$grobs[[15]]$grobs[[1]]$grobs[[25]]$children$guide.label.titleGrob.3249$children$GRID.text.3247$gp$col <- "white"
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$grobs[[15]]$grobs[[1]]$grobs[[26]]$children$guide.label.titleGrob.3252$children$GRID.text.3250$gp$col <- "white" 

#Create subheaders for guide-box (legend).
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$grobs[[15]]$grobs[[1]]$grobs[[2]]$children$guide.title.titleGrob.3228$children$GRID.text.3226$label <- "Subtype      Response Group"

#Adjust grob widths for the firstRow, secondRow, and lastRow (check Final_Waterfall@Grob$grobs$secondRow$grobs[[2]]$widths for the field with the excess space).
Final_Waterfall@Grob$grobs$firstRow$grobs[[2]]$widths[4] <- unit(1,"cm")
Final_Waterfall@Grob$grobs$secondRow$grobs[[2]]$widths[4] <- unit(1,"cm")
Final_Waterfall@Grob$grobs$thirdRow$grobs[[2]]$widths[4] <- unit(1,"cm")
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$widths[4] <- unit(1,"cm")

#Add back #FFC0CBFF (Non-Responders color) to Final_Waterfall.
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$grobs[[6]]$children$geom_rect.rect.3199$gp$fill[45] <- "#FFC0CBFF"
Final_Waterfall@Grob$grobs$lastRow$grobs[[2]]$grobs[[6]]$children$geom_rect.rect.3199$gp$fill[46] <- "#FFC0CBFF"

#Draw waterfall with clinical annotations added.
drawPlot(Final_Waterfall)
#drawPlot(Waterfall(true_final_final_forpresentation_all,clinical=clinical_true_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue")))





#Grab just sample, gene, and mutation columns from all.
true_final_final_forpresentation_all <- true_final_forpresentation_all[,c('sample','gene','mutation')]

#Create Waterfall plot for final_all, and draw it.
true_final_final_forpresentation_all_waterfall <- Waterfall(true_final_final_forpresentation_all,recurrence = 0.13,mutationHierarchy = data.table("mutation"="missense_variant","color"="blue"))
drawPlot(true_final_final_forpresentation_all_waterfall)
