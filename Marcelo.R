

net=read.table("/Users/fla./Desktop/Marcelo/test_detect_res7_max.tsv",header = T)
net$frag1=paste(net$chrom1,net$start1,net$end1,sep="_")
net$frag2=paste(net$chrom2,net$start2,net$end2,sep="_")
head(net)
write.table(net,"/Users/fla./Desktop/Marcelo/net_hic_chromosight_def.txt",col.names = T,row.names = F,quote=F,sep="\t")


#chromosight to washu
net$frag1_washu=paste(net$chrom1,net$start1,net$end1,sep=",")
net$frag2_washu=paste(net$chrom2,net$start2,net$end2,sep=",")
x=net[,c(16,17,11)]
x$frag1_washu=paste0("chr",x$frag1_washu);x$frag2_washu=paste0("chr",x$frag2_washu)
head(x)
write.table(x,"/Users/fla./Desktop/Marcelo/net_hic_chromosight_max_WashU.txt",col.names = F,row.names = F,quote=F,sep="\t")


#gint to washu
hic=read.table("/Users/fla./Desktop/Marcelo/Hug2017_3-4h_all_640k_GInteractions.txt.tsv",header = F)
hic$frag1=paste(hic$V1,hic$V2,hic$V3,sep=",")
hic=hic[-c(1:3)]
hic$frag2=paste(hic$V4,hic$V5,hic$V6,sep=",")
hic=hic[-c(1:3)]
hic$frag1=paste0("chr",hic$frag1);hic$frag2=paste0("chr",hic$frag2)
hic=hic[,c(2,3,1)]
head(hic)
write.table(hic,"/Users/fla./Desktop/Marcelo/hic_Hug2017_3-4h_all_640k_WashU.txt",col.names = F,row.names = F,quote=F,sep="\t")


write.table(hic[1:100000,],"/Users/fla./Desktop/Marcelo/hic_test.txt",col.names = F,row.names = F,quote=F,sep="\t")


#bed to flybase bed
file=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/GSM409073_suHw_NO_1_fdr_consecutive_peak_filtered.bed",skip = 2)
file=file[,1:3]
file$frag=paste(file$V1,file$V2,sep=":")
file$frag2=paste(file$frag,file$V3,sep="..")
file=file[,c(5)]
write.table(file,"/Users/fla./Desktop/Marcelo/Chip-Seq/GSM409073_suHw_NO_1_fdr_consecutive_peak_filtered_pourflybase.bed",col.names = F,row.names = F,quote=F,sep="\t")

#flybase bed to bed
file=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/GSM409073_suHw_NO_1_fdr_consecutive_peak_filtered_dm6_flybase.tsv",sep="\t")
file=file[,1:2]
file=file[grep("\\..",file$V2),]
file=file[-grep("\\?",file$V2),]
file$V2=paste0("chr",file$V2)

file$V2=gsub(":","\t",file$V2)
file$V2=gsub("\\..","\t",file$V2)
head(file)
write.table(file$V2,"/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/suHw_dm6.bed",col.names = F,row.names = F,quote=F,sep="\t")




#chas
library(devtools)
devtools::install_bitbucket("eraineri/chaser",  build_vignettes=TRUE)
library(chaser)
browseVignettes(package="chaser")

###chip
#Negre
beaf=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/BEAF_dm6.bed.txt")
cp190=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CP190_dm6.bed")
ctcf_c=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CTCF_C_dm6.bed")
ctcf_n=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CTCF_N_dm6.bed")
gaf=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/GAF_dm6.bed")
mdj4=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/MDJ4_dm6.bed")
suhw=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/suHw_dm6.bed")


#Maksimento
pita=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Pita_S2_WT_Maksimenko_dm6_peaks.narrowPeak")
zipic=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/ZIPIC_S2_WT_Maksimenko_dm6_peaks.narrowPeak")

#Li no input
# chromator=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Chromator_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# dref=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/DREF_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# fs1h=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Fs1h_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# l3mbt=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/L3mbt_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# rad21=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/RAD21_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# z4=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Z4_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")
# cbp=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CBP_Kc167_WT_Li_dm6_noinput_peaks.narrowPeak")



#Li input
chromator=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Chromator_Kc167_WT_Li_dm6_peaks.narrowPeak")
dref=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/DREF_Kc167_WT_Li_dm6_peaks.narrowPeak")
fs1h=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Fs1h_Kc167_WT_Li_dm6_peaks.narrowPeak")
l3mbt=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/L3mbt_Kc167_WT_Li_dm6_peaks.narrowPeak")
rad21=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/RAD21_Kc167_WT_Li_dm6_peaks.narrowPeak")
z4=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Z4_Kc167_WT_Li_dm6_peaks.narrowPeak")
cbp=read.table("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CBP_Kc167_WT_Li_dm6_peaks.narrowPeak")

### Load networks

net_def=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_def$chrom1=paste0("chr",net_def$chrom1);net_def$chrom2=paste0("chr",net_def$chrom2)

net_max=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_hic_chromosight_max.txt",header = T)
net_max$chrom1=paste0("chr",net_max$chrom1);net_max$chrom2=paste0("chr",net_max$chrom2)

net_opt=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_hic_chromosight_opt.txt",header = T)
net_opt$chrom1=paste0("chr",net_opt$chrom1);net_opt$chrom2=paste0("chr",net_opt$chrom2)

library(chaser)
net_def <- make_chromnet(net_def[,1:6],missing=NA)
net_max <- make_chromnet(net_max[,1:6],missing=NA)
net_opt <- make_chromnet(net_opt[,1:6],missing=NA)

#mitotic
net_mito=read.table("/Users/fla./Desktop/Marcelo/chromosight/mitotic_detect_res7.tsv",header = T)
net_mito$chrom1=paste0("chr",net_mito$chrom1);net_mito$chrom2=paste0("chr",net_mito$chrom2)
net_mito <- make_chromnet(net_mito[,1:6],missing=NA)

#nc12
net_c12=read.table("/Users/fla./Desktop/Marcelo/chromosight/nc12_detect_res7.tsv",header = T)
net_c12$chrom1=paste0("chr",net_c12$chrom1);net_c12$chrom2=paste0("chr",net_c12$chrom2)
net_c12 <- make_chromnet(net_c12[,1:6],missing=NA)

#loop distances
net_c12_short <- make_chromnet(net_c12_short[,1:6],missing=NA)
net_c12_long <- make_chromnet(net_c12_long[,1:6],missing=NA)

#nc12 TAD borders - Hug et al 
net_c12_TAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc12_TADborders_chromosight.txt",header = T)
net_c12_TAD <- make_chromnet(net_c12_TAD[,1:6],missing=NA)
net_c12_TADnoTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_c12_TAD-noTADborders_chromosight.txt",header = T)
net_c12_TADnoTAD <- make_chromnet(net_c12_TADnoTAD[,1:6],missing=NA)
net_c12_noTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc12_noTADborders_chromosight.txt",header = T)
net_c12_noTAD <- make_chromnet(net_c12_noTAD[,1:6],missing=NA)

#nc13
net_c13=read.table("/Users/fla./Desktop/Marcelo/chromosight/nc13_detect_res7.tsv",header = T)
net_c13$chrom1=paste0("chr",net_c13$chrom1);net_c13$chrom2=paste0("chr",net_c13$chrom2)
net_c13 <- make_chromnet(net_c13[,1:6],missing=NA)

#loop distances
net_c13_short <- make_chromnet(net_c13_short[,1:6],missing=NA)
net_c13_long <- make_chromnet(net_c13_long[,1:6],missing=NA)

#nc13 TAD borders - Hug et al 
net_c13_TAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc13_TADborders_chromosight.txt",header = T)
net_c13_TAD <- make_chromnet(net_c13_TAD[,1:6],missing=NA)
net_c13_TADnoTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_c13_TAD-noTADborders_chromosight.txt",header = T)
net_c13_TADnoTAD <- make_chromnet(net_c13_TADnoTAD[,1:6],missing=NA)
net_c13_noTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc13_noTADborders_chromosight.txt",header = T)
net_c13_noTAD <- make_chromnet(net_c13_noTAD[,1:6],missing=NA)

#nc14
net_c14=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_nc14_chromosight_def.tsv",header = T)
net_c14$chrom1=paste0("chr",net_c14$chrom1);net_c14$chrom2=paste0("chr",net_c14$chrom2)
net_c14 <- make_chromnet(net_c14[,1:6],missing=NA)


#TAD net 
net_c14_TAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/TAD_borders_networks/net_nc14_TADborders_chromosight_def.txt",header = T)
net_c14_TAD_def$chrom1=paste0("chr",net_c14_TAD_def$chrom1);net_c14_TAD_def$chrom2=paste0("chr",net_c14_TAD_def$chrom2)
net_c14_TAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_TADborders_chromosight_max.txt",header = T)
net_c14_TAD_max$chrom1=paste0("chr",net_c14_TAD_max$chrom1);net_c14_TAD_max$chrom2=paste0("chr",net_c14_TAD_max$chrom2)
net_c14_TAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_TADborders_chromosight_opt.txt",header = T)
net_c14_TAD_opt$chrom1=paste0("chr",net_c14_TAD_opt$chrom1);net_c14_TAD_opt$chrom2=paste0("chr",net_c14_TAD_opt$chrom2)
net_c14_TAD_def <- make_chromnet(net_c14_TAD_def[,1:6],missing=NA)
net_c14_TAD_max <- make_chromnet(net_c14_TAD_max[,1:6],missing=NA)
net_c14_TAD_opt <- make_chromnet(net_c14_TAD_opt[,1:6],missing=NA)

#no TAD net 
net_c14_noTAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/TAD_borders_networks/net_nc14_noTADborders_chromosight_def.txt",header = T)
net_c14_noTAD_def$chrom1=paste0("chr",net_c14_noTAD_def$chrom1);net_c14_noTAD_def$chrom2=paste0("chr",net_c14_noTAD_def$chrom2)
net_c14_noTAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_noTADborders_chromosight_max.txt",header = T)
net_c14_noTAD_max$chrom1=paste0("chr",net_c14_noTAD_max$chrom1);net_c14_noTAD_max$chrom2=paste0("chr",net_c14_noTAD_max$chrom2)
net_c14_noTAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_noTADborders_chromosight_opt.txt",header = T)
net_c14_noTAD_opt$chrom1=paste0("chr",net_c14_noTAD_opt$chrom1);net_c14_noTAD_opt$chrom2=paste0("chr",net_c14_noTAD_opt$chrom2)
net_c14_noTAD_def <- make_chromnet(net_c14_noTAD_def[,1:6],missing=NA)
net_c14_noTAD_max <- make_chromnet(net_c14_noTAD_max[,1:6],missing=NA)
net_c14_noTAD_opt <- make_chromnet(net_c14_noTAD_opt[,1:6],missing=NA)

#TAD/noTAD
net_c14_TADnoTAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_TAD-noTADborders_chromosight_def.txt",header = T)
net_c14_TADnoTAD_def$chrom1=paste0("chr",net_c14_TADnoTAD_def$chrom1);net_c14_TADnoTAD_def$chrom2=paste0("chr",net_c14_TADnoTAD_def$chrom2)
net_c14_TADnoTAD_def <- make_chromnet(net_c14_TADnoTAD_def[,1:6],missing=NA)

net_c14_TADnoTAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_TAD-noTADborders_chromosight_max.txt",header = T)
net_c14_TADnoTAD_max$chrom1=paste0("chr",net_c14_TADnoTAD_max$chrom1);net_c14_TADnoTAD_max$chrom2=paste0("chr",net_c14_TADnoTAD_max$chrom2)
net_c14_TADnoTAD_max <- make_chromnet(net_c14_TADnoTAD_max[,1:6],missing=NA)

net_c14_TADnoTAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_nc14_TAD-noTADborders_chromosight_opt.txt",header = T)
net_c14_TADnoTAD_opt$chrom1=paste0("chr",net_c14_TADnoTAD_opt$chrom1);net_c14_TADnoTAD_opt$chrom2=paste0("chr",net_c14_TADnoTAD_opt$chrom2)
net_c14_TADnoTAD_opt <- make_chromnet(net_c14_TADnoTAD_opt[,1:6],missing=NA)

#loop distances
net_c14_short=net_c14[which(net_c14$dist<90000),]
net_c14_long=net_c14[which(net_c14$dist>=90000),]
net_c14_short <- make_chromnet(net_c14_short[,1:6],missing=NA)
net_c14_long <- make_chromnet(net_c14_long[,1:6],missing=NA)

#nc14 TAD borders - Hug et al 
net_c14_TAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc14_TADborders_chromosight.txt",header = T)
net_c14_TAD <- make_chromnet(net_c14_TAD[,1:6],missing=NA)
net_c14_TADnoTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_c14_TAD-noTADborders_chromosight.txt",header = T)
net_c14_TADnoTAD <- make_chromnet(net_c14_TADnoTAD[,1:6],missing=NA)
net_c14_noTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_nc14_noTADborders_chromosight.txt",header = T)
net_c14_noTAD <- make_chromnet(net_c14_noTAD[,1:6],missing=NA)

#nc14 Intra/Inter TADs
net_c14_same_TAD=read.table("/Users/fla./Desktop/New/Intra_Inter_TADs/net_c14_same_TAD.txt",header = T)
net_c14_same_TAD <- make_chromnet(net_c14_same_TAD[,1:6],missing=NA)

net_c14_diff_TAD=read.table("/Users/fla./Desktop/New/Intra_Inter_TADs/net_c14_diff_TAD.txt",header = T)
net_c14_diff_TAD <- make_chromnet(net_c14_diff_TAD[,1:6],missing=NA)



#3-4h
net_34h=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_34h$chrom1=paste0("chr",net_34h$chrom1);net_34h$chrom2=paste0("chr",net_34h$chrom2)
net_34h <- make_chromnet(net_34h[,1:6],missing=NA)

#3-4h TAD net 
net_34h_TAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/TAD_borders_networks/net_3-4h_TADborders_chromosight_def.txt",header = T)
net_34h_TAD_def$chrom1=paste0("chr",net_34h_TAD_def$chrom1);net_34h_TAD_def$chrom2=paste0("chr",net_34h_TAD_def$chrom2)
net_34h_TAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_3-4h_TADborders_chromosight_max.txt",header = T)
net_34h_TAD_max$chrom1=paste0("chr",net_34h_TAD_max$chrom1);net_34h_TAD_max$chrom2=paste0("chr",net_34h_TAD_max$chrom2)
net_34h_TAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_3-4h_TADborders_chromosight_opt.txt",header = T)
net_34h_TAD_opt$chrom1=paste0("chr",net_34h_TAD_opt$chrom1);net_34h_TAD_opt$chrom2=paste0("chr",net_34h_TAD_opt$chrom2)
net_34h_TAD_def <- make_chromnet(net_34h_TAD_def[,1:6],missing=NA)
net_34h_TAD_max <- make_chromnet(net_34h_TAD_max[,1:6],missing=NA)
net_34h_TAD_opt <- make_chromnet(net_34h_TAD_opt[,1:6],missing=NA)

#3-4h no TAD net 
net_34h_noTAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/TAD_borders_networks/net_3-4h_noTADborders_chromosight_def.txt",header = T)
net_34h_noTAD_def$chrom1=paste0("chr",net_34h_noTAD_def$chrom1);net_34h_noTAD_def$chrom2=paste0("chr",net_34h_noTAD_def$chrom2)
net_34h_noTAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_3-4h_noTADborders_chromosight_max.txt",header = T)
net_34h_noTAD_max$chrom1=paste0("chr",net_34h_noTAD_max$chrom1);net_34h_noTAD_max$chrom2=paste0("chr",net_34h_noTAD_max$chrom2)
net_34h_noTAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_3-4h_noTADborders_chromosight_opt.txt",header = T)
net_34h_noTAD_opt$chrom1=paste0("chr",net_34h_noTAD_opt$chrom1);net_34h_noTAD_opt$chrom2=paste0("chr",net_34h_noTAD_opt$chrom2)
net_34h_noTAD_def <- make_chromnet(net_34h_noTAD_def[,1:6],missing=NA)
net_34h_noTAD_max <- make_chromnet(net_34h_noTAD_max[,1:6],missing=NA)
net_34h_noTAD_opt <- make_chromnet(net_34h_noTAD_opt[,1:6],missing=NA)

#3-4h TAD/noTAD
net_34h_TADnoTAD_def=read.table("/Users/fla./Desktop/Marcelo/TADcalling/TAD_borders_networks/net_34h_TAD-noTADborders_chromosight_def.txt",header = T)
net_34h_TADnoTAD_def$chrom1=paste0("chr",net_34h_TADnoTAD_def$chrom1);net_34h_TADnoTAD_def$chrom2=paste0("chr",net_34h_TADnoTAD_def$chrom2)
net_34h_TADnoTAD_def <- make_chromnet(net_34h_TADnoTAD_def[,1:6],missing=NA)

net_34h_TADnoTAD_max=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_34h_TAD-noTADborders_chromosight_max.txt",header = T)
net_34h_TADnoTAD_max$chrom1=paste0("chr",net_34h_TADnoTAD_max$chrom1);net_34h_TADnoTAD_max$chrom2=paste0("chr",net_34h_TADnoTAD_max$chrom2)
net_34h_TADnoTAD_max <- make_chromnet(net_34h_TADnoTAD_max[,1:6],missing=NA)

net_34h_TADnoTAD_opt=read.table("/Users/fla./Desktop/Marcelo/TADcalling/net_34h_TAD-noTADborders_chromosight_opt.txt",header = T)
net_34h_TADnoTAD_opt$chrom1=paste0("chr",net_34h_TADnoTAD_opt$chrom1);net_34h_TADnoTAD_opt$chrom2=paste0("chr",net_34h_TADnoTAD_opt$chrom2)
net_34h_TADnoTAD_opt <- make_chromnet(net_34h_TADnoTAD_opt[,1:6],missing=NA)

#3-4h TAD borders - Hug et al 
net_34h_TAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_n34h_TADborders_chromosight.txt",header = T)
net_34h_TAD <- make_chromnet(net_34h_TAD[,1:6],missing=NA)
net_34h_TADnoTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_34h_TAD-noTADborders_chromosight.txt",header = T)
net_34h_TADnoTAD <- make_chromnet(net_34h_TADnoTAD[,1:6],missing=NA)
net_34h_noTAD=read.table("/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_n34h_noTADborders_chromosight.txt",header = T)
net_34h_noTAD <- make_chromnet(net_34h_noTAD[,1:6],missing=NA)

#3-4h  loop distances
net_34h=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_34h$chrom1=paste0("chr",net_34h$chrom1);net_34h$chrom2=paste0("chr",net_34h$chrom2)
net_34h$dist=abs(net_34h$start2-net_34h$start1)
net_34h_short=net_34h[which(net_34h$dist<270000),]
net_34h_long=net_34h[which(net_34h$dist>=270000),]
net_34h_short <- make_chromnet(net_34h_short[,1:6],missing=NA)
net_34h_long <- make_chromnet(net_34h_long[,1:6],missing=NA)

#34h Intra/Inter TADs
net_34h_same_TAD=read.table("/Users/fla./Desktop/New/Intra_Inter_TADs/net_34h_same_TAD.txt",header = T)
net_34h_same_TAD <- make_chromnet(net_34h_same_TAD[,1:6],missing=NA)
net_34h_diff_TAD=read.table("/Users/fla./Desktop/New/Intra_Inter_TADs/net_34h_diff_TAD.txt",header = T)
net_34h_diff_TAD <- make_chromnet(net_34h_diff_TAD[,1:6],missing=NA)

# Zelda depleted
net_zld=read.table("/Users/fla./Desktop/Marcelo/chromosight/sh_zld_detect_res7.tsv",header = T)
net_zld$chrom1=paste0("chr",net_zld$chrom1);net_zld$chrom2=paste0("chr",net_zld$chrom2)
net_zld <- make_chromnet(net_zld[,1:6],missing=NA)

#triploi
net_trip=read.table("/Users/fla./Desktop/Marcelo/chromosight/triploid_detect_res7.tsv",header = T)
net_trip$chrom1=paste0("chr",net_trip$chrom1);net_trip$chrom2=paste0("chr",net_trip$chrom2)
net_trip <- make_chromnet(net_trip[,1:6],missing=NA)

#Select netwok
net=net_trip


#Load the features
#Negre
net <- load_features(net, beaf, type="bed3",missingv=0,auxfun="mean",featnames = "BEAF-32")
net <- load_features(net, cbp, type="bed3",missingv=0,auxfun="mean",featnames = "CBP")
net <- load_features(net, chromator, type="bed3",missingv=0,auxfun="mean",featnames = "Chromator")

net <- load_features(net, cp190, type="bed3",missingv=0,auxfun="mean",featnames = "CP190")
net <- load_features(net, ctcf_c, type="bed3",missingv=0,auxfun="mean",featnames = "CTCF_C")
net <- load_features(net, ctcf_n, type="bed3",missingv=0,auxfun="mean",featnames = "CTCF_N")
net <- load_features(net, dref, type="bed3",missingv=0,auxfun="mean",featnames = "DREF")
net <- load_features(net, fs1h, type="bed3",missingv=0,auxfun="mean",featnames = "Fs1h")

net <- load_features(net, gaf, type="bed3",missingv=0,auxfun="mean",featnames = "GAF")
net <- load_features(net, l3mbt, type="bed3",missingv=0,auxfun="mean",featnames = "L3mbt")

net <- load_features(net, mdj4, type="bed3",missingv=0,auxfun="mean",featnames = "MDJ4")
net <- load_features(net, pita, type="bed3",missingv=0,auxfun="mean",featnames = "Pita")
net <- load_features(net, rad21, type="bed3",missingv=0,auxfun="mean",featnames = "RAD21")

net <- load_features(net, suhw, type="bed3",missingv=0,auxfun="mean",featnames = "suHw")
#Maksimento
net <- load_features(net, z4, type="bed3",missingv=0,auxfun="mean",featnames = "Z4")

net <- load_features(net, zipic, type="bed3",missingv=0,auxfun="mean",featnames = "ZIPIC")
#Li



#Select the TAD borders
#mitotic
# net <- load_features(net, TAD_mitotic, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")
# #nc12
# net <- load_features(net, TAD_c12, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")
# #nc13
# net <- load_features(net, TAD_c13, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")
# #nc14
# net <- load_features(net, TAD_c14, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")
# #3-4h
# net <- load_features(net, TAD_34h, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")
# 

crosschas_net <- chas(net,"crosschas")
pheatmap(crosschas_net,display_numbers = T,number_color = "black", angle_col=45,fontsize_number = 11)

### Cross chas
library(RColorBrewer)
library(pheatmap)
breaksList = seq(0, .25, by = .001)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
pheatmap(crosschas_net_short,display_numbers = T,number_color = "black",color=color, angle_col=45,fontsize_number = 11,breaks =breaksList)



random_hic_net=randomize(net,50,dist.match=T)
crosschas_list=list()
for (i in 1:length(random_hic_net)){
  crosschas_list[[i]]=chas(random_hic_net[[i]],"crosschas")
}

#order random crosschas distribution by features pair
crosschas_list_ordered=list()
  for(i in 1:length(crosschas_list)){
    for(k in 1:ncol(net$features)){
      for (j in 1:ncol(net$features)){
        crosschas_list_ordered[[paste(colnames(crosschas_list[[i]])[k],rownames(crosschas_list[[i]])[j],sep="_")]]=c(crosschas_list_ordered[[paste(colnames(crosschas_list[[i]])[1],rownames(crosschas_list[[i]])[2],sep="_")]],crosschas_list[[i]][j,k])
    }
  }
}

#order real crosschas by features pair
crosschas_ordered=list()
  for(k in 1:ncol(net$features)){
    for (j in 1:ncol(net$features)){
      crosschas_ordered[[paste(colnames(crosschas_net)[k],rownames(crosschas_net)[j],sep="_")]]=crosschas_net[j,k]
  }
}




crosschas_list_zscore=list()
for (i in names(crosschas_list_ordered)){
  crosschas_list_zscore[[i]]=(crosschas_ordered[[i]]-mean(crosschas_list_ordered[[i]]))/sd(crosschas_list_ordered[[i]])
}

crosschas_list_zscore=matrix(unlist(crosschas_list_zscore),nrow=ncol(net$features),ncol=ncol(net$features),dimnames=list(colnames(net$features),colnames(net$features)))

library(pheatmap)
breaksList = seq(0, .25, by = .001)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
pheatmap(crosschas_list_zscore,display_numbers = T,number_color = "black",color=color, angle_col=45,fontsize_number = 11, main="CrossChAs Z-Scores")#,breaks =breaksList)



### ChAs
chas_net <- chas(net,"chas")

random_hic_net=randomize(net,50,dist.match=T)
rm(chas_random_hic_net)
chas_random_hic_net=vector()
for (i in 1:length(random_hic_net)){
  chas_random_hic_net=c(chas_random_hic_net,chas(random_hic_net[[i]],"chas"))
}

value=chas_random_hic_net
sample=names(chas_random_hic_net)
# sample2=tstrsplit(sample,"_")[[4]]

rm(df_random_hic_net)
df_random_hic_net=data.frame()
df_random_hic_net=matrix(c(value,sample),length(chas_random_hic_net),ncol=2)
colnames(df_random_hic_net)=c("Value","Sample")
df_random_hic_net=as.data.frame(df_random_hic_net)
df_random_hic_net$Value=as.numeric(as.character(df_random_hic_net$Value))

value=chas_net
sample=names(chas_net)

rm(df_chas_hic_net)
df_chas_hic_net=data.frame()
df_chas_hic_net=matrix(c(value,sample),length(chas_net),ncol=2)
colnames(df_chas_hic_net)=c("Value","Sample")
df_chas_hic_net=as.data.frame(df_chas_hic_net)
df_chas_hic_net$Value=as.numeric(as.character(df_chas_hic_net$Value))

df_hic_net=rbind(df_random_hic_net,df_chas_hic_net)


#levels hic_net
all_color=c("#92C5DE","#4393C3","#2166AC","#FEE0D2","red","blue","cyan","black","pink","orange","purple","grey","gold3","orchid","hotpink","burlywood")#,"slateblue2")#,"red4")
all_color=rainbow(16)

color=c(rep(all_color,nrow(df_random_hic_net)/nrow(df_chas_hic_net)),rep("#CC0000",nrow(df_chas_hic_net)))
labels=vector()
labels=c(rep("",nrow(df_hic_net)-(length(unique(df_hic_net$Sample)))),as.character(df_hic_net$Sample[(nrow(df_hic_net)-length(unique(df_hic_net$Sample))+1):nrow(df_hic_net)]))

#Make the plot !
#df_hic_net$Sample=factor(df_hic_net$Sample, levels=names(chas_net))
  #df_mono$Genes <-factor(df_mono$Genes, levels=order)
  
# library(ggplot2)
# rm(p)
# p <- ggplot(data=df_hic_net, aes(x=Sample, y=Value, group=Sample),show.legend=F)
# p <- p + ylim(-.2,.4)
# p <- p + ylab("ChAs")
# p <- p + xlab("Features")
# p <- p + geom_jitter(colour=color,width = 0.1, size=2)
# p <- p + geom_text(aes(label=labels),hjust=0, vjust=0,size=3)
# p <- p + theme_classic()  # Black and white theme
# p <- p + ggtitle("Chromosight - nc13 - def")
# print(p)





####compute the z-score
test=df_hic_net
nb_features=length(unique(df_hic_net$Sample))
nrow(df_hic_net)/51
test_random=test[1:(nb_features*50),]
test_chas=test[-(1:(nb_features*50)),]
nrow(test_random)
nrow(test_chas)
mean_random <- aggregate(. ~ Sample, data = test_random, mean)
max_random <- aggregate(. ~ Sample, data = test_random, max)
min_random <- aggregate(. ~ Sample, data = test_random, min)
sd_random <- aggregate(. ~ Sample, data = test_random, sd)


rm(i,zscore)
zscore=vector()
# for(i in mean_random$Sample){
#   zscore=c(zscore,ifelse(test_chas[which(test_chas$Sample==i),1] > max_random[which(max_random$Sample==i),2],(test_chas[which(test_chas$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2],0))
#   zscore=ifelse(test_chas[which(test_chas$Sample==i),1] < min_random[which(min_random$Sample==i),2],(test_chas[which(test_chas$Sample==i),1]+mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2],zscore)
# }

# for(i in mean_random$Sample){
#   if (test_chas[which(test_chas$Sample==i),1] > max_random[which(max_random$Sample==i),2]){
#     zscore=c(zscore,(test_chas[which(test_chas$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
#   }else if (test_chas[which(test_chas$Sample==i),1] < min_random[which(min_random$Sample==i),2]){
#     zscore=c(zscore(test_chas[which(test_chas$Sample==i),1]+mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
#   }else{
#     zscore=c(zscore,0)
#   }}


for(i in mean_random$Sample){
    zscore=c(zscore,(test_chas[which(test_chas$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscore)=mean_random$Sample

b=barplot(zscore,col=all_color, ylab="z-score",main="Z-score - triploid - whole genome",cex.names=.7,ylim=c(-3.5,25),font=2,xaxt="n")
pos=ifelse(zscore>0,round(zscore,1),"")
neg=ifelse(zscore<0,round(zscore,1),"")
text(b, zscore, labels=pos, cex= .9,pos=3)
text(b, zscore, labels=neg, cex= .9,pos=1)
text(cex=1.1, x=b-.5, y=-8, names(zscore), xpd=T, srt=45)


#ChAs vs Abundance plot
library(colorspace)
#all_color=c("#92C5DE","#4393C3","#2166AC","#FEE0D2","red","blue","cyan","black","pink","orange","purple","grey","gold3","orchid","hotpink","burlywood")#,"red4","slateblue2")
nodes=vector()
net$features=as.data.frame(net$features)
for (i in 1:ncol(net$features)){
  nodes=c(nodes,print(length(net$features[which(net$features[,i]>0),i])))
}
names(nodes)=colnames(net$features)
# plot(nodes,chas_net,col=all_color,pch=16,ylim=c(-.2,0.4),xlab="Nodes",ylab="Chromatin Assortativity",main="Chromosight - 34h - same TADs")
# text(nodes, chas_net, labels=colnames(net$features), cex= 0.7,pos=3)
# text(nodes, chas_net, labels=colnames(net$features), cex= 0.7,pos=3,col=all_color)

# Z-Score vs Abundance
plot(nodes,zscore,col=all_color,pch=16,ylim=c(-3.5,20),xlab="Nodes",ylab="ChAs Z-Score",main="Chromosight - triploid - Whole Genome")
#text(nodes, zscore, labels=colnames(net$features), cex= 0.7,pos=3)
text(nodes, zscore, labels=colnames(net$features), cex= 0.7,pos=3,col=lighten(all_color,-0.4))
abline(h=0,lty=2,col="grey")




#add frag column on files
i="/Users/fla./Desktop/Marcelo/chromosight/mitotic_detect_res7.tsv"
test=read.table(i,header=T)
test$frag1=paste(paste0("chr",test$chrom1),test$start1,test$end1,sep="_")
test$frag2=paste(paste0("chr",test$chrom2),test$start2,test$end2,sep="_")
write.table(test,i,col.names=T,row.names=F,quote=F,sep="\t")




i="BEAF-32"
x=c(1:5)
y_zscore=c(zscore_mitotic_def[which(names(zscore_mitotic_def)==i)],
           zscore_nc12_def[which(names(zscore_nc12_def)==i)],
           zscore_nc13_def[which(names(zscore_nc13_def)==i)],
           zscore_nc14_def[which(names(zscore_nc14_def)==i)],
           zscore_34h_def[which(names(zscore_34h_def)==i)])


#ChAs evolution across stages
plot(x, y_zscore_z4, type = "b",col="red",xlab="Stages",ylab="z-score",lwd=2)
lines(x, y_zscore_chromator, type = "b", col = "orange",lwd=2)
lines(x, y_zscore_dref, type = "b", col = "purple",lwd=2)
lines(x, y_zscore_beaf32, type = "b", col = "darkgreen",lwd=2)
legend("topleft", legend = c("Z4","Chromator","DREF","BEAF-32"),
       col = c( "red", "orange","purple","darkgreen"),lty = 1,lwd=2)

plot(x, y_zscore_tad_borders, type = "b",col="orange",xlab="Stages",ylab="z-score",lwd=2,ylim=c(0,0.3171422))
legend("topleft", legend = c("TAD borders"),
       col = c("orange"),lty = 1,lwd=2)




#barplots
barplot=barplot(c(4600,991,6317,6546,3106,2488,3917,8682,3842,12115,3009,7105,7873,3383,6330,34217),
names=c("BEAF-32","CBP","Chromator","CP190","CTCF_C","CTCF_N","DREF","Fs1h","GAF","L3mbt","Mdj4","Pita","RAD21","suHw","Z4","ZIPIC"),
col=all_color,ylab="Peaks",
ylim=c(0,40000),xlab="Features",cex.names=0.8,main="Features density - Whole Genome")
text(x = barplot, y = c(4600,991,6317,6546,3106,2488,3917,8682,3842,12115,3009,7105,7873,3383,6330,34217), label = c(4600,991,6317,6546,3106,2488,3917,8682,3842,12115,3009,7105,7873,3383,6330,34217), pos = 3, cex = .8, col = "black")





#Peaks density  plot at TAD borders
net_c14_TAD_def
frag=unique(GRanges(c(net_c14_TAD_def$chrom1,net_c14_TAD_def$chrom2),IRanges(c(net_c14_TAD_def$start1,net_c14_TAD_def$start2),c(net_c14_TAD_def$end1,net_c14_TAD_def$end2))))

library(genomation)
all_color=c("#92C5DE","#4393C3","#2166AC","#FEE0D2","red","blue","cyan","black","pink","orange","purple","grey","gold3","orchid","hotpink","burlywood")#,"red4","slateblue2")

beaf=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/BEAF_dm6.bed.txt")
cp190=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CP190_dm6.bed")
ctcf_c=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CTCF_C_dm6.bed")
ctcf_n=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CTCF_N_dm6.bed")
gaf=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/GAF_dm6.bed")
mdj4=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/MDJ4_dm6.bed")
suhw=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/suHw_dm6.bed")
pita=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Pita_S2_WT_Maksimenko_dm6_peaks.narrowPeak")
zipic=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/ZIPIC_S2_WT_Maksimenko_dm6_peaks.narrowPeak")
chromator=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Chromator_Kc167_WT_Li_dm6_peaks.narrowPeak")
dref=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/DREF_Kc167_WT_Li_dm6_peaks.narrowPeak")
fs1h=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Fs1h_Kc167_WT_Li_dm6_peaks.narrowPeak")
l3mbt=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/L3mbt_Kc167_WT_Li_dm6_peaks.narrowPeak")
rad21=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/RAD21_Kc167_WT_Li_dm6_peaks.narrowPeak")
z4=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/Z4_Kc167_WT_Li_dm6_peaks.narrowPeak")
cbp=readBed("/Users/fla./Desktop/Marcelo/Chip-Seq/dm6 bed/CBP_Kc167_WT_Li_dm6_peaks.narrowPeak")



TAD_c12=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c12_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c12$V1=paste0("chr",TAD_c12$V1)
TAD_c12=GRanges(TAD_c12$V1,IRanges(TAD_c12$V2,TAD_c12$V3))
TAD_c13=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c13_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c13$V1=paste0("chr",TAD_c13$V1)
TAD_c13=GRanges(TAD_c13$V1,IRanges(TAD_c13$V2,TAD_c13$V3))
TAD_c14=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c14_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c14$V1=paste0("chr",TAD_c14$V1)
TAD_c14=GRanges(TAD_c14$V1,IRanges(TAD_c14$V2,TAD_c14$V3))
TAD_34h=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_3-4h_all_5k_TAD_fdr_def_boundaries.bed")
TAD_34h$V1=paste0("chr",TAD_34h$V1)
TAD_34h=GRanges(TAD_34h$V1,IRanges(TAD_34h$V2,TAD_34h$V3))

chip=list(beaf,cbp,chromator,cp190,ctcf_c,ctcf_n,dref,fs1h,gaf,l3mbt,mdj4,pita,rad21,suhw,z4,zipic)
vec=vector()
for (i in 1:length(chip)){
ol=findOverlaps(frag,chip[[i]])
vec=c(vec,length(unique(subjectHits(ol))))
}

names(vec)=c("BEAF-32","CBP","Chromator","CP190","CTCF_C","CTCF_N","DREF","Fs1h","GAF","L3mbt","MDJ4","Pita","RAD21","suHw","Z4","ZIPIC")
barplot=barplot(vec,col=all_color,ylim=c(0,70),xlab="Features",cex.names=0.8,ylab="Peaks",main="Features density - nc14 TAD borders in network\n n=38")
text(x = barplot, y = vec, label = vec, pos = 3, cex = 1, col = "black")



#loops density across development and resolutions
#plot1
x=c(1:10)
y_34h=c(2,4,22,130,502,1811,4702,5031,7189,24805)
y_c14=c(3,5,18,111,494,1551,2709,2153,3206,17567)
y_c13=c(1,3,27,127,409,935,1377,1524,1948,10979)
y_c12=c(0,3,23,113,364,760,1031,1107,1495,10041)
y_mitotic=c(0,6,36,105,248,398,415,456,656,6320)

plot(x, y_34h, type = "l",xlab="Resolution",ylab="Loops")
lines(x, y_c14, type = "l", col = "red")
lines(x, y_c13, type = "l", col = "green")
lines(x, y_c12, type = "l", col = "orange")
lines(x, y_mitotic, type = "l", col = "blue")

legend("topleft", legend = c("3-4h", "nc14", "nc13","nc12","mitotic"),
       col = c("black", "red", "green","orange","blue"),lty = 1)

#plot2
x=c(1:5)
y0=c(0,1,3,2,0)
y1=c(3,3,5,4,6)
y2=c(23,27,18,22,36)
y3=c(113,127,111,130,105)
y4=c(364,409,494,502,248)
y5=c(760,935,1551,1811,398)
y6=c(1031,1377,2709,4702,415)
y7_def=c(1107,1524,2153,5031,456)
y7_max=c(1495,1948,3206,7189,656)
y7_opt=c(10041,10979,17567,24805,6320)

plot(x, y7_opt, type = "l",xlab="State",ylab="Loops",ylim=c(0,25000))
lines(x, y1, type = "l", col = "red",ylim=c(0,25000))
lines(x, y2, type = "l", col = "green",ylim=c(0,25000))
lines(x, y3, type = "l", col = "orange",ylim=c(0,25000))
lines(x, y4, type = "l", col = "blue",ylim=c(0,25000))
lines(x, y5, type = "l", col = "purple",ylim=c(0,25000))
lines(x, y6, type = "l", col = "pink",ylim=c(0,25000))
lines(x, y7_def, type = "l", col = "brown",ylim=c(0,25000))
lines(x, y7_max, type = "l", col = "grey",ylim=c(0,25000))
lines(x, y0, type = "l", col = "cyan",ylim=c(0,25000))

legend("topleft", legend = c("res_5kb_opt", "res_5kb_max", "res_5kb_def","res_10kb","res_20kb","res_40kb","res_80kb","res_160kb","res_320kb","res_640kb"),
       col = c("black", "grey", "brown","pink","purple","blue","orange","green","red","cyan"),lty = 1)




#TAD borders/noTAD borders networks
#load TAD borders to GRanges
TAD_mitotic=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_mitotic_all_5k_TAD_fdr_def_boundaries.bed")
TAD_mitotic$V1=paste0("chr",TAD_mitotic$V1)
TAD_c12=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c12_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c12$V1=paste0("chr",TAD_c12$V1)
TAD_c13=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c13_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c13$V1=paste0("chr",TAD_c13$V1)
TAD_c14=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_c14_all_5k_TAD_fdr_def_boundaries.bed")
TAD_c14$V1=paste0("chr",TAD_c14$V1)
TAD_34h=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_3-4h_all_5k_TAD_fdr_def_boundaries.bed")
TAD_34h$V1=paste0("chr",TAD_34h$V1)

#load networks
net_34h=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_34h$chrom1=paste0("chr",net_34h$chrom1);net_34h$chrom2=paste0("chr",net_34h$chrom2)
net_34h$frag1=paste(net_34h$chrom1,net_34h$start1,net_34h$end1,sep="_");net_34h$frag2=paste(net_34h$chrom2,net_34h$start2,net_34h$end2,sep="_")

frag_34h=GRanges(c(net_34h$chrom1,net_34h$chrom2),IRanges(c(net_34h$start1,net_34h$start2),c(net_34h$end1,net_34h$end2)))
ol=findOverlaps(frag_34h,TAD_34h)

#TAD borders
frag_34h_TAD=as.data.frame(unique(frag_34h[queryHits(ol)]))
frag_34h_TAD=paste(frag_34h_TAD$seqnames,frag_34h_TAD$start,frag_34h_TAD$end,sep="_")
net_34h_TAD=net_34h[which(net_34h$frag1 %in% frag_34h_TAD & net_34h$frag2 %in% frag_34h_TAD),]
write.table(net_34h_TAD,"/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_n34h_TADborders_chromosight.txt",col.names = T,row.names=F,sep="\t",quote=F)

#noTAD borders
frag_34h_noTAD=as.data.frame(unique(frag_34h[-queryHits(ol)]))
frag_34h_noTAD=paste(frag_34h_noTAD$seqnames,frag_34h_noTAD$start,frag_34h_noTAD$end,sep="_")
net_34h_noTAD=net_34h[which(net_34h$frag1 %in% frag_34h_noTAD & net_34h$frag2 %in% frag_34h_noTAD),]
write.table(net_34h_noTAD,"/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_n34h_noTADborders_chromosight.txt",col.names = T,row.names=F,sep="\t",quote=F)

#TAD borders / no TAD borders
net_34h_TADnoTAD=net_34h[which(net_34h$frag1 %in% frag_34h_TAD & net_34h$frag2 %in% frag_34h_noTAD),]
net_34h_TADnoTAD=rbind(net_34h_TADnoTAD,net_34h[which(net_34h$frag1 %in% frag_34h_noTAD & net_34h$frag2 %in% frag_34h_TAD),])
write.table(net_34h_TADnoTAD,"/Users/fla./Desktop/New/TAD_borders_networks_Hug/net_34h_TAD-noTADborders_chromosight.txt",col.names = T,row.names=F,sep="\t",quote=F)


#numbers TAD borders across development

#Plot
x=c(1:5)
y0=c(861,890,930,925,987)

plot(x, y0, type = "b",col="red",xlab="Stage",ylab="TAD borders",ylim=c(800,1000),names=c("cc","cc","cc","cc","cc"))

legend("topleft", legend = c("res_5kb_opt", "res_5kb_max", "res_5kb_def","res_10kb","res_20kb","res_40kb","res_80kb","res_160kb","res_320kb","res_640kb"),
       col = c("black", "grey", "brown","pink","purple","blue","orange","green","red","cyan"),lty = 1)



## loop distances subnetworks
#nc12
net_c12=read.table("/Users/fla./Desktop/Marcelo/chromosight/nc12_detect_res7.tsv",header = T)
net_c12$chrom1=paste0("chr",net_c12$chrom1);net_c12$chrom2=paste0("chr",net_c12$chrom2)

net_c12$dist=abs(net_c12$start2-net_c12$start1)
plot(density(net_c12$dist),main="Loop distances distribution - nc12 - def",xlim=c(-100000,2000000))
median(net_c12$dist)
abline(v=250000,col="red")
text(354737, 8.351791e-06, labels="threshold = 250kb", cex= 1,pos=4,col="red")

net_c12_short=net_c12[which(net_c12$dist<250000),]
net_c12_long=net_c12[which(net_c12$dist>=250000),]

#nc13
net_c13=read.table("/Users/fla./Desktop/Marcelo/chromosight/nc13_detect_res7.tsv",header = T)
net_c13$chrom1=paste0("chr",net_c13$chrom1);net_c13$chrom2=paste0("chr",net_c13$chrom2)

net_c13$dist=abs(net_c13$start2-net_c13$start1)
plot(density(net_c13$dist),main="Loop distances distribution - nc13 - def",xlim=c(-100000,2000000))
median(net_c13$dist)
abline(v=250000,col="red")
text(354737, 8.351791e-06, labels="threshold = 250kb", cex= 1,pos=4,col="red")

net_c13_short=net_c13[which(net_c13$dist<250000),]
net_c13_long=net_c13[which(net_c13$dist>=250000),]

#nc14
net_c14=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_nc14_chromosight_def.tsv",header = T)
net_c14$chrom1=paste0("chr",net_c14$chrom1);net_c14$chrom2=paste0("chr",net_c14$chrom2)

net_c14$dist=abs(net_c14$start2-net_c14$start1)
plot(density(net_c14$dist),main="Loop distances distribution - nc14 - def",xlim=c(-100000,2000000))
median(net_c14$dist)
abline(v=250000,col="red")
text(354737, 5.351791e-06, labels="threshold = 250kb", cex= 1,pos=4,col="red")

net_c14_short=net_c14[which(net_c14$dist<250000),]
net_c14_long=net_c14[which(net_c14$dist>=250000),]

#3-4h
net_34h=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_34h$chrom1=paste0("chr",net_34h$chrom1);net_34h$chrom2=paste0("chr",net_34h$chrom2)

net_34h$dist=abs(net_34h$start2-net_34h$start1)
plot(density(net_34h$dist),main="Loop distances distribution - 3-4h - def",xlim=c(-100000,2000000))
median(net_34h$dist)
abline(v=250000,col="red")
text(354737, 1.351791e-06, labels="threshold = 250kb", cex= 1,pos=4,col="red")

net_34h_short=net_34h[which(net_34h$dist<250000),]
net_34h_long=net_34h[which(net_34h$dist>=250000),]


#intra / inter TADs
TAD_borders=read.table("/Users/fla./Desktop/Marcelo/consensus_boundaries_5kb.tsv",header = T)
TAD_c12=TAD_borders[which(!(is.na(TAD_borders$c12))),]
TAD_c13=TAD_borders[which(!(is.na(TAD_borders$c13))),]
TAD_c14=TAD_borders[which(!(is.na(TAD_borders$c14))),]

TAD_34h=TAD_borders[which(!(is.na(TAD_borders$X3.4h))),]

TAD_34h=read.table("/Users/fla./Desktop/Marcelo/TADcalling/hicFindTADs/Hug2017_3-4h_all_5k_TAD_fdr_def_boundaries.bed")
TAD_34h$V1=paste0("chr",TAD_34h$V1)
TAD_34h$mid=(TAD_34h$V2+TAD_34h$V3)/2

#load networks
net_34h=read.table("/Users/fla./Desktop/Marcelo/chromosight/net_3-4h_chromosight_def.txt",header = T)
net_34h$chrom1=paste0("chr",net_34h$chrom1);net_34h$chrom2=paste0("chr",net_34h$chrom2)
net_34h$mid1=(net_34h$start1+net_34h$end1)/2;net_34h$mid2=(net_34h$start2+net_34h$end2)/2


net_34h_chrX=net_34h[which(net_34h$chrom1=="chrX"),]
TAD_34h_chrX=TAD_34h[which(TAD_34h$V1=="chrX"),]

vec=vector();vec2=vector()
vec=rep(0,nrow(net_34h_chrX))
vec2=rep(0,nrow(net_34h_chrX))

for(j in 1:nrow(TAD_34h)){
  vec=ifelse(net_34h_chrX$mid1>TAD_34h$mid[j] & net_34h_chrX$mid1<TAD_34h$mid[j+1],j,vec)
  vec2=ifelse(net_34h_chrX$mid2>TAD_34h$mid[j] & net_34h_chrX$mid2<TAD_34h$mid[j+1],j,vec2)
  }

net_34h_chrX$TAD1=vec
net_34h_chrX$TAD2=vec2
table(net_34h_chrX$TAD1==net_34h_chrX$TAD2)



unique(net_34h$chrom1)

net_34h_same_TAD=rbind(net_34h_chr2L[which(net_34h_chr2L$TAD1==net_34h_chr2L$TAD2),],
                       net_34h_chr2R[which(net_34h_chr2R$TAD1==net_34h_chr2R$TAD2),],
                       net_34h_chr3L[which(net_34h_chr3L$TAD1==net_34h_chr3L$TAD2),],
                       net_34h_chr3R[which(net_34h_chr3R$TAD1==net_34h_chr3R$TAD2),],
                       net_34h_chr4[which(net_34h_chr4$TAD1==net_34h_chr4$TAD2),],
                       net_34h_chrX[which(net_34h_chrX$TAD1==net_34h_chrX$TAD2),])

net_34h_diff_TAD=rbind(net_34h_chr2L[which(net_34h_chr2L$TAD1!=net_34h_chr2L$TAD2),],
                       net_34h_chr2R[which(net_34h_chr2R$TAD1!=net_34h_chr2R$TAD2),],
                       net_34h_chr3L[which(net_34h_chr3L$TAD1!=net_34h_chr3L$TAD2),],
                       net_34h_chr3R[which(net_34h_chr3R$TAD1!=net_34h_chr3R$TAD2),],
                       net_34h_chr4[which(net_34h_chr4$TAD1!=net_34h_chr4$TAD2),],
                       net_34h_chrX[which(net_34h_chrX$TAD1!=net_34h_chrX$TAD2),])

table(net_34h_same_TAD$TAD1==net_34h_same_TAD$TAD2)
table(net_34h_diff_TAD$TAD1==net_34h_diff_TAD$TAD2)

write.table(net_34h_same_TAD,"/Users/fla./Desktop/New/Intra_Inter_TADs/net_34h_same_TAD.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(net_34h_diff_TAD,"/Users/fla./Desktop/New/Intra_Inter_TADs/net_34h_diff_TAD.txt",col.names=T,row.names=F,quote=F,sep="\t")

rm(net_c14_same_TAD,net_c14_diff_TAD)
rm(net_c14_chr2L,net_c14_chr2R,net_c14_chr3L,net_c14_chr3R,net_c14_chr4,net_c14_chrX)



