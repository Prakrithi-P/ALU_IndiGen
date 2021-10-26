# ALU_IndiGen
An Indian population specific Alu insertion map 

## split bed file in bins for correlation analysis
sort -k1,1 -k2,2n introns_hg38_pos.bed | mergeBed > Introns_hg38_pos.bed #similary for existing gene annotations
#**intron density per 1MB/chr** #similarly for genic alu density
../bedops/bin/bedmap --echo --echo-overlap-size --fraction-map 1 --delim '\t' windows2.bed Introns_hg38_pos.bed > windows_with_intron_lengths_merged.bed
awk ‘{print $4}’ windows_with_intron_lengths_merged.bed > to_add
awk 'BEGIN {FS=OFS=";"} {sum=0;n=0
for(i=1;i<=NF;i++)
{sum+=$i;++n} print sum 
}' to_add > sum_introns
paste sum_introns windows_with_intron_lengths_merged.bed
g<-c%>% group_by(chr) %>% summarise(mean(Introns),mean(Alus))

##GC content calculation
python3 ./GC_analysis.py -i /home/mmlab33/TagTaste/Homo_sapiens_assembly38.fasta -o gc -w 1000000 -s 1000000

#correlation
d<-read.csv("chr1", sep="\t",header=T)
ggscatter(d,x="GC.Content...",y="Alu.Density",add="reg.line",conf.int=F, cor.coef=T, cor.method = "pearson",title = "Alu Density ~ Intron Density", xlab="Introns(bp)",ylab = "Alu Insertions (bp)")

## chr wise density
b<-read.csv("closest_dist_bins_count_n",sep="\t",header=T)
d_s<-split(b$n,b$binned)
binsum<-lapply(d_s,sum)
write.table(binsum, "binsum", sep="\t", row.names = F, quote = F)

b$binned<-cut(b$kb,9)
d_s<-split(b$n,b$binned)
binsum<-lapply(d_s,sum)
write.table(binsum, "binsum", sep="\t", row.names = F, quote = F)

Count entries within bins (per 10mb of chr)
library(matrixStats)
dd<-as.numeric(scan("D:/ALU_INDIGEN/1000g/priv_genes", character(), quote=""))
bx<-as.numeric(scan("D:/ALU_INDIGEN/1000g/brks", character(), quote=""))
binCounts(dd, idxs = NULL,bx,right=F)


#Venn diagram showing Alus in IndiGen, SGDP and HGSVC
draw.triple.venn(area1 = 22098,                          # common
                 area2 = 12044,
                 area3 = 5102,
                 n12 =3560 ,
                 n23 =4394 ,
                 n13 = 3505,
                 n123 = 3288,
                 fill = c("grey", "#ff8080" ,"orange"),col=c("black","black","black"),fontfamily = "sans",
                 cex=3,
                 cat.cex=1.5,cat.fontfamily = "sans", cat.fontface = "bold",
                 print.mode=c("raw"),
                 category = c("IndiGen","1000 Genomes", "SAS(1000 Genomes)"))


Pca
p3<-ggplot(d3)+ geom_point(aes(x=PC1,y=PC2,color=Population))+scale_color_manual(values=c("violet","yellow","green","blue","black","orange"))+ggtitle("Top 1% insertions (37) showing high Fst")+xlab("PC1 (25.1%)")+ylab("PC2 (19.1%)")+theme(text=element_text(size=16,  family="sans",face="bold"))+ theme(axis.text.y=element_text(size=14),axis.text.x =element_text(size=14))+theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_blank())+guides(color= guide_legend(override.aes = list(size=5)))



