#RNA-seq 2-872 ; ChIp-seq 880-1590 ;个性化内容 1600+
#检验数据完整性&解压数据(fastq-dump)
md5sum 文件名
cat
md5sum *.gz
cat *.md5
#如果用的是自己测序，就是这个
#http://jvenn.toulouse.inra.fr/app/example.html剩下的丢在韦恩图网站里，可以检验是否一致，可以加*
gunzip *.gz
#实验室数据解压就好，会把后缀gz的解压好
#或者公共数据要改.sra，再用如下操作，用的是fastq-dump软件
fastq-dump --split-3 *.sra
#双端数据
fastq-dump *.sra
#单端数据



#除了特别的软件，默认所有软件都装在python3里了
vim copy.sh
#!/bin/sh
set -e
cp  /mnt/USB1/2021.11.16-11.26/LCZ-hos1-ChIPseq/* /mnt/data5/wjc/aaworkflow/LCZ-hos1-ChIPseq
echo "copy done"
nohup sh copy.sh > copynohup 2>&1 
md5sum *.gz
cat *.md5
#然后丢到韦恩图网站去看数据完整性
gunzip *.gz



conda activate python3
vim qc.sh
#!/bin/sh
set -e
mkdir QC
fastqc -o QC -f fastq *.fastq
echo "qc done"
nohup sh qc.sh > qcnohup 2>&1 




trim_galore --phred33 --fastqc -q 20 -o 01_cleanData_trimGalore_trim1 --paired ${id}_R1.fastq.gz ${id}_R2.fastq.gz -j 4
#输出到一个叫01_cleanData_trimGalore_trim1的文件夹，并进行fastqc



#所有软件都在python3环境中找到
#conda activate python3
#要么就在base和基础环境中找到,自己也可以配置环境
#conda activate base
#conda deactivate



#查看挂载情况
#lsblk- f
#lsblk
#挂载方法
#sudo mount /dev/sd某 想要的位置
#密码是lyh
#解除挂载
##umount /dev/ sd 某 想要的位置
#umount 想要的位置



#可以计算reads数（针对的是fq文件）
cat ${id}.fastq | grep -c '^+$'
#这个不一定有用，也可以直接查看质检报告
#这方法也不一定有效，也可以用行数除以4，
cat xxxx.fastqc | wc -l
结果除以4
#可以计算大致的测序深度，拟南芥0.118G=120M，水稻0.37G=370M，人3G，酿酒酵母12M=0.012G
#测序深度说的都是十进制文件，在电脑上的大小才是二进制
#1G=10亿，1亿=0.1G，1K=1000，10K=1万，100K=10万，1000K=100万=1M，1M=100万，10M=1000万，100M=1亿=0.1G，1000M=10亿=1G
#K千，M百万，G十亿
#合并两个fq文件使用cat命令，cat 1.fq 2.fq > xxx.fq



#构建参考基因组索引
#hisat2索引
#bowtie索引
#bowtie2索引
hisat2-build -p 4 genome.fa genome
bowtie-build --threads 2 xxxxx.fa xxxx
bowtie2-build  xxxx.fa xxxx



gffread Saccharomyces_cerevisiae.R64-1-1.52.gff3 -T -o jiaomu.gtf
#也可以使用UCSC的工具完成bed12的制作
gtfToGenePred xxx.gtf yyyy.gengpred
genePredToBed yyyy.genepred zzzz.bed



#取q20变为bam
samtools view -q 20 -bS ${id}.sam > ${id}_q20.bam
#有的比对结果不好，所以需要把比对质量低的去除掉，bS是指的是输入输出文件格式，Sb也可以
#sort&index
samtools sort ${id}_q20.bam -o ${id}_q20_s.bam
samtools index ${id}_q20_s.bam
#计算率
samtools view -F 4 -c ${id}.sam
samtools view -F 4 -c ${id}_q20.bam
samtools view -F 4 -c ${id}_q20_s.bam
过滤掉没有比对上的就是4，把剩下的统计出来，不管是单端还是双端，只要没比对上的就是4
samtools rmdup -s ${id}_q20_s.bam ${id}_q20_s_rmdup.bam
#chipseq才有
#俺的代码samtools rmdup -s WT_input_q20_s.bam WT_input_q20_s_rmdup.bam
#去除pcr重复，s对应单末端，默认对应双末端 S是将双末端当作单末端处理
#要记住一点的就是sam和bam文件是不分单双端的，你自己知道才行。
#这一步在sort之后和index之前就行




#比对得到sam
hisat2 -p 16 --dta --rna-strandness RF -x /mnt/data5/wjc/genome_scripts/genome/arb/ensemble/index/hisat2/hisat2_index -1 ${id}1.fastq -2 ${id}2.fastq -S ${id}.sam
hisat2 -p 8 -x /bios-analysis10/omics2021_post/21210700069/jiyinzuzuoye/zuoyeyuanshiwenjian/yeastgenome/yeastgenome/indexyeast -U weiba.fastq -S bidui.sam
#双端和单端--hisat2 ，rnastrandness RF链特异性，FR也可以，单链的用F或者R，有了这个不会对比对结果产生什么影响，但是会加上一列dx，单端选择R还是F有所区别
#但是不会对判断链特异性有任何影响，其实对啥都没影响，反正没影响就对了，只是在那一列有影响罢了，但是好像对可变剪切有所影响
#dta结果适合用于stringtie进行转录本的拼接
bowtie2 -p 8 -N 1 -x 参考基因组索引 -U xxxxxx.fastq -S xxxx.sam
#这个是bowtie2单端的，双端的就是1和2，单端才用U
#总之用到的话就去问或者查，实际问题需要实际考虑
#N是允许的错配数目，默认是0，这个数目指的不是碱基数目
#max # mismatches in seed alignment; can be 0 or 1 (0)
#总之seed是算法上面的事情，不是read
bowtie -m 1 -q /mnt/data1/genome/tair10/bowtie_index ${id}.fastq -S ${id}.sam
#bowtie -m 1 -p 8 -q /bios-analysis10/omics2021_post/21210700069/genome/rice/bowtieindex/bowtieindex WT_input.fastq -S WT_input.sam
#比对,是俺的代码，bowtie1是不能少了这个q的，双端就是1和2
bowtie2 -p 8 -x /mnt/data5/wjc/genome_scripts/genome/arb/ensemble/index/bowtie2/bowtie2_index -1 ${id}_R1_001.fastq -2 ${id}_R2_001.fastq -S ${id}.sam
#双端
#单端用-U
#RNA-seq用hisat2，chipseq用bowtie2




###################自己写的脚本
#!/bin/sh
set -e
ls -1 *_R1.fastq |awk -F "_R1.fastq" '{print $1}' | while read id
do
echo "${id}"
hisat2 -p 8 --dta --rna-strandness RF -x /mnt/data5/wjc/genome_scripts/genome/arb/ensemble/index/hisat2/hisat2_index -1 ${id}_R1.fastq -2 ${id}_R2.fastq -S ${id}.sam
echo "map done"
samtools view -q 20 -bS ${id}.sam > ${id}_q20.bam
echo "q20 bam done"
samtools sort ${id}_q20.bam -o ${id}_q20_s.bam
echo "sort done"
samtools index ${id}_q20_s_rmdup.bam
echo "index done"
samtools view -F 4 -c ${id}.sam
samtools view -F 4 -c ${id}_q20.bam
samtools view -F 4 -c ${id}_q20_s.bam
echo "reads count done"
rm ${id}_q20.bam  
rm ${id}.sam
done
mkdir Fastq
mv *.fastq Fastq
echo "all done"



featureCounts -T 4 -p -s 2 -t exon -g gene_id -F GTF -a /mnt/data1/genome/tair10/new_TAIR10_GFF3_genes_transposons.gtf -o ${id}.count ${id}.bam
#p是双端测序的意思，单端测序就不要加了，s是链特异性的意思，2是dutp（比对到reverse），-g取得是gene名字，输出的counts就会弄出来，默认就是这个，可以不加，
#-F是指定GTF，gff也不能用，反正ens的不能用，T线程，默认都是1，t是带有外显子的行才能被提取，默认就是外显子
sed '1d' ${id}.count | awk '{print $1"\t"$7}' > ${id}.rawcount
#1d是删除第一行的意思

awk '{if($3=="gene") print $0}' Arabidopsis_thaliana.TAIR10.52.gff3 | awk -F ":" '{print $1"\t"$2}' | awk -F ";" '{print $1}' |awk '{ print $1"\t"$4"\t"$5"\t"$10"\t"$6"\t"$7}' > Arabidopsis_thaliana.TAIR10.52.col6.bed
#制作6列bed

conda activate python3 
#检验数据相关性
/mnt/data5/XWH/LCZ_RNAseq/Correlation
vim multiBamSummary.sh
#!/bin/sh
set -e
multiBamSummary bins \
-bs 1000 \
--bamfiles /mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Col-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Col-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Col-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Clf-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Clf-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Clf-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1Clf-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1Clf-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_16Hos1Clf-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Col-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Col-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Col-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Clf-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Clf-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Clf-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1-3_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1Clf-1_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1Clf-2_211105S_combined__q20_s.bam \
/mnt/data5/XWH/LCZ_RNAseq/Bam/trim1_22Hos1Clf-3_211105S_combined__q20_s.bam \
--labels 16col_1 16col_2 16col_3 16Clf-1 16Clf-2 16Clf-3 16Hos1-1 16Hos1-2 16Hos1-3 16Hos1Clf-1 16Hos1Clf-2 16Hos1Clf-3 \
22col_1 22col_2 22col_3 22Clf-1 22Clf-2 22Clf-3 22Hos1-1 22Hos1-2 22Hos1-3 22Hos1Clf-1 22Hos1Clf-2 22Hos1Clf-3 \
-p 8 \
-out readCounts1.npz \
--outRawCounts readCounts1.tab
echo "done"
nohup sh multiBamSummary.sh > multiBamSummarynohup 2>&1 &
#爆慢，等一两个小时吧
conda activate python3
vim plotCorrelation1.sh
#!/bin/sh
set -e
plotCorrelation \
-in readCounts1.npz \
--corMethod pearson \
--skipZeros \
--plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot scatterplot \
--removeOutliers \
-o scatterplot_PearsonCorr1.png   \
--outFileCorMatrix PearsonCorr1.tab
echo "1 done"
plotCorrelation \
-in readCounts1.npz \
--corMethod pearson \
--skipZeros \
--plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap coolwarm \
--plotNumbers \
--removeOutliers \
-o heatmap_PearsonCorr2.png   \
--outFileCorMatrix PearsonCorr2.tab
echo "2 done"
plotCorrelation \
-in readCounts1.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot scatterplot \
--removeOutliers \
-o scatterplot_SpearmanCorr1.png   \
--outFileCorMatrix SpearmanCorr1.tab
echo "3 done"
plotCorrelation \
-in readCounts1.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap coolwarm \
--plotNumbers \
--removeOutliers \
-o heatmap_SpearmanCorr2.png   \
--outFileCorMatrix SpearmanCorr2.tab
#一、差异基因list
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
suppressMessages(library(DESeq2))
#加载这个包有点慢，要等一小会儿
sample1  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16clf_1.rawcount")
sample2  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16clf_2.rawcount")
sample3  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16clf_3.rawcount")
sample4  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16col_1.rawcount")
sample5  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16col_2.rawcount")
sample6  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16col_3.rawcount")
sample7  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1_1.rawcount")
sample8  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1_2.rawcount")
sample9  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1_3.rawcount")
sample10  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1clf_1.rawcount")
sample11  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1clf_2.rawcount")
sample12  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\16hos1clf_3.rawcount")
sample13  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22clf_1.rawcount")
sample14  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22clf_2.rawcount")
sample15  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22clf_3.rawcount")
sample16  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22col_1.rawcount")
sample17  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22col_2.rawcount")
sample18  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22col_3.rawcount")
sample19  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1_1.rawcount")
sample20  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1_2.rawcount")
sample21  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1_3.rawcount")
sample22  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1clf_1.rawcount")
sample23  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1clf_2.rawcount")
sample24  <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\downstream\\rawcount\\22hos1clf_3.rawcount")
rawcount <- data.frame(row.names = sample1[,1], 
                       clf_1_16 = sample1[,2], clf_2_16 = sample2[,2], clf_3_16 = sample3[,2],
                       col_1_16 = sample4[,2], col_2_16 = sample5[,2], col_3_16 = sample6[,2],
                       hos1_1_16 = sample7[,2], hos1_2_16 = sample8[,2], hos1_3_16 = sample9[,2],
                       hos1clf_1_16 = sample10[,2], hos1clf_2_16 = sample11[,2], hos1clf_3_16 = sample12[,2],
                       clf_1_22 = sample13[,2], clf_2_22 = sample14[,2], clf_3_22 = sample15[,2],
                       col_1_22 = sample16[,2], col_2_22 = sample17[,2], col_3_22 = sample18[,2],
                       hos1_1_22 = sample19[,2], hos1_2_22 = sample20[,2], hos1_3_22 = sample21[,2],
                       hos1clf_1_22 = sample22[,2], hos1clf_2_22 = sample23[,2], hos1clf_3_22 = sample24[,2])

#sample1的第一列作为列名，每一个文件的第二列作为一列，构建一个数据框
#列名不能用数字开头，也不要有-
condition <- factor(rep(c("16clf", "16col", "16hos1", "16hos1clf","22clf", "22col", "22hos1", "22hos1clf"), each = 3))
#顺序是有所谓的，必需要和上面好好对应起来，不过不一定要一样，但是名字是在这里取得，必需要和下面一样
#构建矩阵在这里构建就好
dds <- DESeqDataSetFromMatrix(rawcount, DataFrame(condition), design = ~condition)
dds_filter <- dds[rowSums(counts(dds)) > 1,]
dds_DES <- DESeq(dds_filter)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#弹出5行,有点慢
#DESeq2有着自己的标准化，不是根据基因长度来的，什么拟合、离散度各种高级算法
#我也见过师兄根据基因组注释文件进行标准化RPKM，都差不多，目的是为了得到想要的结果
########################################################下面是以16℃野生型比22℃野生型为例子
res1 <- results(dds_DES, contrast = c("condition", "16col", "22col"))
#这个地方决定了是谁比谁，这是分组矩阵，谁在前就是谁比谁，这是包里带的函数
diff1 <- subset(res1, abs(res1$log2FoldChange) >= log2(1.5) & res1$pvalue <= 0.05)
#subset是筛选函数，abs是计算绝对值函数
up1 <- subset(res1, res1$log2FoldChange >= log2(1.5) & res1$pvalue <= 0.05)
down1 <- subset(res1, res1$log2FoldChange <= -log2(1.5) & res1$pvalue <= 0.05)
write.table(res1, file = "16col_vs_22col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#上面的是全部
write.table(diff1, file = "16col_diff_22col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#差异基因
write.table(up1, file = "16col_up_22col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#上调基因
write.table(down1, file = "16col_down_22col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#下调基因
nrow(up1); nrow(down1)
#nrow是统计行数，看一下上下调基因各有多少个
#[1] 2360
#[1] 1627
########################################################下面是以16℃-clf比22℃-clf为例子
res1 <- results(dds_DES, contrast = c("condition", "16clf", "22clf"))
diff1 <- subset(res1, abs(res1$log2FoldChange) >= log2(1.5) & res1$pvalue <= 0.05)
up1 <- subset(res1, res1$log2FoldChange >= log2(1.5) & res1$pvalue <= 0.05)
down1 <- subset(res1, res1$log2FoldChange <= -log2(1.5) & res1$pvalue <= 0.05)
write.table(res1, file = "16clf_vs_22clf_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
write.table(diff1, file = "16clf_diff_22clf_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
write.table(up1, file = "16clf_up_22clf_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
write.table(down1, file = "16clf_down_22clf_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
nrow(up1); nrow(down1)


#二、火山图
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
gene <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16col_vs_22col_p0.05_FC(1.5).xls")
#之所以用的是vs而不是diff是因为画火山图必须把非差异基因都画进去，不能仅仅画差异基因
gene_up <- intersect(which(gene$log2FoldChange >= log2(1.5)),which(gene$pvalue <= 0.05))
#intersect是查找两个数据组交集的函数，which是筛选函数
gene_down <- intersect(which(gene$log2FoldChange <= -log2(1.5)), which(gene$pvalue <= 0.05))
sig <- rep("no", times = nrow(gene))
#重复no，次数是gene变量的行数，实际上就是每个基因重复一个no
sig[gene_up] <- "up"
#在sig数据集中，把上调的基因从no的换成up
sig[gene_down] <- "down"  
#在sig数据集中，把下调的基因从no换成down
sig <- factor(sig,levels=c("up","down","no")) 
#将其变换成因子，重新赋值，实际上就是把引号去掉了，方便与进行后续操作
gene$type <- sig
#把sig重新赋值给gene的type，相当于gene加了一列，这个地方可以区分图例的类型
gene$log10FDR <- -log(gene$pvalue, 10)
#所以一般情况下:我们可以认为Q value = FDR = adjusted p value，这里就用了p作为fdr
#qvalue:衡量错误发现率的指标（False discovery rate，简称FDR，所有检验中假阳性的概率）
#将p值取log以10为底数的负对数，目的是让他们画图时候看起来更好看，边界不会跑的到处都是
#这一列就是加了新的一列，是gene$log10FDR
gene$log10FDR2 <- gene$log10FDR
#再加一列gene$log10FDR，这不是冗余代码，可以控制边界（就是下面的，大于多少的取多少，小于多少的取多少），而且原数据也不会被覆盖
gene$log10FDR2[gene$log10FDR2 > 20] <- 20
#设定边界，p值取10为底的负对数大于20的取20就好
gene$log2FoldChange2 <- gene$log2FoldChange
#和上面的一样，不是冗余代码，用于控制边界
gene$log2FoldChange2[gene$log2FoldChange2 > 4] <- 4
gene$log2FoldChange2[gene$log2FoldChange2 < -4] <- -4
#设定好边界，大于4的取4，小于-4的取-4
pdf("volcanoplot_16col_vs_22col_p0.05_FC(1.5).pdf",width = 7,height = 6)
#一会要保存的图就是这个名字，宽度高度
xline <- c(-log2(1.5),log2(1.5))
yline <- -log(0.05, 10)
#一会儿画图时候用到的俩变量，是画虚线用到的，这个设置的值是在这里就设置好了的，如果没有这个的话自然也就不会有最后的dev.off()，是在R里画
library(ggplot2)
qplot(y = gene$log10FDR2, x = gene$log2FoldChange2, ylab = "-log10(pvalue)", xlab = "log2(fold change)", size=I(1), colour = sig, shape = sig) + 
  ylim(0, 20) + xlim(-4, 4) + scale_color_manual(values = c("up" = "red", "down" = "green", "no" = "black")) + 
  scale_shape_manual(values = c("up" = 19, "down" = 19, "no" = 19)) + 
  geom_hline(yintercept = yline, lty = 2, size = I(0.5), colour = "grey") + geom_vline(xintercept = xline, lty = 2,size = I(0.5), colour = "grey") +
  theme_bw()  + theme(panel.background = element_rect(colour = "black", size = 1, fill = "white"), panel.grid = element_blank())
dev.off()
#y轴画什么，x轴画什么，y轴叫什么，x轴叫什么，点的大小是1，点的颜色和形状是根据sig变量来区分的（即no的怎么样，up的怎么样子，down的怎么样）
#y轴和x轴的边界，scale_color_manual是点的颜色
#scale_shape_manual是点的形状，19是圆黑圈
#geom_hline和geom_vline是两条灰色的平行虚线，具体和谁平行之前都界定好的，样式是2类，粗细0.5，颜色灰色
#theme_bw()去掉背景的灰色，theme(panel.background = element_rect是子边框和各类边框的合并等，总之无需变动，panel.grid = element_blank()不要网格
#dev.off()是保存pdf到工作路径


#三、差异基因相关性
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
gene1 <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16clf_vs_22clf_p0.05_FC(1.5).xls")
gene2 <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16col_vs_22col_p0.05_FC(1.5).xls")
#读取数据，自带row.names，注意谁在上面谁就是横轴的
#之所以用的是vs而不是diff是因为构建all数据框时候行数需要一致，不然不好弄格式，就会报错
#目的是比较两组差异基因的相关性，两组差异基因实际上就是四组数据
#之所以用的是vs不是diff是因为下面的all里，有的基因在这里是差异基因，在另一个数据集里就不是差异基因了，不能让点点只有横纵坐标中的一个
gene1 <- gene1 [order(row.names(gene1),decreasing=F),]
gene2 <- gene2 [order(row.names(gene2),decreasing=F),]
#order是排序，按照的是第一列的大小升序排列(降序=F就是升序)，row.names是所有基因的名字，也即行名，这个是在之前做差异基因时候写进去的row.names=T
#然后是升序排列，这个排序方式非常烧脑，具体细节去b站看看就好，中括号是用来抓取元素的，最终的结果就是根据基因名字的大小从小到大排序
#实际上order的意思，第一个的意思是最小值在我原来数据中的位次，剩下的也一样，都是这样的，然后是次小值在原来数据中的位置
#最终的结果是基因从小到大的排序，这么做的目的是为了让结果更加整齐
all <- data.frame(gene1 [], gene2 [])
#将数据汇总到1个数据框，因为gene2的列名和gene1是一样的，所以新加的列会多一个.1,
#实际上就是在gene2的都变成了.1，也即加了很多列，不是行和列，就是列，竖的多了一倍
all_gene1_up   <- subset(all, all$log2FoldChange >= log2(1.5) & all$pvalue <= 0.05 )
#subset是筛选函数，gene1里上升的，gene2里不关心，且有p值的因素，最终黑色的点还是原来p值不显著的
all_gene1_down <- subset(all, all$log2FoldChange <= -log2(1.5) & all$pvalue <= 0.05 )
#gene1里下降的，gene2里不关心
all_gene2_up   <- subset(all, all$log2FoldChange.1 >= log2(1.5) & all$pvalue.1 <= 0.05 )
#gene2里上升的，gene1里不关心
all_gene2_down <- subset(all, all$log2FoldChange.1 <= -log2(1.5) & all$pvalue.1 <= 0.05 )
#gene2里下降的，gene1里不关心
print(c(nrow(all_gene1_up),nrow(all_gene1_down),nrow(all_gene2_up),nrow(all_gene2_down)))
#读取行数，也就是gene1上调的，gene1下调的，gene2上调的，gene2下调的
allbind <- rbind(all_gene1_up,all_gene1_down,all_gene2_up,all_gene2_down)
#把这些东西都合在一个矩阵里，行数翻了四倍
#rbind根据行进行合并，就是行的叠加，m行的矩阵与n行的矩阵rbind()最后变成m+n行
all <- unique(allbind)
#unique函数把重复的去除掉，这样就就能到新all，前面的all可以统计每个基因单独上下调多少
#这个时候的all已经混成一锅粥了，有gene1上调的，gene1下调的，gene2上调的，gene2下调的，及其各种交集，这里的all是四组差异基因
#之后提取的红色和绿色的点是同时上调/下调的基因
#这里的all比原来的all不一样的，因为原来的all里啥都有，现在的all里只有差异基因了
#但是为什么不一开始就用差异基因呢，就是因为这样做起来能保证行数都是相同的，不会报错
gene_up <- subset(all, all$log2FoldChange >= log2(1.5) & all$pvalue <= 0.05 & all$log2FoldChange.1 >= log2(1.5) & all$pvalue.1 <= 0.05 )
#同时满足很多条件，即在gene1和gene1里都上调
gene_down <- subset(all, all$log2FoldChange <= -log2(1.5) & all$pvalue <= 0.05 & all$log2FoldChange.1 <= -log2(1.5) & all$pvalue.1 <= 0.05 )
#同时满足很多条件，即在gene1和gene1里都下调
print(c(nrow(gene_up), nrow(gene_down)))
#输出同时上调的基因和同时下调的基因
v1=which(all$log2FoldChange >= log2(1.5))
v2=which(all$pvalue <= 0.05)
v3=which(all$log2FoldChange.1 >= log2(1.5))
v4= which(all$pvalue.1 <= 0.05)
#which函数的功能是已知value返回坐标，就是几有几没有，这时候上面按照顺序排列的优势就体现出来了，非常直观
l1 <- list(v1, v2, v3, v4)
#把这四个列出来，有四个大矩阵，不是取交集
gene_up_location <- Reduce(intersect,l1)
#intersect是为了计算变量的交集
#Reduce的意思是对l1中的矩阵做运算，v1和v2取交集，再和v3取，最终全部取完交集，不然的话要写很多代码
#最终得到了都上调的基因的位置（第几个基因上调了）
v1=which(all$log2FoldChange <= -log2(1.5))
v2=which(all$pvalue <= 0.05)
v3=which(all$log2FoldChange.1 <= -log2(1.5))
v4= which(all$pvalue.1 <= 0.05)
l2 <- list(v1, v2, v3, v4)
gene_down_location <- Reduce(intersect,l2)
#这是下调的，和上面的类似
sig <- rep("no", times = nrow(all))
sig[gene_up_location] <- "up"
sig[gene_down_location] <- "down"  
sig <- factor(sig,levels=c("up","down","no")) 
#上面想知道位置本质上还是想重复上面画第一个热图的操作，没什么别的意思，目的就是为了把有些地方换成down和up
#这样这些都上调的overlap就是红色的，而都下调的overlap就是绿色的，就可以画热图了，这个热图依据的就不再是一组p和foldchange，而是两组，比较复杂，可以看到重叠的地方
#这时候展示在图上的点都是差异基因，但是会进一步筛选
all$type <- sig
#在数据加一列，标记是否有显著差异，这个差异是指的是两组是否都上调，是否都下调
xline <- c(-log2(1.5),log2(1.5))
yline <-  c(-log2(1.5),log2(1.5))
#画图时候用到的俩变量，虚线用到的，值设置的值是根据最开始的
cor <- cor.test(all$log2FoldChange, all$log2FoldChange.1, alternative = "two.sided", method = "pearson")
print(cor)
#检验数据相关性
#这个是R自带的，不是包的东西，双侧检验，皮尔森系数，检测相关性，最终会得出一个值，这个值就是相关性系数cor
all$log2FoldChange[all$log2FoldChange > 10] <- 10
all$log2FoldChange[all$log2FoldChange < -10] <- -10
all$log2FoldChange.1[all$log2FoldChange.1 > 10] <- 10
all$log2FoldChange.1[all$log2FoldChange.1 < -10] <- -10
#将log2FoldChange,log2FoldChange.1太大的值赋值到较小的值，这样画图的时候不会有超出范围的（和前面的火山图没区别），火山图可以横着画也可以竖着画
#下面是画Pdf的结果
library(ggplot2)
library(showtext)
showtext_auto(enable = TRUE)
font_add('TNR', 'times.ttf')
p <- ggplot(all, aes(x = log2FoldChange,y = log2FoldChange.1, colour = type)) +
#aes是映射到图，人家就是这样用的，没什么好说的，ggplot的固定用法
#这种加号是ggplot的基本用法罢了，type是前面的（all里的一列，也就是sig），也是图例的文字
#x轴是gene1的变化倍数（16clf/22clf），y是gene2的变化倍数(16col/22col)，横坐标和纵坐标都是Log2FC
#这时候p值已经不能被体现在图中了，是体现在点的颜色里的，当然点的颜色同时也体现了倍数变化
geom_point() +
#这个就是画散点图，括号里没东西的意思就是不做任何调整，就是黑色的小圆点罢了
theme(legend.position='top') +
#控制图例在图中的位置利用theme(legend.position）参数 该参数对应的设置如下"none", "left", "right", "bottom", "top"
theme(panel.grid=element_blank(), panel.background=element_rect(color="black",size = 1, fill = "white")) +
#不要网格,边框上的事情和画普通火山图的qplot是一样的
scale_discrete_manual(values = c("red","green","black"), aesthetics = 'colour') +
#这是一个固定格式，指的是点的颜色，就是这样用的没什么好说的，value指的是点的颜色，color是根据type来决定的
geom_hline(yintercept = yline, lty = 2, size = I(0.5), colour = "grey") + 
#实际上这个玩意和差异的界定没关系，只不过在qplot火山图里头正好设置的边界是一样的
geom_vline(xintercept = xline, lty = 2,size = I(0.5), colour = "grey") +
#这俩和qplot那个是一样的，是虚线，在上面设置好了
labs(x = "Log2FC  16clf/22clf",y = "Log2FC  16col/22col") +
#自定义x轴和y轴的名字，必须加引号，如果忘记加引号，实际上加上这个就是一个闭包，没法单独跑通了
theme(axis.title.x= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
#family的意思是字体族，vjust是0-1之间，决定的x轴下面的数字的位置
theme(axis.text.x =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) +
#决定的是文本，也就是x轴上的数字
theme(axis.title.y= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
#决定的是题目，也就是y轴旁边的数字
theme(axis.text.y =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) + 
#决定的是文本，也就是y轴上的数字
#这些东西决定的是x轴上面和下面字体的颜色、大小和距离上下边框的相对位置
scale_x_continuous(limits = c( -10,10)) +
scale_y_continuous(limits = c( -10,10))
#这俩是修改坐标轴的显示范围
#p or print(p)就是在R中看
#或者运行ggave存下来在你的路径里找
ggsave(p, file="volcanoR_clf16vsclf22_pk_col16vscol22.pdf",width = 5,height = 5) 

###############################################################################################下面的是程序化模式
volcano_correlation <- function(Path,file1,file2,log2FC,p_value,positive_value,negative_value,pdfname,xlab,ylab,xrange,yrange) {
  setwd(Path)
  #这是一个function函数，在后面会设置的，会自动识别，，相当于R的脚本，别弄错输入文档的引号
  library(ggplot2)
  library(showtext)
  showtext_auto(enable = TRUE)
  font_add('TNR', 'times.ttf')
  #加载特定字体，如果没有的话就会报错误，是两种字体，如果报错的话就是下面两行:
  #Error in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  : 字体类别出错
  #In addition: There were 50 or more warnings (use warnings() to see the first 50)
  #如果要是没有这个的话下面就不要有那个family="TNR"，这其实只是一个简单加载字体的包
  gene1 <- read.delim(file1)
  gene2 <- read.delim(file2)
  #加载数据
  gene1 <- gene1 [order(row.names(gene1),decreasing=F),]
  gene2 <- gene2 [order(row.names(gene2),decreasing=F),]
  all <- data.frame(gene1 [], gene2 [])
  all_gene1_up   <- subset(all, all$log2FoldChange >= log2FC & all$pvalue <= p_value )
  #subset是筛选函数
  all_gene1_down <- subset(all, all$log2FoldChange <= -log2FC & all$pvalue <= p_value )
  all_gene2_up   <- subset(all, all$log2FoldChange.1 >= log2FC & all$pvalue.1 <= p_value )
  all_gene2_down <- subset(all, all$log2FoldChange.1 <= -log2FC & all$pvalue.1 <= p_value )
  #选出变化显著的基因，这个log2FC可以在后面设置，p值也是
  print(c(nrow(all_gene1_up),nrow(all_gene1_down),nrow(all_gene2_up),nrow(all_gene2_down)))
  #输出一组向量，第一组比较里上升的差异基因，第一组比较里下调的差异基因，第二组比较里上升的差异基因，第二组比较里下调的差异基因
  #最终在屏幕上出来的结果就是个简单的统计数
  #all_gene1_up,all_gene1_down,all_gene2_up,all_gene2_down, gene_up,gene_down,cor（这个是最终算出来的相关性，在后面有的）
  allbind <- rbind(all_gene1_up,all_gene1_down,all_gene2_up,all_gene2_down)
  all <- unique(allbind)
  gene_up <- subset(all, all$log2FoldChange >= log2FC & all$pvalue <= p_value & all$log2FoldChange.1 >= log2FC & all$pvalue.1 <= p_value )
  gene_down <- subset(all, all$log2FoldChange <= -log2FC & all$pvalue <= p_value & all$log2FoldChange.1 <= -log2FC & all$pvalue.1 <= p_value )
  print(c(nrow(gene_up), nrow(gene_down)))
  ##取出显著上调、下调的overlap的gene
  gene_up_location <- Reduce(intersect,list(v1=which(all$log2FoldChange >= log2FC), 
  v2=which(all$pvalue <= p_value), v3=which(all$log2FoldChange.1 >= log2FC), v4= which(all$pvalue.1 <= p_value) )) 
  #which函数的功能是已知value返回坐标,reduce的意思是前俩运算，再和第三个算，是为了节省中间变量，intersect是为了计算变量的交集,list就是列出来而已
  #最终得到了都上调的基因和都下调的基因的位置
  gene_down_location <- Reduce(intersect,list(v1=which(all$log2FoldChange <= -log2FC), 
  v2=which(all$pvalue <= p_value), v3=which(all$log2FoldChange.1 <= -log2FC), v4= which(all$pvalue.1 <= p_value) )) 
  #查看上下调基因的位置
  sig <- rep("no", times = nrow(all))
  sig[gene_up_location] <- "up"
  sig[gene_down_location] <- "down"  
  sig <- factor(sig,levels=c("up","down","no")) 
  all$type <- sig
  xline <- c(-log2FC,log2FC)
  yline <- c(-log2FC,log2FC)
  cor <- cor.test(all$log2FoldChange, all$log2FoldChange.1, alternative = "two.sided", method = "pearson")
  print(cor)
  all$log2FoldChange[all$log2FoldChange > positive_value] <- positive_value
  all$log2FoldChange[all$log2FoldChange < negative_value] <- negative_value
  all$log2FoldChange.1[all$log2FoldChange.1 > positive_value] <- positive_value
  all$log2FoldChange.1[all$log2FoldChange.1 < negative_value] <- negative_value
  p <- ggplot(all, aes(x = log2FoldChange,y = log2FoldChange.1, colour = type)) +
    geom_point() +
    theme(legend.position='top') +
    theme(panel.grid=element_blank(), panel.background=element_rect(color="black",size = 1, fill = "white")) +
    scale_discrete_manual(values = c("red","green","black"), aesthetics = 'colour') +
    geom_hline(yintercept = yline, lty = 2, size = I(0.5), colour = "grey") + 
    geom_vline(xintercept = xline, lty = 2,size = I(0.5), colour = "grey") +
    labs(x = xlab,y = ylab) +
    theme(axis.title.x= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
    theme(axis.text.x =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) +
    theme(axis.title.y= element_text(size=20, family="TNR",color="black", vjust=0.5, hjust=0.5)) +
    theme(axis.text.y =  element_text(size=18, family="TNR",color="black",vjust=0.5, hjust=0.5)) + 
    scale_x_continuous(limits = c( -xrange,xrange)) +
    scale_y_continuous(limits = c( -yrange,yrange)) 
  ggsave(p, file=pdfname,width = 5,height = 5) 
}
#下面是赋予各种值
volcano_correlation("C:\\Users\\94526\\Desktop\\R工作路径",
"C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16clf_vs_22clf_p0.05_FC(1.5).xls", 
"C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16col_vs_22col_p0.05_FC(1.5).xls",
log2(1.5),0.05,10,-10,"volcanoR_clf16vsclf22_pk_col16vscol22.pdf","Log2FC  16clf/22clf","Log2FC  16col/22col",10,10)
#这个就是设定的参数，和上面的function函数接合，和单独分步骤画的没有任何区别


#四、热图
#以下是16clf-16col上调的基因为例子画热图
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
genelist <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16clf_up_16col_p0.05_FC(1.5).xls")
all <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\22度和16度野生型相比突变体改变量.xls")
#如果变量名包含首字母为数字、#、$等情况时，则会自动加上X.，使变量看上去更像一个字符型变量。
#这里的这个文档是把所有的东西都放在了一起，用的是excel手动粘贴
#除了第一列是基因名，剩下18列，3*6.3分别是Log2FC，6分别是16clf/16wt,16hos1/16wt,16hos1clf/16wt,22clf/22wt,22hos1/22wt,22hos1clf/22wt,wt就是col
#目的是为了画热图时候方便，而且提供了完整的基因清单,用下面的命令来制造排好序的excel
#setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
#gene1 <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16clf_vs_22clf_p0.05_FC(1.5).xls")
#gene2 <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16col_vs_22col_p0.05_FC(1.5).xls")
#gene1 <- gene1 [order(row.names(gene1),decreasing=F),]
#gene2 <- gene2 [order(row.names(gene2),decreasing=F),]
#write.table(gene1, file = "s_16clf_vs_22clf_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#必需要排好序，然后才能粘贴啊！
row.names(all) <- all[,1]
#第一列取出来，然后赋值给行名，目的是把行名给替换成第一列的基因名
#因为是手动粘贴的，所以没有行名，所以要整一个行名
#用R做的话就不用，比较方便
diyi <- row.names(genelist)
dier <- row.names(all)
#把行名拿出来构成矩阵
all1 <- subset(all, is.element(dier,diyi) == T)
#is.element(dier,diyi) == T
#用于检查第一个对象的元素是否存在于第二个对象中,T的意思是在的话就在那个位置写上TRUE，但是本质上还是genelist矩阵，不是只有TRUE和FALSE
#不信你输出一下excel试试
#wwwwjjcc<-is.element(dier,diyi) == T
#write.table(wwwwjjcc, file = "ss.xls", sep = "\t", row.names = T)
#下面是is.element的一个例子，基本上看大的在小的里面有没有，反过来的话就全有了（真子集）
#eg:jcc<-c(1,2,3)
#ity<-c(2,8)
#is.element(jcc,ity) == T
#[1] FALSE  TRUE FALSE
all1 <- all1[,c("X16Clf_vs_16Col_log2FoldChange","X16Hos1_vs_16Col_log2FoldChange","X16Hos1Clf_vs_16Col_log2FoldChange")]
#取出三列有用的信息，都是FC-change
colnames(all1)[which(names(all1) == "X16Clf_vs_16Col_log2FoldChange")] <- "16_clf"
colnames(all1)[which(names(all1) == "X16Hos1_vs_16Col_log2FoldChange")] <- "16_hos1"
colnames(all1)[which(names(all1) == "X16Hos1Clf_vs_16Col_log2FoldChange")] <- "16_hos1clf"
#colnames基本用法是设置行或列标题的名称,name是展示所有元素的列名，which是返回位置的意思，例如colnames(all)[2]是第二个元素的名字
#取出genelist的16度的Clf,Hos1,Hos1Clf并改掉列名
#这么做的目的是，不仅包括了genelist自身的信息，还包含了别的组的信息，比如有的基因在16clf里上调但是在别人里下调，因此不是无用功啊
suppressMessages(library(ComplexHeatmap)) 
#和下面两个分别提供Heatmap函数、提供colorRamp2函数、提供brewer.pal函数
suppressMessages(library(circlize)) 
suppressMessages(library(RColorBrewer))
#####################################################################简单进行热图的绘制，也没聚类，每次画出来的结果都是不一样的
#所有每次存的表格也不同,因为聚类每次都是随机的，虽然结果也差不多，但是从1-10那确实是完全随机的（eg:第一次a、b基因在k1，第二次在k5）
#最终聚类就是把kmeans分了十类，但是很明显这不一定是真的，因为图都乱七八糟的，而且聚类也不怎么准确
heatdata <- all1
k=kmeans(heatdata,10)
#将heatdata数据聚类,这是分成了十份
heatdata$k=k$cluster
#将聚类的值写到heatdata里面，k$cluster是分类结果,就是给heatdata加了一列
htframe_sort=heatdata[order(heatdata$k),]
#根据K值对heatdatas数据重新排序
htmatrix_sort <- as.matrix(htframe_sort[,1:3])
#将数据框转化成矩阵，并取得前三列，因为之后heatmap()函数需要输入的是矩阵
Heatmap(htmatrix_sort ,col=colorRamp2(c(3, log2(1.5), 0, -log2(1.5), -3), brewer.pal(n=5, name='RdBu')),show_row_names=FALSE,
        cluster_rows=FALSE,cluster_columns=FALSE, width = unit(3, "cm"),name = "log2FC",column_title = "热图")
#-3到3的区间被线性插入值用于获取对应的颜色，值大于0.4515的就是红色，值小于-0.4515的就是蓝色
#值大于3的被映射为深红色，小于-3的被映射为深蓝色，这个说的不是图例，是图的颜色，log2(1.5)=0.4515
#n分了几份，前面就要写几个区间值，RdBu是红-蓝渐变
#heatmap绘图,show_row_names是否显示行名称。默认值为TRUE,现在不要行名，不然太多了
#show_column_names是否显示列名称。默认值为TRUE，这里是要有的
#cluster_rows=FALSE,如果为TRUE，则在行上创建聚类的簇；
#cluster_columns=FALSE,如果为TRUE，则在列上创建聚类的簇。
#width是图的宽度，cm是单位，也可以是mm
#name是图例的名字
#column_title是大题目
#输完这一步就能在R里画出热图
write.table(htframe_sort,"htframe_sort_used1.xls",sep="\t")
#写表格

###############################################################下面开始调整参数，上面的就是画个图，大致看一下，下面可以矫正k，重新分类或者不分类
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
suppressMessages(library(ComplexHeatmap)) 
suppressMessages(library(circlize)) 
suppressMessages(library(RColorBrewer))
htframe_sort <- read.delim("htframe_sort_used1.xls")
htframe_sort2 <- rbind( subset(htframe_sort,htframe_sort[,4]==6),
                       subset(htframe_sort,htframe_sort[,4]==7),
                       subset(htframe_sort,htframe_sort[,4]==9),
                       subset(htframe_sort,htframe_sort[,4]==4),
                       subset(htframe_sort,htframe_sort[,4]==5),
                       subset(htframe_sort,htframe_sort[,4]==1),
                       subset(htframe_sort,htframe_sort[,4]==10),
                       subset(htframe_sort,htframe_sort[,4]==8),
                       subset(htframe_sort,htframe_sort[,4]==2),
                       subset(htframe_sort,htframe_sort[,4]==3))
#rbind是合并的意思，就是行数翻了10倍数（相比于单个），实际上就是重新排序一下k，这时候的顺序想怎么来就怎么来，但是这会决定最终的新热图顺序
#之所以现在和后面有顺序，就是为了要在找目的基因的时候好好对应起来
#下面是得到不同k值分类的位置
k6 <- which(htframe_sort2[,4]==6)
k7 <- which(htframe_sort2[,4]==7)
k9 <- which(htframe_sort2[,4]==9)
k4 <- which(htframe_sort2[,4]==4)
k5 <- which(htframe_sort2[,4]==5)
k1 <- which(htframe_sort2[,4]==1)
k10 <- which(htframe_sort2[,4]==10)
k2 <- which(htframe_sort2[,4]==2)
k3 <- which(htframe_sort2[,4]==3)
k8 <- which(htframe_sort2[,4]==8)
#就是赋值，但是为什么是乱序呢？
#因为这是人为自己定的，顺序想怎么做就怎么做，但是为了下面方便嘛
#总之就是大致瞅一眼，画的好看点，红的尽量聚在一块儿，蓝的尽量聚在一块儿。
#实际中不一定按照这个顺序
#k几都是序列号，但是下面的type[]就是把某些gene名字拿出来，根据的是k几
type <- rep("no", times = nrow(htframe_sort2))
#tpye里有gene数量个no
#type是所有的基因名字，只不过用no来代替了
type[k6] <- 1
#type里把k6的那些换成1
type[k7] <- 2
type[k9] <- 3
type[k4] <- 4
type[k5] <- 5
type[k1] <- 6
type[k10] <- 7
type[k8] <- 8
type[k2] <- 9
type[k3] <- 9
htframe_sort2$type <- type
#在数据加一列，重新分类排序,新加的是type列，也就相当于矫正后的k值
#就是在这里可以对k重新分类，可以把两个k聚在一起
write.table(htframe_sort2,"htframe_sort_used2.xls",sep="\t")
#保存再加一列重新分类排序后的矩阵
##################################################################################上面就是为了矫正k，这里有矫正（最后俩合到一起了）
#但是有的时候会不矫正，比如最终的type可能比k要少
#res参数是ppi中的标称分辨率，如果没有指定，取为72 ppi来设置文本和行宽的大小。ppi是像素密度，表示的是每英寸所拥有的像素数，是图像分辨率的单位。
tiff("heatmap1.tiff",width=1000,height=2000,res=300)
#保存图像设置，要是没有这个就会在R里画图
htmatrix_sort2 <- as.matrix(htframe_sort2[,1:3])
h1 <- Heatmap(htmatrix_sort2 ,col=colorRamp2(c(3, log2(1.5), 0, -log2(1.5), -3), brewer.pal(n=5, name='RdBu')),show_row_names=FALSE,
cluster_rows=FALSE,cluster_columns=FALSE, width = unit(3, "cm"),name="log2FC",column_title = "热图") 
#就重新画个热图呗，上面是画热图的代码，和普通的热图都是一样的
ha1_row = rowAnnotation(df = data.frame(Type =htframe_sort2$type),
                        col = list(Type = c("1" = "#99FFCC","2"="red",'3'="#0099CC",
                                            '4' =  "#CC99FF",'5'="#FF6699",'6'="pink",
                                            '7' =  "orange",'8'="light blue",'9'="light green",'10'="lawngreen",'11'="gold",
                                            '12'="wheat",'13'="maroon1",'14'="cornflowerblue")), width = unit(0.2, "cm"))
draw(ha1_row + h1)
#不必深究，只是设置颜色的，这个Type是图例的名字
#但是画热图是要有顺序的，是根据新k值的前后顺序，注意不是大小顺序，最开始初步画热图的时候是根据大小顺序的，但是根据新的k值前后顺序就可以把想要的聚在一起
#第一次画热图的时候没必要反复试，打开excel结合大致热图在重新画新热图时候调一调就行了
#但是这有个bug，就是10会排在1后面，而不是9后面，暂时没法修改，就是最终画图例时候不要超过9个吧
dev.off()
#复制一大串，在R里的老热图，工作路径中的是新热图+图例
#顺序都是根据自己的意愿自己调整的
#第一次画热图得到的结果是随机的，因为聚类结果是随机的，但是第二次就不是随机的了

#五、GO分析
setwd(dir = "C:\\Users\\94526\\Desktop")
suppressMessages(library(clusterProfiler))
#提供enrichGO函数，还可以提供enrichKEGG函数，具体使用方法可以上网去学习
suppressMessages(library(topGO))
#提供画图函数，如果不加载的画就需要用下面的方法自己画了，用ggplot
suppressMessages(library(org.At.tair.db))
#拟南芥的GO包
#上面的包用bioconductor或者install.packages安装
genelist <- read.delim("C:\\Users\\94526\\Desktop\\RNAseq_example\\Rplot\\16clf_diff_22clf_p0.05_FC(1.5).xls")
eGOBP <- enrichGO(gene = row.names(genelist), OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1, qvalueCutoff = 1, readable = F)
#这些都是由clusterProfiler提供的
#指定gene，指定
eGOBP1 <- enrichGO(gene = row.names(genelist), OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = F)
write.table(data.frame(eGOBP),"16clf_vs_22clf_p0.05_FC(1.5)_genes GO Enrichment_BP.xls",sep="\t",row.names=F)
write.table(data.frame(eGOBP1),"try.xls",sep="\t",row.names=F)
#这里的目的是设定好边界就可以筛选的更加精确，不过两种方法没什么大的区别
#意思就是q值和p值在xx范围内的才要
barplot(eGOBP1, showCategory=20,title="Enrichment_BP")
dotplot(eGOBP1, showCategory=20, title="Enrichment_BP")
#需要搞一个叫Rgraphviz的包,下一步自动载入
plotGOgraph(eGOBP1)
#画网络图
#对于拟南芥来说，gene就是差异基因对应的向量，keyType指定基因ID的类型，默认为ENTREZID, 该参数的取值可以参考keytypes(org.Hs.eg.db)的结果，
#建议采用ENTREZID, OrgDb指定该物种对应的org包的名字，ont代表GO的3大类别，BP, CC, MF，也可以选择ALL; pAdjustMethod指定多重假设检验矫正的方法
#这里默认pAdjustMethod="BH",所以这里没有写出来,cutoff指定对应的阈值，readable=TRUE代表将基因ID转换为gene symbol

####下面是用R自己画，可以根据包里的东西调整边界之类的
setwd(dir = "C:\\Users\\94526\\Desktop\\R工作路径")
library(ggplot2)
enrich <- read.table("C:\\Users\\94526\\Desktop\\R工作路径\\16clf_vs_22clf_p0.05_FC(1.5)_genes GO Enrichment_BP.xls", header=T,sep="\t")
enrich1 <- enrich[order(enrich$pvalue),]
#按照pvalue进行排序
enrich2 <- enrich1[1:20,]
#选择前20个
count <- as.numeric(unlist(strsplit(enrich2$GeneRatio,"/2409",fixed=T)))
#这个数字是excel待着第一行的行数
enrich3 <- data.frame(enrich2[,2],count,enrich2[,5])
colnames(enrich3) <- c("terms","count","PValue")
enrich3 <- enrich3[order(enrich3$count),]
enrich3$terms <- factor(enrich3$terms,levels=enrich3$terms)
p = ggplot(enrich3,aes(count,terms))
p = p + geom_point()+theme(axis.text.x = element_text(colour="black",size=1))
p=p + geom_point(aes(size=count))
pbubble = p+ geom_point(aes(color=-1*log10(PValue),size=count))
pr = pbubble+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(-log[10](pvalue)), size="Number of genes",
             x="Number of genes",y="")
pr + theme_bw() + theme(axis.text.x =element_text(size=12), axis.text.y=element_text(size=12))


#方法一
setwd(dir = "C:\\Users\\94526\\Desktop\\GO_WangFfan")
#R中的路径用/或者\\表示，输出的文件到这里取找
suppressMessages(library(clusterProfiler))
#提供enrichGO函数，还可以提供enrichKEGG函数，具体使用方法可以上网去学习，也可以用library(clusterProfiler)来加载，只是屏幕上会蹦出一些没有用的信息罢了
suppressMessages(library(topGO))
#提供画图函数，如果不加载的画就需要用下面的方法自己画了，用ggplot
suppressMessages(library(org.At.tair.db))
#拟南芥的GO包
#上面的包用bioconductor或者install.packages安装
#install.packages('topGO')
#上面如果不适合R版本的话，使用下面的方法
#install.packages('BiocManager')
#包和只需要下载一次就好了，BiocManager
#BiocManager::install('topGO')
#在R里命名尽量不要使用中文字符或者- + =等
genelist <- read.delim("C:\\Users\\94526\\Desktop\\GO_WangFfan\\genelist.txt",header=T)
#读取基因list，格式无所谓，这里是用列名字，如果只是用一列基因list单独的话，总之是要把基因这一列读成一列变量
#header参数：默认为FALSE即数据框的列名为V1,V2...,设置为TRUE时第一行作为列名
#因为编码可能会存在问题，编码格式有所不同，所以推荐使用txt文档作为读取数据的方式
eGOBP <- enrichGO(gene = genelist$gene, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1, qvalueCutoff = 1, readable = F)
#这些都是由clusterProfiler提供的
#指定gene，上面的是全部的生物学过程，下面的是明显富集到的生物学过程
eGOBP1 <- enrichGO(gene = genelist$gene, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = F)
write.table(data.frame(eGOBP),"Wangfan GO Enrichment_BP.xls",sep="\t",row.names=F)
#这里的目的是设定好边界就可以筛选的更加精确，不过两种方法没什么大的区别
#意思就是q值和p值在xx范围内的才要，因为富集也是有富集得明显和不明显的区别的
barplot(eGOBP1, showCategory=20,title="Enrichment_BP")
#用R的包自带的函数画一个，这样不能调整边界等等，可以大致看看
dotplot(eGOBP1, showCategory=20, title="Enrichment_BP")
#需要搞一个叫Rgraphviz的包,下一步自动载入
plotGOgraph(eGOBP1)
#画网络图
#对于拟南芥来说，gene就是差异基因对应的向量，keyType指定基因ID的类型，默认为ENTREZID, 该参数的取值可以参考keytypes(org.Hs.eg.db)的结果，
#建议采用ENTREZID, OrgDb指定该物种对应的org包的名字，ont代表GO的3大类别，BP, CC, MF，也可以选择ALL; pAdjustMethod指定多重假设检验矫正的方法
#这里默认pAdjustMethod="BH",所以这里没有写出来,cutoff指定对应的阈值，readable=TRUE代表将基因ID转换为gene symbol

######################下面的步骤是利用上面的数据自己画图，而不是使用R包绘图，这样可以很好地画边界#############################################
library(ggplot2)
#install.packages('(ggplot2')
setwd(dir = "C:\\Users\\94526\\Desktop\\GO_WangFfan")
enrich <- read.table("C:\\Users\\94526\\Desktop\\GO_WangFfan\\Wangfan GO Enrichment_BP.xls", header=T,sep="\t")
#sep="\t" 表示以tab（制表符）为分隔符
enrich1 <- enrich[order(enrich$pvalue),]
#按照pvalue进行排序,而非padj
enrich2 <- enrich1[1:20,]
#选择前20个
count <- enrich2$Count
#或者可以使用count <- as.numeric(unlist(strsplit(enrich2$GeneRatio,"/2409",fixed=T)))
enrich3 <- data.frame(enrich2[,2],count,enrich2[,5])
#把富集名称和p值取出来，还有count
colnames(enrich3) <- c("terms","count","PValue")
#加行名，分别是terms，count ，p值
enrich3$terms <- factor(enrich3$terms,levels=enrich3$terms)
#变成因子，没有具体的意义
p = ggplot(enrich3,aes(count,terms))
p = p + geom_point()+theme(axis.text.x = element_text(colour="black",size=1))
p=p + geom_point(aes(size=count))
pbubble = p+ geom_point(aes(color=-1*log10(PValue),size=count))
pr = pbubble+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(-log[10](pvalue)), size="Number of genes",
             x="Number of genes",y="")
pr + theme_bw() + theme(axis.text.x =element_text(size=12), axis.text.y=element_text(size=12)) 
#画图的参数，暂时不用改，具体就是一些横纵坐标的具体信息了，颜色大小什么的
####点越红表示越显著，点越大表示富集到该通路的基因越多，也是横纵标，而这种富集分析是没有纵坐标的


#########方法二,DAVID GO（利用一个网站富集，而R包只负责画图）
#有时候R包并不能很好地富集到想要的东西，这个时候就要换一种方式了
#https://david-d.ncifcrf.gov/
#https://david.ncifcrf.gov/
#看看哪个能用，都是davidgo的网站，有时候服务器就会崩溃
#1.选择start analysis
#2.把genelist粘贴到方框里面
#3.下拉选择TAIR_ID作为Select Identifier的选项
#4.选择genelist然后submit list
#5在Analyze above gene list with one of DAVID tools里面选择第一个大选项functional annotation tools
#6.点到Gene_Ontology (3 selected)默认里面可以下载bp/cc/mf的chat或者在下面的functional annotation chart大选项里选择，ctrl+a和ctrl+s保存成txt格式，打开就能看到一个list了
library(ggplot2)
setwd(dir = "C:\\Users\\94526\\Desktop\\GO_WangFfan")
enrich <- read.table("C:\\Users\\94526\\Desktop\\GO_WangFfan\\GO.txt", header=T,sep="\t")
enrich1 <- enrich[order(enrich$PValue),]
enrich2 <- enrich1[1:20,]
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
#俩处理数据的包，不用纠结，安装就好了
enrich2 <- enrich2 %>% separate(Term, c("ID","term"), "[~]")
enrich3 <- data.frame(enrich2$term, enrich2$Count,enrich2$PValue)
#将Term分成ID和term
colnames(enrich3) <- c("terms","count","PValue")
enrich3 <- enrich3[order(enrich3$count),]
enrich3$terms <- factor(enrich3$terms,levels=enrich3$terms)
p = ggplot(enrich3,aes(count,terms))
p = p + geom_point()+theme(axis.text.x = element_text(colour="black",size=1))
p=p + geom_point(aes(size=count))
pbubble = p+ geom_point(aes(color=-1*log10(PValue),size=count))
pr = pbubble+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(-log[10](pvalue)), size="Number of genes",
             x="Number of genes",y="")
pr + theme_bw() + theme(axis.text.x =element_text(size=12), axis.text.y=element_text(size=12))
#和之前的一样







#除了特别的软件，默认所有软件都装在python3里了
vim copy.sh
#!/bin/sh
set -e
cp  /mnt/USB1/2021.11.16-11.26/LCZ-hos1-ChIPseq/* /mnt/data5/wjc/aaworkflow/LCZ-hos1-ChIPseq
echo "copy done"
nohup sh copy.sh > copynohup 2>&1 

md5sum *.gz
cat *.md5
#然后丢到韦恩图网站去看数据完整性

gunzip *.gz

conda activate python3
vim qc.sh
#!/bin/sh
set -e
mkdir QC
fastqc -o QC -f fastq *.fastq
echo "qc done"
nohup sh qc.sh > qcnohup 2>&1 

vim trim1andqc.sh
#!/bin/sh
set -e
ls -1 *R1_001.fastq |awk -F "R1_001." '{print $1}' | while read id
do
echo "${id}"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-m 20 -o trim1_${id}R1_001.fastq -p trim1_${id}R2_001.fastq \
${id}R1_001.fastq ${id}R2_001.fastq
done
echo "trim1 done"
mkdir trim1QC
fastqc -o trim1QC -f fastq trim1*.fastq
echo "trim1 qc done"
mkdir Fastq
mv [^trim]*.fastq Fastq
echo "all done"
#下面是运行代码
nohup sh trim1andqc.sh > trim1andqcnohup 2>&1 

##bowtie2进行比对 双端
vim workflow1.sh
#!/bin/sh
set -e
ls -1 *R1_001.fastq |awk -F "_R1_001." '{print $1}' | while read id
do
echo "${id}"
bowtie2 -p 8 -x /mnt/data1/genome/tair10/bowtie2_index/bowtie2_index -1 ${id}_R1_001.fastq -2 ${id}_R2_001.fastq -S ${id}.sam
echo "map done"
samtools view -q 20 -bS ${id}.sam > ${id}_q20.bam
echo "q20 bam done"
samtools sort ${id}_q20.bam -o ${id}_q20_s.bam
echo "sort done"
samtools rmdup ${id}_q20_s.bam ${id}_q20_s_rmdup.bam
echo "rmdup done"
samtools index ${id}_q20_s_rmdup.bam
echo "index done"
samtools view -F 4 -c ${id}.sam
samtools view -F 4 -c ${id}_q20.bam
samtools view -F 4 -c ${id}_q20_s.bam
samtools view -F 4 -c ${id}_q20_s_rmdup.bam
echo "reads count done"
rm ${id}_q20_s.bam
rm ${id}_q20.bam 
rm ${id}.sam  
done
mkdir Fastq
mv *.fastq Fastq
echo "all done"

nohup sh workflow1.sh > workflow1nohup 2>&1 

vim workflow2.sh
#!/bin/sh
set -e
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -p 1 -b ${id}.bam -o ${id}_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"
bamCoverage -p 1 -b ${id}.bam -o ${id}_RPGC.bigwig --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 119000000
echo "RPGC bigwig done"
bedtools bamtobed -i ${id}.bam > ${id}.bed
echo "bed done"
awk '{print "chr"$0}'  ${id}.bed > ${id}_chr.bed
echo "chr done"
done
mkdir Bam
mv *.bam Bam
mv *.bam.bai Bam
mkdir Bed 
mv *.bed Bed
mkdir RPKMbigwig
mv *RPKM.bigwig RPKMbigwig
mkdir RPGCbigwig
mv *RPGC.bigwig RPGCbigwig
echo "all done"

nohup sh workflow2.sh > workflow2nohup 2>&1 


/mnt/data5/wjc/aaworkflow/chipseqtest/matrix
vim matrix.sh
#!/bin/sh
set -e
computeMatrix scale-regions \
-S /mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Col-input_211120N_S1_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Col-H2A-Z_211120N_S25_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Clf-input_211120N_S3_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Clf-H2A-Z_211120N_S27_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Hos1-input_211120N_S2_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_22Hos1-H2A-Z_211120N_S26_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Col-input_211120N_S4_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Col-H2A-Z_211120N_S28_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Hos1-input_211120N_S5_L003_q20_s_rmdup_RPKM.bigwig \
/mnt/data5/wjc/aaworkflow/chipseqtest/RPKMbigwig/trim1_16Hos1-H2A-Z_211120N_S29_L003_q20_s_rmdup_RPKM.bigwig \
-R /mnt/data5/wjc/genome_scripts/genome/arb/ensemble/Ara_final_col6.bed \
-b 1000 -a 1000 -m 3000 --binSize 10  --sortRegions keep  --missingDataAsZero  -p 8 -out all-H2AZ_matrix.mat.gz
echo "all-H2AZ matrix done"
plotHeatmap -m all-H2AZ_matrix.mat.gz -out all_H2AZ_heatmap.png
plotProfile -m all-H2AZ_matrix.mat.gz -out all_H2AZ_profile.png
echo "plot done"
gzip -d all-H2AZ_matrix.mat.gz
sed '1d' all-H2AZ_matrix.mat > all-H2AZ_matrix.mat_nohead
echo "all done"
#没解压之前才能这么用呢,其实这个画图功能还是挺强大的，但是我们实验室还是选择以R画图为主，这个参数好多都没有发掘出来，我都放到deeptools的收藏夹里面了
nohup sh matrix.sh > matrixnohup 2>&1 &
#-S输入的是bw文件，-R：后面跟gene.bed文件，该文件可以从基因注释文件（gff3格式）转化而来
#--binSize： bin大小。默认值是10。
#-b:所选参考点的上游距离。-a:所选参考点的下游距离。默认值是0。
#sortRegions:输出文件是否应该显示排序的区域。默认情况下不对区域进行排序。如果需要输出顺序与输入区域匹配，则指定“keep”。
#out输出文件，regionBodyLength, 也就是-m，是基因body的长度
#--missingDataAsZero如果设置，丢失的数据(NAs)将被视为零。 默认情况是忽略这些情况，它们将在热图中被描述为黑色区域。 
#输出的矩阵是每个基因上面不同bin的得分情况，总之对象是每一个基因，最终的图是整个基因组上每个基因的平均分布，是整体的一个概览
#单线程的话大概需要4个小时以上，16线程就是30min左右，根据情况自己弄吧，一开始是用16线程计算，后续使用较少的收尾，最后用单线程写进去，所以16线程并不是单线程速度的16倍
#plotheatmap画图的（包含了plotProfile）

#得出的结果可以用下面的方法进行画图，如果用R来画图的话，需要把这个mat文件解压，并且删去第一行，下载到本地windows

####下面是macs2找peaks
/mnt/data5/wjc/aaworkflow/chipseqtest/macs2peaks
vim macs2peaks.sh
#!/bin/sh
set -e
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup.bam  -g 118000000 -n 16clf_H2AZ -B --nomodel -f BAMPE
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Clf-H2A-Z_211120N_S27_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Clf-input_211120N_S3_L003_q20_s_rmdup.bam  -g 118000000 -n 22clf_H2AZ -B --nomodel -f BAMPE
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Col-H2A-Z_211120N_S28_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Col-input_211120N_S4_L003_q20_s_rmdup.bam  -g 118000000 -n 16col_H2AZ -B --nomodel -f BAMPE
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Col-H2A-Z_211120N_S25_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Col-input_211120N_S1_L003_q20_s_rmdup.bam  -g 118000000 -n 22col_H2AZ -B --nomodel -f BAMPE
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Hos1-H2A-Z_211120N_S29_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Hos1-input_211120N_S5_L003_q20_s_rmdup.bam  -g 118000000 -n 16hos1_H2AZ -B --nomodel -f BAMPE
macs2 callpeak -t /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Hos1-H2A-Z_211120N_S26_L003_q20_s_rmdup.bam \
-c /mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_22Hos1-input_211120N_S2_L003_q20_s_rmdup.bam  -g 118000000 -n 22hos1_H2AZ -B --nomodel -f BAMPE
echo "macs2peaks done"
#-t实验组输入-c对照组输入-g有效基因组大小-n输出文件的前缀名 -B: 会保存更多的信息在bedGraph文件中，如fragment pileup, control lambda, -log10pvalueand-log10qvalue scores
#Rscript NAME_model.r得到双峰的pdf文件
nohup sh macs2peaks.sh > macs2peaksnohup 2>&1
#下面是自己琢磨的使用范例
macs2 callpeak -t chip.bam -c input.bam -g 118000000 -n output -B (--broad) --nomodel -f BAMPE
#双末端就不需要建模拟合什么的，不要那个所谓的双峰模型了
#生成的几个文件其实都是peak的结果，其中summit和narrow配套，gappedpeak和broad匹配，这算是两个文件
#xls表格算是一个文件，如果有的话，R脚本算是一个文件
#最后还有两个分别是ip和input的bdg文件，除了excel表都可以用igv打开，反正就是不同的方法，不同的软件，call出的peak都是不太相同的，但是大体上也是差不太多的，毕竟比对后bam一样


#针对双端数据，代码为：
/mnt/data5/wjc/aaworkflow/chipseqtest/macs2peaks/PE
vim macs2callpeakPE.sh
#!/bin/sh
set -e
IP_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup.bam
Input_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup.bam
macs2 callpeak -t $IP_bam -c $Input_bam -f BAMPE -g 1.18e8  -q 0.05 --nomodel  -B  -n 16Clf-H2A-Z_q0.05-regular-nomodel --outdir 16Clf-H2A-Z_regularnomodel
echo "16Clf-H2A-Z regular nonmodel q0.05"
macs2 callpeak -t $IP_bam -c $Input_bam -f BAMPE -g 1.18e8 --broad --broad-cutoff 0.05 -q 0.05  --nomodel  -B -n 16Clf-H2A-Z_broadcutoff0.05q0.05-broad-nomodel --outdir 16Clf-H2A-Z_broadnomodel
echo "16Clf-H2A-Z broad nonmodel broadcutoff0.05 q0.05"
echo "all done"
nohup sh macs2callpeakPE.sh > macs2callpeakPEnohup 2>&1 


#针对单端数据，代码为：
/mnt/data5/wjc/aaworkflow/chipseqtest/macs2peaks/SE
#就把上面的数据当作单端处理了，但是人家实际上是双端测序
vim macs2callpeakSE.sh
#!/bin/sh
set -e
IP1_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup.bam
Input_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup.bam
macs2 callpeak -t $IP1_bam -c $Input_bam -f BAM -g 1.18e8  -q 0.05  -B  -n 16Clf-H2A-Z_q0.05-regular-model --outdir 16Clf-H2A-Z_regularmodel
echo "16Clf-H2A-Z regular model q0.05"
macs2 callpeak -t $IP1_bam -c $Input_bam -f BAM -g 1.18e8 --broad --broad-cutoff 0.05 -q 0.05  -B  -n 16Clf-H2A-Z_broadcutoff0.05q0.05-broad-model --outdir 16Clf-H2A-Z_broadmodel
echo "16Clf-H2A-Z broad model broadcutoff0.05 q0.05"
echo "all done"
nohup sh macs2callpeakSE.sh > macs2callpeakSEnohup 2>&1 


####下面是Sicer找peaks
/mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ
conda activate python2
#必须在这个下面进行，否则就报错误,还要是python2才行
vim sicer.sh
#!/bin/sh
set -e
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Clf-H2A-Z_q20_s_rmdup_chr.bed \
16Clf-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/16clf newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Clf-H2A-Z_q20_s_rmdup_chr.bed \
16Clf-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/16clf newtair10 1 200 200 0.9 600 0.05
#
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Col-H2A-Z_q20_s_rmdup_chr.bed \
16Col-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/16col newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Col-H2A-Z_q20_s_rmdup_chr.bed \
16Col-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/16col newtair10 1 200 200 0.9 600 0.05
#
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Hos1-H2A-Z_q20_s_rmdup_chr.bed \
16Hos1-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/16hos1 newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 16Hos1-H2A-Z_q20_s_rmdup_chr.bed \
16Hos1-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/16hos1 newtair10 1 200 200 0.9 600 0.05
#
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Clf-H2A-Z_q20_s_rmdup_chr.bed \
22Clf-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/22clf newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Clf-H2A-Z_q20_s_rmdup_chr.bed \
22Clf-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/22clf newtair10 1 200 200 0.9 600 0.05
#
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Col-H2A-Z_q20_s_rmdup_chr.bed \
22Col-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/22col newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Col-H2A-Z_q20_s_rmdup_chr.bed \
22Col-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/22col newtair10 1 200 200 0.9 600 0.05
#
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Hos1-H2A-Z_q20_s_rmdup_chr.bed \
22Hos1-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrow/22hos1 newtair10 1 200 200 0.9 200 0.05
sh /mnt/data1/software/SICER_V1.1/SICER/SICER.sh /mnt/data5/wjc/aaworkflow/chipseqtest/gaimingbed 22Hos1-H2A-Z_q20_s_rmdup_chr.bed \
22Hos1-input_q20_s_rmdup_chr.bed /mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broad/22hos1 newtair10 1 200 200 0.9 600 0.05
echo "callpeak done"
#
nohup sh sicer.sh > sicernohup 2>&1 
#######
["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] \
["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]
#物种就是Osmsu?或者newtair10，["redundancy threshold"]都是1，["window size (bp)"]是分辨率就是200
#剩下的也都不用改变，其实宽峰和窄峰就是却决于[gap size"]的值是200还是600
#要提前建立好文件夹


##提取差异变化显著的,只有sicer才用这个，因为macs2直接就制定了fc和q啊p什么的
##由于后面注释时染色体名字格式为1，而不是chr1，因此将chr删除
/mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/narrowpeak
#这个里面放了所有的summary文件
vim narrowpaeak.sh
#!/bin/sh
set -e
awk '{if($7>=1.5) print $0}' 16Clf-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Clf_H2AZ_FC1.5_FDR0.001.bed
#提取差异变化显著的，实际上绝大部分都是合格的
awk -F "chr" '{print $2}' sig_16Clf_H2AZ_FC1.5_FDR0.001.bed > sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
#删掉chr
wc -l sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
#统计行数，就是peaks数目
echo "16clfnarrowpeak done"
awk '{if($7>=1.5) print $0}' 16Col-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Col_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_16Col_H2AZ_FC1.5_FDR0.001.bed > sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
echo "16colnarrowpeak done"
awk '{if($7>=1.5) print $0}' 16Hos1-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Hos1_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_16Hos1_H2AZ_FC1.5_FDR0.001.bed > sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
echo "16hos1narrowpeak done"
awk '{if($7>=1.5) print $0}' 22Clf-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Clf_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Clf_H2AZ_FC1.5_FDR0.001.bed > sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
echo "22clfnarrowpeak done"
awk '{if($7>=1.5) print $0}' 22Col-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Col_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Col_H2AZ_FC1.5_FDR0.001.bed > sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
echo "22colnarrowpeak done"
awk '{if($7>=1.5) print $0}' 22Hos1-H2A-Z_q20_s_rmdup_chr-W200-G200-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Hos1_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Hos1_H2AZ_FC1.5_FDR0.001.bed > sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
echo "22hos1narrowpeak done"
##
nohup sh narrowpaeak.sh > narrowpaeaknohup 2>&1 
#
nohup: ignoring input
20686 sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
16clfnarrowpeak done
21442 sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
16colnarrowpeak done
19121 sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
16hos1narrowpeak done
18985 sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
22clfnarrowpeak done
20267 sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
22colnarrowpeak done
19811 sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
22hos1narrowpeak done
#
/mnt/data5/wjc/aaworkflow/chipseqtest/sicer/H2AZ/broadpeak
####
vim broadpaeak.sh
#!/bin/sh
set -e
awk '{if($7>=1.5) print $0}' 16Clf-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Clf_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_16Clf_H2AZ_FC1.5_FDR0.001.bed > sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
#
awk '{if($7>=1.5) print $0}' 16Col-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Col_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_16Col_H2AZ_FC1.5_FDR0.001.bed > sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
#
awk '{if($7>=1.5) print $0}' 16Hos1-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_16Hos1_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_16Hos1_H2AZ_FC1.5_FDR0.001.bed > sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
#
awk '{if($7>=1.5) print $0}' 22Clf-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Clf_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Clf_H2AZ_FC1.5_FDR0.001.bed > sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
#
awk '{if($7>=1.5) print $0}' 22Col-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Col_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Col_H2AZ_FC1.5_FDR0.001.bed > sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
#
awk '{if($7>=1.5) print $0}' 22Hos1-H2A-Z_q20_s_rmdup_chr-W200-G600-islands-summary | awk '{if($8<=0.001) print $0}' >  sig_22Hos1_H2AZ_FC1.5_FDR0.001.bed
awk -F "chr" '{print $2}' sig_22Hos1_H2AZ_FC1.5_FDR0.001.bed > sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
wc -l sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
#
nohup sh broadpaeak.sh > broadpaeaknohup 2>&1
#########
nohup: ignoring input
13982 sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
15344 sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed
13258 sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed
12987 sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed
13950 sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed
13773 sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed


########################################################################################################
#######################################
#######################################################################
##################################################################
#annotation to genebody
peakanno=function(path,filename,samplename) {
  setwd(path)
  suppressMessages(library("ChIPpeakAnno"));
  suppressMessages(library("TxDb.Athaliana.BioMart.plantsmart28"));
  txdb=TxDb.Athaliana.BioMart.plantsmart28;
  annDatagene=genes(txdb);
  peaks <- read.delim(filename,header=T)
  peaks_ranges <- GRanges(seqnames=peaks[,1],ranges=IRanges(start=peaks[,2],end=peaks[,3],names=paste("peak",rep(1:nrow(peaks)),sep="")),strand="*",seqinfo=seqinfo(txdb))
  annoOvgene <- annotatePeakInBatch(peaks_ranges, AnnotationData=annDatagene,output="overlapping",maxgap=0L)
  write.table(annoOvgene, paste(samplename,"_annoOvgene.xls",sep=""), sep="\t",row.names=FALSE)
  annogene <- data.frame(genename=unique(as.character(subset(annoOvgene,annoOvgene$feature != "NA")$feature)))
  write.table(annogene, paste(samplename,"_annoOvgene_genelist.xls",sep=""), sep="\t",row.names=FALSE)
  print(nrow(annogene))  
}
#H2AZ
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed","16Clf_H2AZ_FC1.5_FDR0.001_peaks")
#19399
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_16Col_H2AZ_FC1.5_FDR0.001_nochr.bed","16Col_H2AZ_FC1.5_FDR0.001_peaks")
#19255
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_16Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed","16Hos1_H2AZ_FC1.5_FDR0.001_peaks")
#18152
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_22Clf_H2AZ_FC1.5_FDR0.001_nochr.bed","22Clf_H2AZ_FC1.5_FDR0.001_peaks")
#18215
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_22Col_H2AZ_FC1.5_FDR0.001_nochr.bed","22Col_H2AZ_FC1.5_FDR0.001_peaks")
#18964
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_22Hos1_H2AZ_FC1.5_FDR0.001_nochr.bed","22Hos1_H2AZ_FC1.5_FDR0.001_peaks")
#18605
###############################################这个是macs2的#####################################
peakanno("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\macs",
"C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\macs2peak\\PE\\16Clf-H2A-Z_regularnomodel\\16Clf-H2A-Z_q0.05-regular-nomodel_peaks.narrowPeak"
,"16Clf_H2AZ_mac2narrow_peaks")
#18839

#下面是手动，以H2AZ-16clf-narrow为例
setwd("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow")
suppressMessages(library("ChIPpeakAnno"))
suppressMessages(library("TxDb.Athaliana.BioMart.plantsmart28"))
txdb=TxDb.Athaliana.BioMart.plantsmart28
annDatagene=genes(txdb)
#把每一个基因都拉出来了
peaks <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\sicer\\H2AZ\\narrowpeak\\sig_16Clf_H2AZ_FC1.5_FDR0.001_nochr.bed",header=T)
#macs2也是callpeak结果，header应该是F
seqinfo=seqinfo(txdb)
#如果没有txdb包的话可以自己做一个
seqnames=peaks[,1]
ranges=IRanges(start=peaks[,2],end=peaks[,3],names=paste("peak",rep(1:nrow(peaks)),sep=""))
#把所有的peaks都弄出来了，输入一次就知道了，paste就是简单的连接
peaks_ranges <- GRanges(seqnames=peaks[,1],ranges=IRanges(start=peaks[,2],end=peaks[,3],names=paste("peak",rep(1:nrow(peaks)),sep="")),strand="*",seqinfo=seqinfo(txdb))
annoOvgene <- annotatePeakInBatch(peaks_ranges, AnnotationData=annDatagene,output="overlapping",maxgap=0L)
write.table(annoOvgene, paste("16Clf_H2AZ_FC1.5_FDR0.001_peaks","_annoOvgene.xls",sep=""), sep="\t",row.names=FALSE)
#paste是连接在一起，就是把名字连接在了一起而已
annogene <- data.frame(genename=unique(as.character(subset(annoOvgene,annoOvgene$feature != "NA")$feature)))
#其实这里就已经是把peak关联到的gene拿出来了，但实际有的peak是关联不到基因的，也就是NA，所以直接去掉就好了
write.table(annogene, paste("16Clf_H2AZ_FC1.5_FDR0.001_peaks","_annoOvgene_genelist.xls",sep=""), sep="\t",row.names=FALSE)
print(nrow(annogene))  
#最终能得到一个关联表格和一个genelist表格


#水稻不是很方便，需要自己制造基因组文件，不过发现用原来的chr1文件就行（因为是自己输入的嘛）
#水稻annotation to genebody
peakanno=function(path,filename,samplename) {
  setwd(path)
  suppressMessages(library("ChIPpeakAnno"));
  genebody <- read.delim("E:/复旦/实验室/col6 bed/rice_gene.bed",header=F)
  x <- Seqinfo(seqnames=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","ChrUn","ChrSy"),seqlengths=c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856,633585,592136),isCircular=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),genome="Oryza_sativa")
  annDatagene <- GRanges(seqnames=genebody[,1],ranges=IRanges(start=genebody[,2],end=genebody[,3],names=genebody[,4]),strand=genebody[,6],seqinfo=x)
  peaks <- read.delim(filename,header=F)
  peaks_ranges <- GRanges(seqnames=peaks[,1],ranges=IRanges(start=peaks[,2],end=peaks[,3],names=paste("peak",rep(1:nrow(peaks)),sep="")),strand="*",seqinfo=x)
  annoOvgene <- annotatePeakInBatch(peaks_ranges, AnnotationData=annDatagene,output="overlapping",maxgap=0L)
  write.table(annoOvgene, paste(samplename,"_annoOvgene.xls",sep=""), sep="\t",row.names=FALSE)
  annogene <- data.frame(genename=unique(as.character(subset(annoOvgene,annoOvgene$feature != "NA")$feature)))
  write.table(annogene, paste(samplename,"_annoOvgene_genelist.xls",sep=""), sep="\t",row.names=FALSE)
  print(nrow(annogene))  
}
#区别就是自己制造一个annDatagene而已，注意区分带不带header和有无chr吧
peakanno("E:/复旦/实验室/课题/甲基化/me数据/rice/H3K36me/SIcer/Broad","sig_WT_H3K36me1_FC1.5_FDR0.001.bed","WT_H3K36me1_FC1.5_FDR0.001_peaks")
#18336


#H2AZ IP and Input all genelist RPKMmatrix画pattern
#自动化
ProfilePlot.Body <- function(Path,colMeanDataList,DatanameList,ColorList,mainname,mainnamex,mainnamey,filename,Ylim) {
  setwd(Path)
  pdf(filename)
  Datanumber <- length(colMeanDataList)
  Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
  Xlim <- c(-1000,4000)
  flag <- par(no.readonly=TRUE)
  par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
  library(showtext)
  showtext_auto(enable = TRUE)
  font_add('TNR', 'times.ttf')
  if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:6){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
    for (i in 7:Datanumber){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=2,lwd=2,col=ColorList[i])  
    }
  }
  legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,2,2,2,2,2,2),lwd=2,text.font=1,bty='n',ncol=1)    
  axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
  axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
  text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)
  dev.off() 
} 
all_H2AZ <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\matrix\\all-H2AZ_matrix.mat_nohead", header = F)
row.names(all_H2AZ) <- all_H2AZ[,4]
#把基因名弄成行名
Col_22_H2AZ_input <- all_H2AZ[, 7:506]
Col_22_H2AZ_IP <- all_H2AZ[, 507:1006]
Clf_22_H2AZ_input <- all_H2AZ[, 1007:1506]
Clf_22_H2AZ_IP <- all_H2AZ[, 1507:2006]
Hos1_22_H2AZ_input <- all_H2AZ[, 2007:2506]
Hos1_22_H2AZ_IP <- all_H2AZ[, 2507:3006]
Col_16_H2AZ_input <- all_H2AZ[, 3007:3506]
Col_16_H2AZ_IP <- all_H2AZ[, 3507:4006]
Clf_16_H2AZ_input <- all_H2AZ[, 4007:4506]
Clf_16_H2AZ_IP <- all_H2AZ[, 4507:5006]
Hos1_16_H2AZ_input <- all_H2AZ[, 5007:5506]
Hos1_16_H2AZ_IP <- all_H2AZ[, 5507:6006]
ProfilePlot.Body("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\pattern", list(colMeans(Col_22_H2AZ_IP)/100,colMeans(Clf_22_H2AZ_IP)/100,colMeans(Hos1_22_H2AZ_IP)/100, 
                                            colMeans(Col_16_H2AZ_IP)/100,colMeans(Clf_16_H2AZ_IP)/100,colMeans(Hos1_16_H2AZ_IP)/100, 
                                            colMeans(Col_22_H2AZ_input)/100,colMeans(Clf_22_H2AZ_input)/100, 
                                            colMeans(Hos1_22_H2AZ_input)/100,colMeans(Col_16_H2AZ_input)/100, 
                                            colMeans(Clf_16_H2AZ_input)/100,colMeans(Hos1_16_H2AZ_input)/100), 
                 c("22 Col H2AZ","22 Clf H2AZ","22 Hos1 H2AZ","16 Col H2AZ","16 Clf H2AZ","16 Hos1 H2AZ", 
                 "22 Col input","22 Clf input","22 Hos1 input","16 Col input","16 Clf input","16 Hos1 input"), 
                 c("black","red","orange","green","blue","purple","black","red","orange","green","blue","purple"), 
                 "H2AZ",1500,3.8,"H2AZ IP and Input all genelist RPKMmatrix画pattern.pdf",c(0.5,3.5))

############################################################分步手动###################################################################
setwd("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\pattern")
all_H2AZ <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\matrix\\all-H2AZ_matrix.mat_nohead", header = F)
row.names(all_H2AZ) <- all_H2AZ[,4]
Col_22_H2AZ_input <- all_H2AZ[, 7:506]
Col_22_H2AZ_IP <- all_H2AZ[, 507:1006]
Clf_22_H2AZ_input <- all_H2AZ[, 1007:1506]
Clf_22_H2AZ_IP <- all_H2AZ[, 1507:2006]
Hos1_22_H2AZ_input <- all_H2AZ[, 2007:2506]
Hos1_22_H2AZ_IP <- all_H2AZ[, 2507:3006]
Col_16_H2AZ_input <- all_H2AZ[, 3007:3506]
Col_16_H2AZ_IP <- all_H2AZ[, 3507:4006]
Clf_16_H2AZ_input <- all_H2AZ[, 4007:4506]
Clf_16_H2AZ_IP <- all_H2AZ[, 4507:5006]
Hos1_16_H2AZ_input <- all_H2AZ[, 5007:5506]
Hos1_16_H2AZ_IP <- all_H2AZ[, 5507:6006]
#分别指定谁是谁，1000+1000+2000的和除以10
colMeanDataList <- list(colMeans(Col_22_H2AZ_IP)/100,colMeans(Clf_22_H2AZ_IP)/100,colMeans(Hos1_22_H2AZ_IP)/100,
                                            colMeans(Col_16_H2AZ_IP)/100,colMeans(Clf_16_H2AZ_IP)/100,colMeans(Hos1_16_H2AZ_IP)/100,
                                            colMeans(Col_22_H2AZ_input)/100,colMeans(Clf_22_H2AZ_input)/100,
                                            colMeans(Hos1_22_H2AZ_input)/100,colMeans(Col_16_H2AZ_input)/100,
                                            colMeans(Clf_16_H2AZ_input)/100,colMeans(Hos1_16_H2AZ_input)/100)
#所有基因的列取平均值后除以100，也就在这里决定最后图例的顺序
Datanumber <- length(colMeanDataList)
#length向量有12个组
Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
#把-1000到-10分成100份，把0到2990分成300份，把3000到3900分成100份，那么也就是分成了500份
Xlim <- c(-1000,4000)
Ylim <- c(0.5,3.5)
#xlim是-1000和4000
flag <- par(no.readonly=TRUE)
#我个人感觉应该是防止画不同图出错，但是好像没有大用
par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
#pin：以英寸表示图形的宽和高，mgp指坐标轴与图之前的距离
library(showtext)
showtext_auto(enable = TRUE)
font_add('TNR', 'times.ttf')
#字体包
DatanameList <- c("22 Col H2AZ","22 Clf H2AZ","22 Hos1 H2AZ","16 Col H2AZ","16 Clf H2AZ","16 Hos1 H2AZ", 
                 "22 Col input","22 Clf input","22 Hos1 input","16 Col input","16 Clf input","16 Hos1 input")
ColorList <- c("black","red","orange","green","blue","purple","black","red","orange","green","blue","purple")

##############################################################################################修改参数区域
if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:6){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
    for (i in 7:Datanumber){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=2,lwd=2,col=ColorList[i])  
    }
  }
#####################################################################################修改参数区域
#这是整个画图的核心，主要用的就是R的plot函数，l是绘制线，x轴是500份中的一个点，y轴是得分除以100那个
#如果向量只有一个的话（很明显这里有12个），就正常画就行
#但是如果有超过1个的话，那么第一个还是那么画的，但是2-6个就不一样了，7个及其以后也不一样，ity是虚线还是是实线的意思，lwd是线条的宽度
#第一个为什么和后面的不一样，因为需要用第一个画出来坐标轴的信息，后面直接添加直线就好了
#cex.lab=1.3是字体的大小,bty是边框的意思，tck是刻度线的长度，有正负数，但是实际测试好像没哦什么用，因为后面还有一个tck哈哈哈
#xaxt="n",yaxt="n"表示不显示x轴和y轴的标签
#####################################附加信息修改区域########################
mainname <- "H2AZ"
mainnamex <- 1500
mainnamey <- 3.8
legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,2,2,2,2,2,2),lwd=2,text.font=1,bty='n',ncol=1)
#图例，第一块可以使用"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"，ncol是图例平均分成几列，bty是有无边框
axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
#添加轴线，cex：指定符号的大小。side一个整数，表示在图形的哪边绘制坐标轴（1=下，2=左，3=上，4=右）
text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)
#表示在横纵坐标xxx的地方加上txt文字，font表示字体样式，但是已经指定了TNR，就属实显得比较多余了，xpd就是显示有没有
手动输出文件，就叫做"H2AZ IP and Input all genelist RPKMmatrix画pattern.pdf"
#####################################附加信息修改区域########################


#H2AZ IP 减 Input all genelist RPKMmatrix画pattern
ProfilePlot.Body <- function(Path,colMeanDataList,DatanameList,ColorList,mainname,mainnamex,mainnamey,filename,Ylim) {
  setwd(Path)
  pdf(filename)
  Datanumber <- length(colMeanDataList)
  Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
  Xlim <- c(-1000,4000)
  flag <- par(no.readonly=TRUE)
  par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
  library(showtext)
  showtext_auto(enable = TRUE)
  font_add('TNR', 'times.ttf')
  if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:Datanumber){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
  }
  legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2),lwd=2,text.font=1,bty='n',ncol=1)    
  axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
  axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
  text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)
  dev.off() 
} 
all_H2AZ <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\matrix\\all-H2AZ_matrix.mat_nohead", header = F)
row.names(all_H2AZ) <- all_H2AZ[,4]
Col_22_H2AZ_IP <- all_H2AZ[, 507:1006] - all_H2AZ[, 7:506]
Clf_22_H2AZ_IP <- all_H2AZ[, 1507:2006] - all_H2AZ[, 1007:1506]
Hos1_22_H2AZ_IP <- all_H2AZ[, 2507:3006] - all_H2AZ[, 2007:2506]
Col_16_H2AZ_IP <- all_H2AZ[, 3507:4006] - all_H2AZ[, 3007:3506]
Clf_16_H2AZ_IP <- all_H2AZ[, 4507:5006] - all_H2AZ[, 4007:4506]
Hos1_16_H2AZ_IP <- all_H2AZ[, 5507:6006] - all_H2AZ[, 5007:5506]
ProfilePlot.Body("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\pattern",
list(colMeans(Col_22_H2AZ_IP)/100,colMeans(Clf_22_H2AZ_IP)/100,colMeans(Hos1_22_H2AZ_IP)/100,colMeans(Col_16_H2AZ_IP)/100,
                                  colMeans(Clf_16_H2AZ_IP)/100,colMeans(Hos1_16_H2AZ_IP)/100),
                 c("22 Col H2AZ","22 Clf H2AZ","22 Hos1 H2AZ","16 Col H2AZ","16 Clf H2AZ","16 Hos1 H2AZ"),
                 c("black","red","orange","green","blue","purple"), "H2AZ",1500,1.75,"H2AZ IP 减 Input all genelist RPKMmatrix画pattern.pdf",c(-1,1.5))


#H2AZ IP and Input all H2AZ FC1.5 FDR0.001 peak genelist RPKMmatrix画pattern
ProfilePlot.Body <- function(Path,colMeanDataList,DatanameList,ColorList,mainname,mainnamex,mainnamey,filename,Ylim) {
  setwd(Path)
  pdf(filename)
  Datanumber <- length(colMeanDataList)
  Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
  Xlim <- c(-1000,4000)
  flag <- par(no.readonly=TRUE)
  par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
  library(showtext)
  showtext_auto(enable = TRUE)
  font_add('TNR', 'times.ttf')
  if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:6){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
    for (i in 7:Datanumber){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=2,lwd=2,col=ColorList[i])  
    }
  }
  legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,2,2,2,2,2,2),lwd=2,text.font=1,bty='n',ncol=1)    
  axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
  axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
  text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)
  dev.off() 
} 
all_H2AZ <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\matrix\\all-H2AZ_matrix.mat_nohead", header = F)
row.names(all_H2AZ) <- all_H2AZ[,4]
target <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow\\allgenlistpeak.txt", header = F)
all_H2AZ <- subset(all_H2AZ,is.element(row.names(all_H2AZ),target[,1])==T)  
Col_22_H2AZ_input <- all_H2AZ[, 7:506]
Col_22_H2AZ_IP <- all_H2AZ[, 507:1006]
Clf_22_H2AZ_input <- all_H2AZ[, 1007:1506]
Clf_22_H2AZ_IP <- all_H2AZ[, 1507:2006]
Hos1_22_H2AZ_input <- all_H2AZ[, 2007:2506]
Hos1_22_H2AZ_IP <- all_H2AZ[, 2507:3006]
Col_16_H2AZ_input <- all_H2AZ[, 3007:3506]
Col_16_H2AZ_IP <- all_H2AZ[, 3507:4006]
Clf_16_H2AZ_input <- all_H2AZ[, 4007:4506]
Clf_16_H2AZ_IP <- all_H2AZ[, 4507:5006]
Hos1_16_H2AZ_input <- all_H2AZ[, 5007:5506]
Hos1_16_H2AZ_IP <- all_H2AZ[, 5507:6006]
####################################
ProfilePlot.Body("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\pattern",
                list(colMeans(Col_22_H2AZ_IP)/100,colMeans(Clf_22_H2AZ_IP)/100,colMeans(Hos1_22_H2AZ_IP)/100,
                colMeans(Col_16_H2AZ_IP)/100,colMeans(Clf_16_H2AZ_IP)/100,colMeans(Hos1_16_H2AZ_IP)/100,
                colMeans(Col_22_H2AZ_input)/100,colMeans(Clf_22_H2AZ_input)/100,colMeans(Hos1_22_H2AZ_input)/100,
                colMeans(Col_16_H2AZ_input)/100,colMeans(Clf_16_H2AZ_input)/100,colMeans(Hos1_16_H2AZ_input)/100),
      c("22 Col H2AZ","22 Clf H2AZ","22 Hos1 H2AZ","16 Col H2AZ","16 Clf H2AZ","16 Hos1 H2AZ","22 Col input","22 Clf input",
        "22 Hos1 input","16 Col input","16 Clf input","16 Hos1 input"),
      c("black","red","orange","green","blue","purple","black","red","orange","green","blue","purple"), 
      "H2AZ",1500,4.9,"H2AZ IP and Input all H2AZ FC1.5 FDR0.001 peak genelist RPKMmatrix画pattern.pdf",c(0.5,4.5))

########关键的地方进行分解
#target <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow\\allgenlistpeak.txt", header = F)
#手动去重把所有的peak基因都合在一起，用excel什么的都可以，最后再粘贴到txt里面就好了
#all_H2AZ <- subset(all_H2AZ,is.element(row.names(all_H2AZ),target[,1])==T)  
#row.names(all_H2AZ)
#剩下的就是把所有的基因名字列出来，其实就是老套路取交集了，不懂得话研究RNAseq那里面得就好


#H2AZ IP 减 Input all H2AZ FC1.5 FDR0.001 peak genelist RPKMmatrix画pattern
ProfilePlot.Body <- function(Path,colMeanDataList,DatanameList,ColorList,mainname,mainnamex,mainnamey,filename,Ylim) {
  setwd(Path)
  pdf(filename)
  Datanumber <- length(colMeanDataList)
  Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
  Xlim <- c(-1000,4000)
  flag <- par(no.readonly=TRUE)
  par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
  library(showtext)
  showtext_auto(enable = TRUE)
  font_add('TNR', 'times.ttf')
  if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:Datanumber){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
  }
  legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2),lwd=2,text.font=1,bty='n',ncol=1)    
  axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
  axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
  text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)
  dev.off() 
} 
all_H2AZ <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\upstream\\matrix\\all-H2AZ_matrix.mat_nohead", header = F)
row.names(all_H2AZ) <- all_H2AZ[,4]
target <- read.delim("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\peakanno\\H2AZ\\narrow\\allgenlistpeak.txt", header = F)
all_H2AZ <- subset(all_H2AZ,is.element(row.names(all_H2AZ),target[,1])==T)  
Col_22_H2AZ_IP <- all_H2AZ[, 507:1006] - all_H2AZ[, 7:506]
#Col_22_H2AZ_IP_enrich <- Col_22_H2AZ_IP[rowSums(Col_22_H2AZ_IP)>0,]
#把行求和大于0的拉出来，为什么要这样做，是剔除无效值，理论上IP比input高才行，但是这不是必需的，两种都画一下发现其实区别也不大
Clf_22_H2AZ_IP <- all_H2AZ[, 1507:2006] - all_H2AZ[, 1007:1506]
Hos1_22_H2AZ_IP <- all_H2AZ[, 2507:3006] - all_H2AZ[, 2007:2506]
Col_16_H2AZ_IP <- all_H2AZ[, 3507:4006] - all_H2AZ[, 3007:3506]
Clf_16_H2AZ_IP <- all_H2AZ[, 4507:5006] - all_H2AZ[, 4007:4506]
Hos1_16_H2AZ_IP <- all_H2AZ[, 5507:6006] - all_H2AZ[, 5007:5506]
###############################################################
ProfilePlot.Body("C:\\Users\\94526\\Desktop\\Chipseq_example\\downstream\\pattern",list(colMeans(Col_22_H2AZ_IP)/100,colMeans(Clf_22_H2AZ_IP)/100,
                                    colMeans(Hos1_22_H2AZ_IP)/100,colMeans(Col_16_H2AZ_IP)/100,colMeans(Clf_16_H2AZ_IP)/100,colMeans(Hos1_16_H2AZ_IP)/100),
                 c("22 Col H2AZ","22 Clf H2AZ","22 Hos1 H2AZ","16 Col H2AZ","16 Clf H2AZ","16 Hos1 H2AZ"),
                 c("black","red","orange","green","blue","purple"), "H2AZ",1500,3.35,
                 "H2AZ IP 减 Input all H2AZ FC1.5 FDR0.001 peak genelist RPKMmatrix画pattern.pdf",c(-0.5,3.0))
#如果要画单个的，代码会简单一些
#画单个的相减或者是allgene最简单，直接变一变代码就好了十分简单，可以参照水稻的那个来，或者直接从批量方法去改也行
#但是总体来说，画单个的，还是水稻的那个比较好









all_H3K27 <- read.delim("5weeksH3K27_matrix.mat_nohead", header = F)
row.names(all_H3K27) <- all_H3K27[,4]


target <- read.delim("F:\\WF_ChIP\\annotationresults\\5weeksH3K27\\5weeksH3K27bingji.txt", header = T)
all_H2AZ <- subset(all_H3K27,is.element(row.names(all_H3K27),target[,1])==T)  


Col_H3K27 <- all_H3K27[, 7:506] - all_H3K27[, 507:1006]
elf_H3K27 <- all_H3K27[, 1007:1506] - all_H3K27[, 1507:2006]
elf_ino80_H3K27 <- all_H3K27[, 2007:2506] - all_H3K27[, 2507:3006]
elf_63_H3K27 <- all_H3K27[, 3007:3506] - all_H3K27[, 3507:4006]
elf_63_ino80_H3K27 <- all_H3K27[, 4007:4506] - all_H3K27[, 4507:5006]
ino80_H3K27 <- all_H3K27[, 5007:5506] - all_H3K27[, 5507:6006]
ref6_H3K27 <- all_H3K27[, 6007:6506] - all_H3K27[, 6507:7006]
#Col_H3K27_enrich <- Col_H3K27[rowSums(Col_H3K27)>0,]
#elf_H3K27_enrich <- elf_H3K27[rowSums(elf_H3K27)>0,]
#####


colMeanDataList <- list(colMeans(Col_H3K27)/100,colMeans(elf_H3K27)/100,colMeans(elf_ino80_H3K27)/100, colMeans(elf_63_H3K27)/100,
                        colMeans(elf_63_ino80_H3K27)/100,colMeans(ino80_H3K27)/100,colMeans(ref6_H3K27)/100)


Datanumber <- length(colMeanDataList)

Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))

Xlim <- c(-1000,4000)
Ylim <- c(-1,0.5)

flag <- par(no.readonly=TRUE)
#我个人感觉应该是防止画不同图出错，但是好像没有大用
par(pin=c(5,5),lwd=2,mgp = c(2.5, 0.6, 0))
#pin：以英寸表示图形的宽和高，mgp指坐标轴与图之前的距离
library(showtext)
showtext_auto(enable = TRUE)
font_add('TNR', 'times.ttf')


DatanameList <- c("Col_H3K27","elf_H3K27","elf_ino80_H3K27","elf_63_H3K27","elf_63_ino80_H3K27","ino80_H3K27", "ref6_H3K27")
ColorList <- c("black","red","orange","green","blue","purple","yellow")



if(Datanumber==1) {
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
  }else{
    plot(x=Dot,y=colMeanDataList[[1]],tck=0.02,bty="o",type="l", lty=1, lwd=2,xaxt="n",yaxt="n",col=ColorList[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="Normalized reads density",main="",family = 'TNR', cex.lab=1.3)
    for (i in 2:7){
      lines(x=Dot,y=colMeanDataList[[i]],type="l",lty=1,lwd=2,col=ColorList[i])  
    }
  }
  

mainname <- "H3K27"
mainnamex <- 1500
mainnamey <- 0.6
legend("topright",legend = DatanameList,col=ColorList,lty=c(1,1,1,1,1,1,1),lwd=2,text.font=1,bty='n',ncol=1)
#图例，第一块可以使用"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"，ncol是图例平均分成几列，bty是有无边框
axis(side=1,at=c(-1000,0,3000,4000),labels=c("-1 kb","TSS","TTS","1 kb"),font=1,tck=0.02,lwd=2,cex.axis =1.3)
axis(side=2,lwd=2,tck=0.02,font=1,las=1,cex.axis =1.3)
#添加轴线，cex：指定符号的大小。side一个整数，表示在图形的哪边绘制坐标轴（1=下，2=左，3=上，4=右）
text(x=mainnamex,y=mainnamey, mainname, family = 'TNR',font = 1,cex=1.4,xpd=T)



Peakover <- function(Afilename,Bfilename,peak1name,peak2name,mainname){
 #packages
 suppressMessages(library("ChIPpeakAnno"))
 suppressMessages(library("GenomicRanges"))
 #Afilename-filename of Afile; Bfilename-filename of Bfile; peak1name-name of peaks in Afile; peak2name-name of peaks in Bfile; mainname-prefiex of name in outputfile
 Afile <- read.delim(Afilename,header=F)
 Bfile <- read.delim(Bfilename,header=F)
 Agranges <- GRanges(seqnames=Afile[,1],ranges=IRanges(start=Afile[,2],end=Afile[,3],names=paste("A",rep(1:nrow(Afile)),sep="")))
 Bgranges <- GRanges(seqnames=Bfile[,1],ranges=IRanges(start=Bfile[,2],end=Bfile[,3],names=paste("B",rep(1:nrow(Bfile)),sep="")))
 result <<- findOverlapsOfPeaks(Agranges,Bgranges,connectedPeaks="merge")
 pdf(paste(mainname,".pdf",sep=""))
 makeVennDiagram(result,NameOfPeaks=c(peak1name,peak2name),height=3000,width=3000,col="transparent",fill=c("red","green"),alpha=c(0.5,0.5),main=mainname)
 dev.off()
 write.table(data.frame(result$mergedPeaks)[,1:3],paste(mainname,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F) 
 }



Peakover("C:\\Users\\94526\\Desktop\\3MATERIALS\\R\\peaks\\DongH31.narrowPeak","C:\\Users\\94526\\Desktop\\3MATERIALS\\R\\peaks\\LINH31.narrowPeak",
"Donglab_H3.1_peaks","Linlab_H3.1_peaks","Peaks overlap H3.1")

Peakover("C:\\Users\\94526\\Desktop\\3MATERIALS\\R\\peaks\\DongH33.narrowPeak","C:\\Users\\94526\\Desktop\\3MATERIALS\\R\\peaks\\LINH33.narrowPeak",
"Donglab_H3.3_peaks","Linlab_H3.3_peaks","Peaks overlap H3.3")




#针对双端数据，代码为：
/mnt/data5/wjc/aaworkflow/chipseqtest/macs2peaks/PE
vim macs2callpeakPE.sh
#!/bin/sh
set -e
IP_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup.bam
Input_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup.bam
macs2 callpeak -t $IP_bam -c $Input_bam -f BAMPE -g 1.18e8  -q 0.05 --nomodel  -B --SPMR -n 16Clf-H2A-Z_q0.05-regular-nomodel --outdir 16Clf-H2A-Z_regularnomodel
echo "16Clf-H2A-Z regular nonmodel q0.05"
macs2 callpeak -t $IP_bam -c $Input_bam -f BAMPE -g 1.18e8 --broad --broad-cutoff 0.05 -q 0.05  --nomodel  -B  -n 16Clf-H2A-Z_broadcutoff0.05q0.05-broad-nomodel --outdir 16Clf-H2A-Z_broadnomodel
echo "16Clf-H2A-Z broad nonmodel broadcutoff0.05 q0.05"
echo "all done"
nohup sh macs2callpeakPE.sh > macs2callpeakPEnohup 2>&1 


#针对单端数据，代码为：
/mnt/data5/wjc/aaworkflow/chipseqtest/macs2peaks/SE
vim macs2callpeakSE.sh
#!/bin/sh
set -e
IP1_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-H2A-Z_211120N_S30_L003_q20_s_rmdup.bam
Input_bam=/mnt/data5/wjc/aaworkflow/chipseqtest/Bam/trim1_16Clf-input_211120N_S6_L003_q20_s_rmdup.bam
macs2 callpeak -t $IP1_bam -c $Input_bam -f BAM -g 1.18e8  -q 0.05  -B --SPMR -n 16Clf-H2A-Z_q0.05-regular-model --outdir 16Clf-H2A-Z_regularmodel
echo "16Clf-H2A-Z regular model q0.05"
macs2 callpeak -t $IP1_bam -c $Input_bam -f BAM -g 1.18e8 --broad --broad-cutoff 0.05 -q 0.05  -B  -n 16Clf-H2A-Z_broadcutoff0.05q0.05-broad-model --outdir 16Clf-H2A-Z_broadmodel
echo "16Clf-H2A-Z broad model broadcutoff0.05 q0.05"
echo "all done"
nohup sh macs2callpeakSE.sh > macs2callpeakSEnohup 2>&1 




#############以上所有的代码都不要加 --SPMR       ！！！！！！！！！！！！！！！！！！！！！！！！！！！
#原因如下


# -B的意思是输出两个bedgraph文件（bdg），实际上和bigwig文件差不多，其中一个是treat组，另一个是control组,但是bdg文件并不是经过标准化的，是原始的堆积值
#如果不加-B，就没有这两个beddgraph文件，这两个bdg文件的第一列是染色体号，第二列和第三列是起始和终止位点，第四列是值得注意的点！
#第四列是原始的这个位点堆积了多少条reads的值（后面有详细说明）
#因为
#######--SPMR的意思实际是什么呢?
#同一组IP和input的测序深度实际上是不同的，那怎么办呢？经过仔细研究了一番，发现输出的两个bdg文件的第四列有蹊跷
#为什么呢？因为IP组永远都是整数，而input组永远都是小数，其实是因为它会把input组的那个原始的reads堆积值根据IP组乘以或者除以某个系数来实现这两组
#测序深度带的统一，我试了一下，又拿另一个组的IP2和input2去做了一下，发现也是这样的，而且这个缩放的系数是不一样的，也证明了确实是这样
#我后来还读了软件的说明，是这样的
#实际上，如果你加了--SPMR，那么他就会直接消除测序深度的影响，不像上面的那个只是缩放两组，而是用百万reads来直接进行类似于RPKM后的操作
#那么，如果加了SPMR，那么也就证明会把所有的组的IP和input都缩放到同一个水平，我导入了IGV发现的确是这样
#如果每家SPMR，那么其实就是类似于原始的堆积（pileup），而加了SPMR也就是给标准化了，也就是所有的东西测序深度的影响就消除啦！
#但是为什么不能加SPMR呢？第一是因为如果仔细看软件的说明的话，就会发现SPMR不能和后续的bdgdiff一起用
#其实是因为标准化以后的第四列太小了，基本上都是小于1的数，所以大家其实都差不多大，算法就没有办法去找到差异peaks了
#我试过，加了SPMR只能找出来几个差异peaks，而不加的话则很正常，我还查了资料，就是这样的



#######说了一大堆，其实就是不要加--SPMR才能做后续的差异peak分析
######刘翻和小浩用的差异peaks分析和我的不一样，首先是因为我们用macs2 call的peaks，所以我喜欢用macs2的bdgdiff去做差异peaks分析
###小浩喜欢用sicer callpeak，所以他喜欢用sicer去call peak，而刘翻使用的是别的软件，叫diffbind

macs2 bdgdiff  --t1 t1IP.bdg --t2 t2IP.bdg --c1 c1input.bdg --c2 c2input,bdg  --d1 DEPTH1 --d2 DEPTH2 --outdir 自己搞的文件夹xxxxx --o-prefix 直接起的名字xxxxx
#d1实际上就hi第一组IP的深度，因为input前面说被IP给抵消了，d2就是第二组IP的深度，要是懒得看的话，教你一个小妙招，直接去看bam文件的大小，然后
#比方说，第一个IP500mb，第二个IP800mb，那么d1就是5，d2就是8，其实在callpeak的时候也能看到，自己算也行

[-C CUTOFF] [-l MINLEN] [-g MAXGAP]
#可选选项，第一个是过滤，因为最后输出的那个差异性也有个类似于p值得东西，取对数还是什么的，反正最后默认是保留3以上的，
#3以下的就是不显著的，也就是common
#第二个是少于多少的长度就不要了，默认值是200，其实网上有人说甲基化什么的用147什么的，而且可以用macs2 predictd -i xxxx.bam去看
#但是我觉得不用，直接不要管，默认参数就好了啊
#第三个更没用，默认默认

#(--o-prefix 名字xxxxx | -o OFILE OFILE OFILE)
#其实输出的文件名称又两种选择，第一种是会生成三个后缀差不多的bed，1是第一组比第二组高的，2是低的
#另外就是可以直接指定名字，随你的便

#最后得到的差异peak可以再去关联基因，再画一个pattern也行，或者去关联基因，去做GO也行


##########################双端片段长度
#!/bin/sh
set -e
#bam files
WTH2AZ_IP1=WT1_H2AZ_repeat1_q20_s_rm.bam              
WTH3_IP1=WT1_H3_repeat1_q20_s_rm.bam                  
WT_in1=WT1_input_repeat1_q20_s_rm.bam                 
muH2AZ_IP1=ino80_H2AZ_repeat1_q20_s_rm.bam            
muH3_IP1=ino80_H3_repeat1_q20_s_rm.bam                
mu_in1=ino80_input_repeat1_q20_s_rm.bam               
WTH2AZ_IP2=WT1_H2AZ_repeat2_q20_s_rm.bam              
WTH3_IP2=WT1_H3_repeat2_q20_s_rm.bam                
WT_in2=WT1_input_repeat2_q20_s_rm.bam              
muH2AZ_IP2=ino80_H2AZ_repeat2_q20_s_rm.bam            
muH3_IP2=ino80_H3_repeat2_q20_s_rm.bam                
mu_in2=ino80_input_repeat2_q20_s_rm.bam   
#H2AZ
bamPEFragmentSize  --bamfiles  ${WTH2AZ_IP1} ${WTH2AZ_IP2}  ${muH2AZ_IP1} ${muH2AZ_IP2}  \
--histogram H2AZ_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel WTH2AZ_IP1 WTH2AZ_IP2  muH2AZ_IP1 muH2AZ_IP2 --maxFragmentLength 500  
echo "1done"
#H3
bamPEFragmentSize  --bamfiles  ${WTH3_IP1} ${WTH3_IP2}  ${muH3_IP1} ${muH3_IP2}  \
--histogram H3_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel  WTH3_IP1 WTH3_IP2 muH3_IP1 muH3_IP2  --maxFragmentLength 500 
echo "2done"
#input
bamPEFragmentSize  --bamfiles  ${WT_in1} ${WT_in2}  ${mu_in1} ${mu_in2} \
--histogram input_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel WT_in1 WT_in2 mu_in1 mu_in2   --maxFragmentLength 500 
echo "3done"      

setwd(dir = "C:\\Users\\94526\\Desktop\\rawcount")
suppressMessages(library(DESeq2))
sample1  <- read.delim("asf1ab-1-210906S_combined_q20_s.rawcount")
sample2  <- read.delim("asf1ab-2-210906S_combined_q20_s.rawcount")
sample3  <- read.delim("asf1ab-3-210906S_combined_q20_s.rawcount")
sample4  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
sample5  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
sample6  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
rawcount <- data.frame(row.names = sample1[,1], 
                       asf1ab_1 = sample1[,2], asf1ab_2 = sample2[,2], asf1ab_3 = sample3[,2],
                       col_1 = sample4[,2], col_2 = sample5[,2], col_3 = sample6[,2])
#在这里就能看到差异表达的本质
condition <- factor(rep(c("asf1ab", "col"), each = 3))
dds <- DESeqDataSetFromMatrix(rawcount, DataFrame(condition), design = ~condition)
dds_filter <- dds[rowSums(counts(dds)) > 1,]
#行求和大于0的得到保留，因为许多基因是在几个样本中都不表达的
#dds实际上也就是存放原始数据而已
dds_DES <- DESeq(dds_filter)
#算的比较久，在这里完成计算
res1 <- results(dds_DES, contrast = c("condition", "asf1ab", "col"))
#决定16clf比16col
#这里是输出了所有的差异基因


genebed <- read.delim("C:\\Users\\94526\\Desktop\\rawcount\\TAIR10_col6.bed",header=F)
#这个genebed的数量和count的基因数目正好一致，也就是说有多少gene_id就有多少genebed基因，注意配套就好了
row.names(genebed) <- genebed[,4]
#gene名充当行名
genebed <- genebed[as.character(row.names(data.frame(res1))),]
#只保留genelist中的差异基因，其他的不要了,因为有许多基因的表达量是0 0 0 0 0 0，所以genelist中的基因是所有基因的真子集
#data<-na.omit(data)
#去掉含有NA的行，如果出了这个问题请注意以下
#data<-subset(data,value!="NA")
#去掉value行的NA
gene<- GRanges(seqnames=genebed[,1],ranges=IRanges(start=genebed[,2],end=genebed[,3]),names=genebed[,4],strand=genebed[,6])
#构建GRanges对象，所有东西都是genebed里面来的,在这里就会有基因的长度，基因从哪里开始的，到哪里去结束
#做ChIPseq的peakoverlap和annotation都能用到这个GRanges
#很多bioconductor都是用这个包去做的，在加载DESeq2的时候会协助加载这个包GenomicRanges
gene_sp <- split(gene, factor(elementMetadata(gene)$names,levels=as.character(row.names(data.frame(res1)))))
#factor将其变换成因子，重新赋值，实际上就是把引号去掉了，方便与进行后续操作
#GenomicRanges包里的专用方法，把每一个基因去分开了
#其实到这里，只是借助了res1的基因list，并不涉及res1里面的内容
rowRanges(dds_filter) <- gene_sp
#实际上，我也并不是很理解这个的意思，都是DESeq2去提供的
FPKM <- data.frame(fpkm(dds_filter))
#这个fpkm是DESeq2提供的，就是这样用的
write.table(FPKM,"FPKM.xls",sep="\t")
#就像LCZ说得一样，在进行差异分析时候仅仅是消除测序深度的影响，并不会涉及到FPKM或者是TPM，也就是不在乎基因的长度，因为一对一去比较
#但是DESeq2的确可以做fpkm


diff1 <- subset(res1, abs(res1$log2FoldChange) >= log2(1.5) & res1$pvalue <= 0.05)
#变化倍数大于1.5倍的，p值小于0.05的
up1 <- subset(res1, res1$log2FoldChange >= log2(1.5) & res1$pvalue <= 0.05)
down1 <- subset(res1, res1$log2FoldChange <= -log2(1.5) & res1$pvalue <= 0.05)
write.table(res1, file = "asf1ab_vs_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#全部
write.table(diff1, file = "asf1ab_diff_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#差异基因
write.table(up1, file = "asf1ab_up_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#上调基因
write.table(down1, file = "asf1ab_down_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#下调基因
nrow(up1); nrow(down1)


setwd(dir = "C:\\Users\\94526\\Desktop\\rawcount")
suppressMessages(library(DESeq2))
sample1  <- read.delim("asf1ab-1-210906S_combined_q20_s.rawcount")
sample2  <- read.delim("asf1ab-2-210906S_combined_q20_s.rawcount")
sample3  <- read.delim("asf1ab-3-210906S_combined_q20_s.rawcount")
sample4  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
sample5  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
sample6  <- read.delim("Col-1-210906S_combined_q20_s.rawcount")
rawcount <- data.frame(row.names = sample1[,1], 
                       asf1ab_1 = sample1[,2], asf1ab_2 = sample2[,2], asf1ab_3 = sample3[,2],
                       col_1 = sample4[,2], col_2 = sample5[,2], col_3 = sample6[,2])
#在这里就能看到差异表达的本质
condition <- factor(rep(c("asf1ab", "col"), each = 3))
dds <- DESeqDataSetFromMatrix(rawcount, DataFrame(condition), design = ~condition)
dds_filter <- dds[rowSums(counts(dds)) > 1,]
#行求和大于0的得到保留，因为许多基因是在几个样本中都不表达的
#dds实际上也就是存放原始数据而已
dds_DES <- DESeq(dds_filter)
#算的比较久，在这里完成计算
res1 <- results(dds_DES, contrast = c("condition", "asf1ab", "col"))
#决定16clf比16col
#这里是输出了所有的差异基因


genebed <- read.delim("C:\\Users\\94526\\Desktop\\rawcount\\TAIR10_col6.bed",header=F)
#这个genebed的数量和count的基因数目正好一致，也就是说有多少gene_id就有多少genebed基因，注意配套就好了
row.names(genebed) <- genebed[,4]
#gene名充当行名
genebed <- genebed[as.character(row.names(data.frame(res1))),]
#只保留genelist中的差异基因，其他的不要了,因为有许多基因的表达量是0 0 0 0 0 0，所以genelist中的基因是所有基因的真子集
#data<-na.omit(data)
#去掉含有NA的行，如果出了这个问题请注意以下
#data<-subset(data,value!="NA")
#去掉value行的NA
gene<- GRanges(seqnames=genebed[,1],ranges=IRanges(start=genebed[,2],end=genebed[,3]),names=genebed[,4],strand=genebed[,6])
#构建GRanges对象，所有东西都是genebed里面来的,在这里就会有基因的长度，基因从哪里开始的，到哪里去结束
#做ChIPseq的peakoverlap和annotation都能用到这个GRanges
#很多bioconductor都是用这个包去做的，在加载DESeq2的时候会协助加载这个包GenomicRanges
gene_sp <- split(gene, factor(elementMetadata(gene)$names,levels=as.character(row.names(data.frame(res1)))))
#factor将其变换成因子，重新赋值，实际上就是把引号去掉了，方便与进行后续操作
#GenomicRanges包里的专用方法，把每一个基因去分开了
#其实到这里，只是借助了res1的基因list，并不涉及res1里面的内容
rowRanges(dds_filter) <- gene_sp
#实际上，我也并不是很理解这个的意思，都是DESeq2去提供的
FPKM <- data.frame(fpkm(dds_filter))
#这个fpkm是DESeq2提供的，就是这样用的
write.table(FPKM,"FPKM.xls",sep="\t")
#就像LCZ说得一样，在进行差异分析时候仅仅是消除测序深度的影响，并不会涉及到FPKM或者是TPM，也就是不在乎基因的长度，因为一对一去比较
#但是DESeq2的确可以做fpkm


diff1 <- subset(res1, abs(res1$log2FoldChange) >= log2(1.5) & res1$pvalue <= 0.05)
#变化倍数大于1.5倍的，p值小于0.05的
up1 <- subset(res1, res1$log2FoldChange >= log2(1.5) & res1$pvalue <= 0.05)
down1 <- subset(res1, res1$log2FoldChange <= -log2(1.5) & res1$pvalue <= 0.05)
write.table(res1, file = "asf1ab_vs_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#全部
write.table(diff1, file = "asf1ab_diff_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#差异基因
write.table(up1, file = "asf1ab_up_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#上调基因
write.table(down1, file = "asf1ab_down_col_p0.05_FC(1.5).xls", sep = "\t", row.names = T)
#下调基因
nrow(up1); nrow(down1)


stringtie /mnt/data5/wjc/aaworkflow/asf1ab_col_fas2_hira2_RNAseq/Bam/hira_1_q20_s.bam -A hira_1.txt -e -p 8 \
-G /mnt/data5/wjc/genome_scripts_software/genome/ara/TAIR10/TAIR10_GTF.gtf --rf -o hira_1.gtf 2>stingtie.log
#基础环境的软件有问题，后面的这个2>是输出结果，好用的得很

setwd("C:/Users/94526/Desktop/hira/stingtieFPKMTPM")
hira_1 <- read.delim("hira_1.txt")
row.names(hira_1) <- hira_1[,1]
hira_1_order <- hira_1 [order(row.names(hira_1),decreasing=F),]
write.table(hira_1_order, file = "hira_1.xls", sep = "\t", row.names = F)
#Stringtie结果排序

#处理异常表格
wf <- data.frame(row.names = sample1[,1], rep1= sample1[,2] ,rep2= sample1[,3],rep3= sample1[,4],mean=sample1[,5])
write.table(wf, file = "Col_RNAseq_FPKM.xls", sep = "\t", row.names = T)



#把bigwig文件拆分为正负链
bamCoverage -p 20 -b Col-1_q20_s.bam -o Col-1_q20_s_RPKM_R.bigwig --filterRNAstrand=reverse --binSize 10 --normalizeUsing RPKM



scp -r /data1/wjc/genome/TAIR10 jccity@10.157.34.11:/mnt/data5/wjc
jccity@10.157.34.11's password:
在不同服务器之间传输数据






