#!/bin/sh
set -e
#1.Quality control & trimmings
mkdir QC
fastqc *.gz -o QC
ls -1 *_R1.fastq.gz |awk -F "_R1.fastq.gz" '{print $1}' | while read id
do
echo "${id}"
trim_galore --phred33 --fastqc -q 20  --stringency 3 --length 20 -o trimGalore_trim_s3_l20 --paired ${id}_R1.fastq.gz ${id}_R2.fastq.gz -j 4
done
echo "trim done"
mkdir rawfastq ; mv *.gz rawfastq ; mv QC rawfastq ; cd trimGalore_trim_s3_l20 ; mv *.gz .. ; cd ..
ls -1 *_R1_val_1.fq.gz |awk -F "_R1_val_1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
mv "${id}"_R1_val_1.fq.gz "${id}"_R1.fq.gz
mv "${id}"_R2_val_2.fq.gz "${id}"_R2.fq.gz
done
echo "rename done"
#2.Mapping & data cleaning
ls -1 *_R1.fq.gz |awk -F "_R1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
hisat2 -p 4 --dta --rna-strandness RF -x /mnt/data5/wjc/genome/ara/fixed_TAIR10/hisat2_index/hisat2_index -1 ${id}_R1.fq.gz -2 ${id}_R2.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'chloroplast' | grep -v 'mitochondria' > ${id}_rmptmt.sam
echo "rmptmt done"
samtools view -q 20 -@ 4 -bS ${id}_rmptmt.sam > ${id}_rmptmt_q20.bam
echo "q20 bam done"
samtools sort -@ 4 ${id}_rmptmt_q20.bam -o ${id}_rmptmt_q20_s.bam
echo "sort done"
samtools index -@ 4 ${id}_rmptmt_q20_s.bam
echo "index done"
samtools view -F 4 -@ 4 -c ${id}.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt_q20.bam
samtools view -F 4 -@ 4 -c ${id}_rmptmt_q20_s.bam
echo "reads count done"
rm ${id}.sam ; rm ${id}_rmptmt.sam ; rm ${id}_rmptmt_q20.bam
echo "sam_count"
echo "rmptmt_count"
echo "q20_count"
echo "sort_count"
done
mkdir Fastq4map ; mv *.fq.gz Fastq4map
echo "finalbam done"
#3.Normalization
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -p 4 -b ${id}.bam -o ${id}_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"
bamCoverage -p 4 -b ${id}.bam -o ${id}_BPM.bigwig --binSize 10 --normalizeUsing BPM
echo "BPM_TPM bigwig done"
done
mkdir Bam ; mv *.bam Bam ; mv *.bam.bai Bam ; mkdir RPKMbigwig ; mv *RPKM.bigwig RPKMbigwig ; mkdir BPMbigwig ; mv *BPM.bigwig BPMbigwig
mkdir track ; mv BPMbigwig track ; mv RPKMbigwig track
echo "bigwig done"
#4.Fragmentsize
cd ./Bam ; ls *.bam >bam_name.txt
sed -i ':label;N;s/\n/ /;b label'  bam_name.txt
ls -1 *.bam | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup.bam" '{print $1}' >bam_name1.txt
sed -i ':label;N;s/\n/ /;b label'  bam_name1.txt
cat bam_name.txt |  while read id
do
echo "${id}"
cat bam_name1.txt |  while read slabel
do
echo "${slabel}"
bamPEFragmentSize  --bamfiles ${id} --histogram all_bam_fragsize.pdf  --numberOfProcessors 4 --maxFragmentLength 500 --samplesLabel ${slabel}
done
done
mkdir fragmentsize ; mv *.pdf fragmentsize
echo "fragmentsize done" 
#5.Correlation using Bam
cat bam_name.txt |  while read id
do
echo "${id}"
multiBamSummary bins -bs 1000 -p 4 --bamfiles ${id} -out readCounts1.npz --outRawCounts readCounts1.tab
done
echo "corrlation matrix done"
plotCorrelation -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_PearsonCorr1.png --outFileCorMatrix PearsonCorr1.tab
echo "1 done"
plotCorrelation -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_PearsonCorr2.png   --outFileCorMatrix PearsonCorr2.tab
echo "2 done"
plotCorrelation -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_SpearmanCorr1.png --outFileCorMatrix SpearmanCorr1.tab
echo "3 done"
plotCorrelation -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_SpearmanCorr2.png --outFileCorMatrix SpearmanCorr2.tab
echo "4 done"
plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_PearsonCorr3.png --outFileCorMatrix PearsonCorr3.tab
echo "1z done"
plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_PearsonCorr4.png --outFileCorMatrix PearsonCorr4.tab
echo "2z done"
plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot scatterplot --removeOutliers -o scatterplot_SpearmanCorr3.png --outFileCorMatrix SpearmanCorr3.tab
echo "3z done"
plotCorrelation --zMin 0.95 -in readCounts1.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap coolwarm --plotNumbers --removeOutliers -o heatmap_SpearmanCorr4.png --outFileCorMatrix SpearmanCorr4.tab
echo "4z done"
echo "plot correlation done"
mkdir Correlationplot ; mv *.tab Correlationplot ; mv *.png Correlationplot ; mv *.npz Correlationplot ; mv Correlationplot ..
#5.count
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
featureCounts -T 4 -p -B -C -s 2 -t exon -g gene_id -F GTF -a /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_exon.gtf -o ${id}.count ${id}.bam
echo "count done"
sed '1d' ${id}.count | awk '{print $1"\t"$7}' > ${id}.rawcount
echo "rawcount done"
done
mkdir count ; mkdir rawcount ; mv *.count count ; mv *.rawcount rawcount ; mv *.summary count ; mv count .. ; mv rawcount .. ; cd .. ; mkdir featurecount 
mv count featurecount ; mv rawcount featurecount ; cd featurecount ; cd count ; mkdir nouse_summary ; mv *.summary nouse_summary ; cd ../.. ; cd Bam
echo "count done"
#6.Strandness infer
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
infer_experiment.py -r /home/others/jccity/genome/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed -i ${id}.bam
done
echo "bed6_infer done"
cd .. ; mv trimGalore_trim_s3_l20 Fastq4map ; mkdir fq ; mv Fastq4map fq ; mv rawfastq fq
echo "all of RNAseq done" 