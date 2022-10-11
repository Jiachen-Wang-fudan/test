#!/bin/sh
set -e
#1.quality control & trimming
mkdir QC
fastqc *.gz -o QC
ls -1 *.fastq.gz |awk -F ".fastq.gz" '{print $1}' | while read id
do
echo "${id}"
trim_galore --phred33 --fastqc -q 20  --stringency 3 --length 20 -o trimGalore_trim_s3_l20  ${id}.fastq.gz  -j 4
done
echo "trim done"
mkdir rawfastq ; mv *.gz rawfastq ; mv QC rawfastq ; cd trimGalore_trim_s3_l20 ; mv *.gz .. ; cd ..
ls -1 *_trimmed.fq.gz |awk -F "_trimmed.fq.gz" '{print $1}' | while read id
do
echo "${id}"
mv "${id}"_trimmed.fq.gz   "${id}".fq.gz
done
echo "rename done"
#2.mapping & data cleaning
ls -1 *.fq.gz |awk -F ".fq.gz" '{print $1}' | while read id
do
echo "${id}"
bowtie2 -p 4 -x /mnt/data5/wjc/genome/ara/fixed_TAIR10/bowtie2_index/bowtie2_index -U ${id}.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'chloroplast' | grep -v 'mitochondria' > ${id}_rmptmt.sam
echo "rmptmt done"
samtools view -@ 4 -q 20 -bS ${id}_rmptmt.sam > ${id}_rmptmt_q20.bam
echo "q20 bam done"
samtools sort -@ 4 ${id}_rmptmt_q20.bam -o ${id}_rmptmt_q20_s.bam
echo "sort done"
samtools rmdup -@ 4 -s ${id}_rmptmt_q20_s.bam ${id}_rmptmt_q20_s_rmdup.bam
echo "rmdup done"
samtools index -@ 4 ${id}_rmptmt_q20_s_rmdup.bam
echo "index done"
samtools view -F 4 -@ 4 -c ${id}.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt_q20.bam
samtools view -F 4 -@ 4 -c ${id}_rmptmt_q20_s_rmdup.bam
echo "sam_count"
echo "rmptmt_count"
echo "q20_count"
echo "rmdup_count"
echo "reads_count_done"
rm ${id}_rmptmt_q20_s.bam ; rm ${id}_rmptmt_q20.bam ; rm ${id}_rmptmt.sam ; rm ${id}.sam
done
mkdir Fastq4map ; mv *.fq.gz Fastq4map
echo "finalbam done"
#3.normalization & BAM2BED
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -p 4 -b ${id}.bam -o ${id}_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
bamCoverage -e 150 -p 4 -b ${id}.bam -o ${id}_RPKM_e.bigwig --binSize 10 --normalizeUsing RPKM
bamCoverage -e 200 -p 4 -b ${id}.bam -o ${id}_RPKM_e200.bigwig --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"
bamCoverage -e 150 -p 4 -b ${id}.bam -o ${id}_BPM.bigwig --binSize 10 --normalizeUsing BPM
echo "BPM_TPM bigwig done"
bamCoverage -e 150 -p 4 -b ${id}.bam -o ${id}_RPGC.bigwig --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 119481543
echo "RPGC bigwig done"
bedtools bamtobed -i ${id}.bam > ${id}.bed
echo "bed done"
awk '{print "chr"$0}'  ${id}.bed > ${id}_chr.bed
echo "Add chr done"
rm ${id}.bed
done
mkdir Bam ; mv *.bam Bam ; mv *.bam.bai Bam ; mkdir RPKMbigwig ; mv *RPKM.bigwig RPKMbigwig ; mkdir RPKMbigwig_e
mv *RPKM_e.bigwig RPKMbigwig_e ; mkdir RPKMbigwig_e200 ; mv *RPKM_e200.bigwig RPKMbigwig_e200
mkdir Bed ; mv *.bed Bed ; mkdir BPMbigwig ; mv *BPM.bigwig BPMbigwig ; mkdir RPGCbigwig ; mv *RPGC.bigwig RPGCbigwig
mkdir track ; mv BPMbigwig track ; mv RPGCbigwig track ; mv RPKMbigwig track ; mv RPKMbigwig_e track ; mv RPKMbigwig_e200 track
echo "bigwig done"
#4.matrix-build
##4.1RPKM_e
cd ./track/RPKMbigwig_e
ls *.bigwig >bw_name_e.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name_e.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPKM_e.bigwig" '{print $1}' >bw_name2.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name2.txt
cat bw_name_e.txt |  while read id
do
echo "${id}"
cat bw_name2.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 4 -R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 4 -out all_matrix_e.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix_e.mat.gz -out all_matrix_e_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix_e.mat.gz   -out all_matrix_e_pattern.pdf  \
--perGroup --legendLocation upper-right  \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
mkdir allmatrix_e ; mv *.pdf allmatrix_e ; mv all_matrix_e.mat.gz allmatrix_e
##4.2RPKM
cd ../.. ; cd ./track/RPKMbigwig ; ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPKM.bigwig" '{print $1}' >bw_name3.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name3.txt
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name3.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 4 -R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 4 -out all_matrix.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix.mat.gz -out all_matrix_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix.mat.gz   -out all_matrix_pattern.pdf  \
--perGroup --legendLocation upper-right \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
echo "matrix done"
mkdir allmatrix ; mv *.pdf allmatrix ; mv all_matrix.mat.gz allmatrix
##4.3BPM_e
cd ../.. ; cd ./track/BPMbigwig ; ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_BPM.bigwig" '{print $1}' >bw_name4.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name4.txt
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name4.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 4 -R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 4 -out all_matrix_e_BPM.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix_e_BPM.mat.gz -out all_matrix_e_BPM_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix_e_BPM.mat.gz   -out all_matrix_e_BPM_pattern.pdf  \
--perGroup --legendLocation upper-right \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
echo "matrix done"
mkdir allmatrix ; mv *.pdf allmatrix ; mv all_matrix_e_BPM.mat.gz allmatrix
##4.4RPGC
cd ../.. ; cd ./track/RPGCbigwig ; ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPGC.bigwig" '{print $1}' >bw_name5.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name5.txt
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name5.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 4 -R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 4 -out all_matrix_e_RPGC.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix_e_RPGC.mat.gz -out all_matrix_e_RPGC_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix_e_RPGC.mat.gz   -out all_matrix_e_RPGC_pattern.pdf  \
--perGroup --legendLocation upper-right \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
echo "matrix done"
mkdir allmatrix ; mv *.pdf allmatrix ; mv all_matrix_e_RPGC.mat.gz allmatrix
##4.5RPKM_e200
cd ../.. ; cd ./track/RPKMbigwig_e200 ; ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPKM_e200.bigwig" '{print $1}' >bw_name3.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name3.txt
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name3.txt |  while read slabel
do
echo "${slabel}"
computeMatrix scale-regions \
-S ${id} \
-p 4 -R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero -p 4 -out all_matrix_e200.mat.gz
echo "matrix1 done"
plotHeatmap -m all_matrix_e200.mat.gz -out all_matrix_e200_heatmap.pdf --sortRegions keep --colorList "white,#77AC30"
echo "plot1 done"
plotProfile -m all_matrix_e200.mat.gz   -out all_matrix_e200_pattern.pdf  \
--perGroup --legendLocation upper-right \
--colors grey black blue pink darkred red yellow green purple brown olive orange "#FF00FF" "#CCFB5D" "#46C7C7" "#000080" "#56A5EC" \
--samplesLabel ${slabel}
done
done
echo "plot2 done"
echo "matrix done"
mkdir allmatrix ; mv *.pdf allmatrix ; mv all_matrix_e200.mat.gz allmatrix
echo "all kinds of patterns done"
cd ../.. ; mv trimGalore_trim_s3_l20 Fastq4map ; mkdir fq ; mv Fastq4map fq ; mv rawfastq fq
echo "all of Single-ChIP-seq done" 