#!/bin/sh
set -e
#首先要准备好fastq的压缩文件，该脚本仅仅针对于双端测序chipseq
#拿到数据以后先改名字，改名模板应该改为
#xxxxxx_R1.fastq.gz
#xxxxxx_R2.fastq.gz
#有重复的话填上就行，没有的话就不用
#名字里面不要有短横线-，一律改成下划线_
#相应的matrix应该在bigwig的文件夹下找(track)
#一个log文件，在里面搜索findme，然后就能够找到mapping率
#1.做质控，然后去接头
mkdir QC
fastqc *.gz -o QC
ls -1 *_R1.fastq.gz |awk -F "_R1.fastq.gz" '{print $1}' | while read id
do
echo "${id}"
trim_galore --phred33 --fastqc -q 20  --stringency 3 --length 20 -o trimGalore_trim_s3_l20 --paired ${id}_R1.fastq.gz ${id}_R2.fastq.gz -j 4
done
echo "trim done"

mkdir rawfastq
mv *.gz rawfastq
mv QC rawfastq

cd trimGalore_trim_s3_l20
mv *.gz ..
cd ..

ls -1 *_R1_val_1.fq.gz |awk -F "_R1_val_1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
mv "${id}"_R1_val_1.fq.gz   "${id}"_R1.fq.gz
mv "${id}"_R2_val_2.fq.gz   "${id}"_R2.fq.gz
done
echo "rename done"






#2.比对，去叶绿体、去q20、去重复、统计数目
ls -1 *_R1.fq.gz |awk -F "_R1.fq.gz" '{print $1}' | while read id
do
echo "${id}"
bowtie2 -p 4 -x /mnt/data5/wjc/genome/ara/fixed_TAIR10/bowtie2_index/bowtie2_index -1 ${id}_R1.fq.gz -2 ${id}_R2.fq.gz -S ${id}.sam
echo "map done"
cat ${id}.sam | grep -v 'chloroplast' | grep -v 'mitochondria' >${id}_rmptmt.sam
echo "rmptmt done"
samtools view -q 20 -@ 4 -bS ${id}_rmptmt.sam > ${id}_rmptmt_q20.bam
echo "q20 bam done"
samtools sort -@ 4 -n ${id}_rmptmt_q20.bam -o ${id}_rmptmt_q20_s1.bam
echo "nsort1 done"
samtools fixmate -@ 4 -m ${id}_rmptmt_q20_s1.bam ${id}_rmptmt_q20_s1_fixmate.bam
echo "fixmate done"
samtools sort -@ 4 ${id}_rmptmt_q20_s1_fixmate.bam -o ${id}_rmptmt_q20_s1_fixmate_s2.bam
echo "sort2 done"
samtools markdup -@ 4 -r ${id}_rmptmt_q20_s1_fixmate_s2.bam ${id}_rmptmt_q20_s1_fixmate_s2_markdup.bam
echo "markdup done"
samtools index -@ 4 ${id}_rmptmt_q20_s1_fixmate_s2_markdup.bam
echo "index done"

samtools view -F 4 -@ 4 -c ${id}.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt.sam
samtools view -F 4 -@ 4 -c ${id}_rmptmt_q20.bam

#samtools view -F 4 -c ${id}_rmptmt_q20_s1.bam
#echo "q20_s1_count"
#samtools view -F 4 -c ${id}_rmptmt_q20_s1_fixmate.bam
#echo "q20_s1_fixmate_count"
#samtools view -F 4 -c ${id}_rmptmt_q20_s1_fixmate_s2.bam
#echo "q20_s1_fixmate_s2_count"

samtools view -F 4 -c ${id}_rmptmt_q20_s1_fixmate_s2_markdup.bam

echo "sam_count"
echo "rmptmt_count"
echo "q20_count"
echo "markdup_count"
rm ${id}_rmptmt_q20_s1_fixmate_s2.bam
rm ${id}_rmptmt_q20_s1_fixmate.bam
rm ${id}_rmptmt_q20_s1.bam
rm ${id}_rmptmt_q20.bam
rm ${id}_rmptmt.sam
rm ${id}.sam
echo "findme"
done

mkdir Fastq4map
mv *.fq.gz Fastq4map
echo "finalbam done"





#3.标准化成bigwig，分别是RPKM延申，RPKM不延申，BPM延申、RPGC延申,顺手把bam变成bed做一下
ls -1 *.bam | awk -F "." '{print $1}' | while read id
do
echo "${id}"
bamCoverage -p 4 -b ${id}.bam -o ${id}_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
bamCoverage -e -p 4 -b ${id}.bam -o ${id}_RPKM_e.bigwig --binSize 10 --normalizeUsing RPKM
#bamCoverage -e 200 -p 4 -b ${id}.bam -o ${id}_RPKM_e200.bigwig --binSize 10 --normalizeUsing RPKM
echo "RPKM bigwig done"

bamCoverage -e -p 4 -b ${id}.bam -o ${id}_BPM.bigwig --binSize 10 --normalizeUsing BPM
echo "BPM_TPM bigwig done"
bamCoverage -e -p 4 -b ${id}.bam -o ${id}_RPGC.bigwig --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 119481543
echo "RPGC bigwig done"

bedtools bamtobed -i ${id}.bam > ${id}.bed
echo "bed done"
awk '{print "chr"$0}'  ${id}.bed > ${id}_chr.bed
echo "chr done"

rm ${id}.bed

done

mkdir Bam
mv *.bam Bam
mv *.bam.bai Bam

mkdir RPKMbigwig
mv *RPKM.bigwig RPKMbigwig
mkdir RPKMbigwig_e
mv *RPKM_e.bigwig RPKMbigwig_e
#mkdir RPKMbigwig_e200
#mv *RPKM_e200.bigwig RPKMbigwig_e200

#mv *.sam bamafter

mkdir Bed 
mv *.bed Bed
mv Bed Bam
mkdir BPMbigwig
mv *BPM.bigwig BPMbigwig
mkdir RPGCbigwig
mv *RPGC.bigwig RPGCbigwig

mkdir track
mv BPMbigwig track
mv RPGCbigwig track
mv RPKMbigwig track
mv RPKMbigwig_e track
#mv RPKMbigwig_e200 track
echo "bigwig done"


#4.看fragmentsize
cd ./Bam
ls *.bam >bam_name.txt
sed -i ':label;N;s/\n/ /;b label'  bam_name.txt
ls -1 *.bam | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup.bam" '{print $1}' >bam_name1.txt
sed -i ':label;N;s/\n/ /;b label'  bam_name1.txt
cat bam_name.txt |  while read id
do
echo "${id}"
cat bam_name1.txt |  while read slabel
do
echo "${slabel}"
bamPEFragmentSize  --bamfiles ${id} \
--histogram all_bam_fragsize.pdf  --numberOfProcessors 4 --maxFragmentLength 500   \
--samplesLabel ${slabel}
done
done
mkdir fragmentsize
mv *.pdf fragmentsize
echo "fragmentsize done" 
#rm bam_name.txt
cd ..


#5.构建matrix并且画图，分别是RPKM延申，RPKM不延申，BPM延申、RPGC延申
##################################
#################5.1RPKM_e###################################
cd ./track/RPKMbigwig_e
###############
ls *.bigwig >bw_name_e.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name_e.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPKM_e.bigwig" '{print $1}' >bw_name2.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name2.txt
###############
cat bw_name_e.txt |  while read id
do
echo "${id}"
cat bw_name2.txt |  while read slabel
do
echo "${slabel}"
###############
computeMatrix scale-regions \
-S ${id} \
-p 4 \
-R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
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
#rm bw_name_e.txt
#rm bw_name2.txt
mkdir allmatrix_e
mv *.pdf allmatrix_e
mv all_matrix_e.mat.gz allmatrix_e
cd ../..


###############
#################5.2RPKM###################################
cd ./track/RPKMbigwig
ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPKM.bigwig" '{print $1}' >bw_name3.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name3.txt
###############
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name3.txt |  while read slabel
do
echo "${slabel}"
###############
computeMatrix scale-regions \
-S ${id} \
-p 4 \
-R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
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
#rm bw_name.txt
#rm bw_name3.txt
mkdir allmatrix
mv *.pdf allmatrix
mv all_matrix.mat.gz allmatrix
cd ../..



###############
#################5.3BPM_e###################################
cd ./track/BPMbigwig
ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_BPM.bigwig" '{print $1}' >bw_name4.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name4.txt
###############
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name4.txt |  while read slabel
do
echo "${slabel}"
###############
computeMatrix scale-regions \
-S ${id} \
-p 4 \
-R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
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
#rm bw_name.txt
#rm bw_name4.txt
mkdir allmatrix
mv *.pdf allmatrix
mv all_matrix_e_BPM.mat.gz allmatrix
cd ../..







###############
#################5.4RPGC###################################
cd ./track/RPGCbigwig
ls *.bigwig >bw_name.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name.txt
ls -1 *.bigwig | awk -F "_rmptmt_q20_s1_fixmate_s2_markdup_RPGC.bigwig" '{print $1}' >bw_name5.txt
sed -i ':label;N;s/\n/ /;b label'  bw_name5.txt
###############
cat bw_name.txt |  while read id
do
echo "${id}"
cat bw_name5.txt |  while read slabel
do
echo "${slabel}"
###############
computeMatrix scale-regions \
-S ${id} \
-p 4 \
-R /mnt/data5/wjc/genome/ara/fixed_TAIR10/TAIR10_allgene_nochloroplast.bed \
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
#rm bw_name.txt
#rm bw_name5.txt
mkdir allmatrix
mv *.pdf allmatrix
mv all_matrix_e_RPGC.mat.gz allmatrix
cd ../..
echo "all kinds of patterns done"

#6.复制一些必要的代码，移动一些东西
mkdir code_for_more
cp /home/others/jccity/scripts/ara/macs2_callpeak_paired.sh code_for_more
#macs2进行callpeak
cp /home/others/jccity/scripts/ara/matix_peakcenter.sh code_for_more
#用center构建matrix
cp /home/others/jccity/scripts/ara/manorm.sh code_for_more
#差异peak分析
cp /home/others/jccity/scripts/ara/manorm_sort.r code_for_more
#筛选manorm的结果
cp /home/others/jccity/scripts/ara/GO.r code_for_more
#做GO富集
cp /home/others/jccity/scripts/ara/fragmentsize.sh code_for_more
#看片段大小
cp /home/others/jccity/scripts/ara/annotation_ara.r code_for_more
#将peaks和基因关联，得到peak基因
cp /home/others/jccity/scripts/ara/heatmap.rsh code_for_more
#画change的大热图
cp /home/others/jccity/scripts/ara/0_readme.txt code_for_more
#一个小说明
cp /home/others/jccity/scripts/ara/pattern.r code_for_more
#画pattern的r代码
echo "cp codes done" 


mv trimGalore_trim_s3_l20 Fastq4map
mkdir fq
mv Fastq4map fq
mv rawfastq fq
echo "mv fq done" 


#######还有很多代码我还没整理好，这个东西还在不断更新中，比如diffbind、sicer、bdgdiff等软件，用到再说吧
echo "all of chip done" 