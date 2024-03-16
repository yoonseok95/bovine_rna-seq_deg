#!/bin/sh

#quality control of untrimmed fastq
mkdir fastq_results
mkdir trimmed_fastq
cd /mnt/sdb1/hanwoo_mRNA_DEG/untrimmed_mRNA_fastq

for infile in *_1.fq.gz 
do 
	base=$(basename ${infile} _1.fq.gz)
	fastp \
	-i ${base}_1.fq.gz \ 
	-I ${base}_2.fq.gz \
	-o /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq/${base}_trimmed_1.fq.gz \
	-O /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq/${base}_trimmed_2.fq.gz \
	--detect_adapter_for_pe \
	--overrepresentation_analysis \
	--correction \
	--cut_right \
	--thread 16 \
	--html /mnt/sdb1/hanwoo_mRNA_DEG/fastp_result/${base}_trimmed.fastp.html --json /mnt/sdb1/hanwoo_mRNA_DEG/fastp_result/${base}_trimmed.fastp.json 
done

mv *.html /hanwoo_mRNA_DEG/fastq_results
mv *.json /hanwoo_mRNA_DEG/fastq_results

#building index for mapping with STAR two-pass alignment
cd /mnt/sdb1/hanwoo_mRNA_DEG

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /home/iyoonseok95/b_taurus_ref/umd3.1.1/index_gff \
--genomeFastaFiles GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile genomic.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 149

#building index for mapping with STAR two-pass alignment
cd /mnt/sdb1/hanwoo_mRNA_DEG
mkdir umd3.1.1_star_2pass_mapping/

for infile in /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq/*_1.fq.gz
do 
	base=$(basename ${infile} _trimmed_1.fq.gz) 
	
	STAR --runThreadN 12 \
	--genomeDir /home/iyoonseok95/b_taurus_ref/umd3.1.1/index_gff \
	--sjdbGTFfile /home/iyoonseok95/b_taurus_ref/umd3.1.1/ncbi_dataset/data/GCF_000003055.6/genomic.gff \
	--readFilesIn /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq/${base}_trimmed_1.fq.gz /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq/${base}_trimmed_2.fq.gz \
	--twopassMode Basic \
	--outFileNamePrefix /mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star_2pass_mapping/${base}_umd3.1_star_ \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat
done

#Re-buliding genome index using SJ.out_filtered.tab file
cd /mnt/sdb1/hanwoo_mRNA_DEG
mkdir umd3.1.1_star_2pass_mapping/rebuilidng_filtering_tab
cd umd3.1.1_star_2pass_mapping

for infile in *_SJ.out.tab 
do 
	base=$(basename ${infile} _SJ.out.tab)
	
	cat ${base}_SJ.out.tab \
  	| awk '($5 > 0 && $7 > 2 && $6==0)' \
  	| cut -f1-6 \
	| sort \
	| uniq > rebuilidng_filtering_tab/${base}_SJ.out_filtered.tab
	
	cat rebuilidng_filtering_tab/*.tab \
	| awk '($5 > 0 && $7 > 2 && $6==0)' \
	| cut -f1-6 \
	| sort \
	| uniq > rebuilidng_filtering_tab/hw_highlow_SJ.out_filtered.tab	
done

#buliding genome index using filtered output
cd /mnt/sdb1/hanwoo_mRNA_DEG

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /home/iyoonseok95/b_taurus_ref/umd3.1.1/star_index_filtering_2pass240312 \
--genomeFastaFiles /home/iyoonseok95/b_taurus_ref/umd3.1.1/ncbi_dataset/data/GCF_000003055.6/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /home/iyoonseok95/b_taurus_ref/umd3.1.1/ncbi_dataset/data/GCF_000003055.6/genomic.gff \
--sjdbFileChrStartEnd umd3.1.1_star_2pass_mapping/rebuilidng_filtering_tab/hw_highlow_SJ.out_filtered.tab \
--sjdbOverhang 149

#2nd mapping(map paired-end reads to genome)
cd /mnt/sdb1/hanwoo_mRNA_DEG/trimmed_fastq

for infile in *_1.fq.gz
do 
	base=$(basename ${infile} _trimmed_1.fq.gz) 
	
	STAR --runThreadN 12 \
	--readFilesIn ${base}_trimmed_1.fq.gz ${base}_trimmed_2.fq.gz \
	--genomeDir /home/iyoonseok95/b_taurus_ref/umd3.1.1/star_index_filtering_2pass240312 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix /mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star_2pass_mapping/${base}_umd3.1_star_2nd_filtered \
	--outSAMunmapped Within \
	--readFilesCommand zcat	
done

#indexing bam file with HTseq
cd /mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star_2pass_mapping

for infile in *_umd3.1_star_2nd_filteredAligned.sortedByCoord.out.bam
do 
	base=$(basename ${infile})
	
	samtools index ${base}
done

#filtering genomic_gene.gff file
#cd /mnt/sdb1/hanwoo_mRNA_DEG
#awk '/gene/' ./genomic.gff > genomic_gene.gff

#mRNA quantification
cd /mnt/sdb1/hanwoo_mRNA_DEG/
mkdir umd3.1.1_star2pass_htseq_count/
cd umd3.1.1_star_2pass_mapping

for infile in *_umd3.1_star_2nd_filteredAligned.sortedByCoord.out.bam
do 
	base=$(basename ${infile} _umd3.1_star_2nd_filteredAligned.sortedByCoord.out.bam)
	
	htseq-count \
	--format bam \
	--order pos \
	--mode intersection-strict \
	--stranded reverse \
	--minaqual 1  \
	--type exon \
	--idattr gene \
	--add-chromosome-info ${base}_umd3.1_star_2nd_filteredAligned.sortedByCoord.out.bam ./genomic_genev2.gff > /mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star2pass_htseq_count/${base}_umd3.1.1_star2pass_htseq-count.tsv
done



