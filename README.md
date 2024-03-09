# RNA-seq DEG analysis
  - focusing on the analysing differential expression genes associated with the increasing of intramuscular fat deposition in Korean native cattle

## DEG Analysis pipeline
![workflow](RNA-seq_pipeline.png)

![alt text](RNA-seq_DEG_pipeline.png)

## 1. RNA library construction and NGS

## 2. Data pre-processing
## 2-1. Importing data and download reference genome data
```
$ cd /[directory_path]/fastp_input
$ gunzip *.fq.gz
```

## 2-2. Quality contorl
```
$ for infile in *_1.fq; \
  do base=$(basename ${infile} _1.fq); \
  echo ${base}_1.fq; \
  done

$ for infile in *_1.fq; \
  do base=$(basename ${infile} _1.fq); \
  fastp \
  -i ${base}_1.fq \
  -I ${base}_2.fq \
  -o /[directory_path]/${base}_trimmed_1.fastq \
  -O /[directory_path]/${base}_trimmed_2.fastq \
  --detect_adapter_for_pe \
  --overrepresentation_analysis \
  --correction \
  --cut_right \
  --thread 12 \
  --html /[directory_path]/${base}_trimmed.fastp.html \
  --json /[directory_path]/${base}_trimmed.fastp.json; \
  done

$ mkdir fastq_results/
$ mv *.html /[directory_path]/fastp_results
$ mv *.json /[directory_path]/fastp_results
$ mkdir fastqc_results/
$ fastqc -o /[directory_path]/ *.fastq
```  
  - Viewing the FastQC results #check fastqc_html file out

## 2-3. 

## 3. Estimation of gene expression level
## 3-1. Building the STAR index
```
$ cd ncbi_dataset
$ STAR \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir /[directory_path] \
  --genomeFastaFiles GCF_002263795.2_ARS-UCD1.3_genomic.fna \
  --sjdbGTFfile genomic.gff \
  --sjdbGTFtagExonParentTranscript Parent \
  --sjdbOverhang 149 
```
## 3-2. Mapping trimmed_fastq file to the ARS-UCD1.3 bovine reference genome
```
$ cd 'trimmed data'
$ for infile in *_1.fastq; \
  do base=$(basename ${infile} _trimmed_1.fastq); \
  echo ${base}; \
  done 
$ for infile in *_1.fastq \
  ; do base=$(basename ${infile} _trimmed_1.fastq) \
  ; STAR \
    --runThreadN 12 \
    --genomeDir /[directory_path] \
    --sjdbGTFfile /[directory_path]/genomic.gff \
    --readFilesIn ${base}_trimmed_1.fastq ${base}_trimmed_2.fastq \
    --twopassMode Basic \
    --outFileNamePrefix /[directory_path]/${base}_ucd1.3_star_ \
    --outSAMtype BAM SortedByCoordinate \
  ; done

$ for infile in *_1.fastq; \
  do base=$(basename ${infile} _trimmed_1.fastq);
  STAR \
    --runThreadN 12 \
    --genomeDir /[directory_path] \
    --sjdbGTFfile /[directory_path]/genomic.gff \
    --readFilesIn ${base}_trimmed_1.fastq ${base}_trimmed_2.fastq \
    --twopassMode Basic \
    --outFileNamePrefix /[directory_path]/${base}_ucd1.3_star_ \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat; \
  done

$ samtools view \
  ucd1.3_star_mapping_L_5th_19_ucd1.3_star_Aligned.sortedByCoord.out.bam|more
```
- STAR align하기(파일 형식에 맞는 옵션 사용): fastq파일, #gzip파일
- samtools view: view first few alignment of BAM files

## 3-3. Re-buliding genome index using SJ.out_filtered.tab file
```
$ cd ucd1.3_star_mapping

$ for infile in *_SJ.out.tab; \
  do base=$(basename ${infile} _SJ.out.tab); \
  echo ${base}; \
  done

$ cat hw_L_5th_19_ucd1.3_star_SJ.out.tab|more

$ for infile in *_SJ.out.tab; \
  do base=$(basename ${infile} _SJ.out.tab); \
  cat ${base}_SJ.out.tab \
  | awk '($5 > 0 && $7 > 2 && $6==0)' \
  | cut -f1-6 \
  | sort \
  | uniq > ${base}_SJ.out_filtered.tab; \
  done

$ cat *.tab \
| awk '($5 > 0 && $7 > 2 && $6==0)' \
| cut -f1-6 \
| sort \
| uniq > hw_highlow_SJ.out_filtered.tab

$ cd [directory_path]

$ STAR \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir /[directory_path]/star_index_filtering_2pass \
  --genomeFastaFiles GCF_002263795.2_ARS-UCD1.3_genomic.fna \
  --sjdbGTFfile genomic.gff \
  --sjdbFileChrStartEnd /[direcotry_path]/hw_highlow_SJ.out_filtered.tab \
  --sjdbOverhang 149
```
- Build genome index with filtered output 

## 3-4. Mapping reads to the reference genome 2nd pass
```
$ cd 'trimmed data’

$ for infile in *_1.fastq; \
  do base=$(basename ${infile} _trimmed_1.fastq); \
  echo ${base}; \
  done

$ for infile in *_1.fastq; \
  do base=$(basename ${infile} _trimmed_1.fastq); \
  STAR \
  --runThreadN 12 \
  --readFilesIn ${base}_trimmed_1.fastq ${base}_trimmed_2.fastq \
  --genomeDir /[directory_path]/star_index_filtering_2pass/ \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /[directory_path]/${base}_star2pass_filtered \
  --outSAMunmapped Within \
  --readFilesCommand zcat; \
  done

$ for infile in *_1.fastq; \
  do base=$(basename ${infile} _trimmed_1.fastq); \
  STAR \
  --runThreadN 12 \
  --readFilesIn ${base}_trimmed_1.fastq ${base}_trimmed_2.fastq \
  --genomeDir /[directory_path]/star_index_filtering_2pass/ \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /[directory_path]/${base}_star2pass_filtered \
  --outSAMunmapped Within; \
  done

$ cat sample_5th_9_star2pass_filteredLog.final.out
```

## 3-5. Counting the read mapped to the reference genome sequence
  - Installation of htseq-count tool
```
$ pip install HTseq
```
  - Indexing bam file
```
$ cd [directory_paht]

$ for infile in *_filteredAligned.sortedByCoord.out.bam; \
  do base=$(basename ${infile} _filteredAligned.sortedByCoord.out.bam); \
  echo ${base}; \
  done

$ for infile in *_filteredAligned.sortedByCoord.out.bam; \
  do base=$(basename ${infile} _filteredAligned.sortedByCoord.out.bam); \
  samtools index ${base}_filteredAligned.sortedByCoord.out.bam; \
  done
```
  - Filtering genomic.gff file
```
$ cd /home/biolab302/바탕화면/jiyeon/ncbi_dataset/data/GCF_002263795.2

$ awk '/gene/' ./genomic.gff > genomic_gene.gff

$ cat genomic_gene.gff |more
```
  - mRNA quantification
```
$ cd /home/biolab302/바탕화면/jiyeon/analysis_file/ucd1.3_star_mapping

$ for infile in *_filteredAligned.sortedByCoord.out.bam; \
  do base=$(basename ${infile} _filteredAligned.sortedByCoord.out.bam); \
  echo ${base}; \
  done 

$ for infile in *_filteredAligned.sortedByCoord.out.bam; \
  do base=$(basename ${infile} _filteredAligned.sortedByCoord.out.bam); \
  htseq-count \
  --format bam \
  --order pos \
  --mode intersection-strict \
  --stranded reverse \
  --minaqual 1  \
  --type exon \
  --idattr gene \
  --add-chromosome-info ${base}_filteredAligned.sortedByCoord.out.bam /[directory_path]/genomic_gene.gff > ${base}_htseq-count.tsv; \
  done
```
## 4. Expression analysis of DEGs (in R)
