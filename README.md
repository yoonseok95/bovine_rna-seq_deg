# RNA-seq DEG analysis
  - focusing on the analysing differential expression genes associated with the increasing of intramuscular fat deposition in Korean native cattle

## DEG Analysis pipeline
![workflow](RNA-seq_pipeline.png)

![alt text](RNA-seq_DEG_pipeline.png)

## 1. RNA library construction and NGS
## 1-1. Sample
![alt text](sample.png)

## 2. Data pre-processing
## 2-1. Importing data and download reference genome data
```
$ cd /[directory_path]/fastp_input
$ gunzip *.fq.gz
```
```
add how to download bovine reference genome seq
```
  - Reference sequence
  -. UMD3.1 release 78    
     GenBank : GCA_000003055.5  
     RefSeq : GCF_000003055.6
     Name : Bos_taurus_UMD_3.1.1 (Data Nov 25, 2014) 
  -. bosTau9 (ARS_UCD1.2)
     GenBank : GCA_002263795.2
     RefSeq : GCF_002263795.1
     Name : bosTau9 (NCBI Annotation Release 106)
  -. bosTau9 (ARS_UCD1.3)
     GenBank : GCA_002263795.3
     RefSeq : GCF_002263795.2
     Name : NCBI eukaryotic genome annotation pipeline

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
```
Parameter		Description
--runThreadN		Number of threads (processors) for mapping reads to a genome
--runMode		run mode for STAR. genomeGenerate mode builds genome index
--genomeDir		PATH to the directory where genome indices will be stored
--genomeFastaFiles	reference genome file in FASTA format. You can also provide multiple genome files (separated by space) for index building.
--sjdbGTFfile		GTF file for gene annotation. This is an optional parameter
--sjdbOverhang		length of reads around the splice junctions. The ideal values should read length - 1 (or max read length - 1). 
                        For example, if your read length is 150, the value should be 149. In most cases, 
                        the default value of 100 also works.
```

## 3-2. Mapping trimmed_fastq file to the ARS-UCD1.3 bovine reference genome
## 3-2-1 Description: STAR(Spliced Transcripts Alignment to a Reference) Aligner
  -. Seed searching : For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):
  -.Seed1(mapped to the reference genome), Seed2(Searching again for only unmapped portion of the read)
  -. Clustering, stitching and scoring : scoring based on mismatches, indels, gaps, etc. 
## 3-2-2 Two-pass alignment of RNA-seq reads with STAR
  -.The 1-pass mapping mode generates all required data essential for many downstream analyses such as differential gene expression analysis  
  -.robustly and accurately identify novel splice junction for differential splicing analysis and variant discovery

### 1) Mapping
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
```
Parameter   description
--runThreadN		Number of threads (processors) for mapping reads to genome
--readFilesIn		Read files for mapping to the genome. In case of paired-end reads, provide read1 and read2 files. 
                        If there are multiple samples, separate files by a comma. For example, for paired-end reads, 
                        --readFilesIn S1read1.fastq,S2read1.fastq S1read2.fastq,S2read2.fastq
--genomeDir		PATH to the directory containing built genome indices
--outSAMtype 		'BAM SortedByCoordinate' Output coordinate sorted BAM file which is useful for many downstream analyses. This is optional.
--outSAMunmapped Within	Output unmapped reads from the main SAM file in SAM format. This is optional
--outFileNamePrefix	Provide output file prefix name
```
- STAR align하기(파일 형식에 맞는 옵션 사용): fastq파일, #gzip파일
- samtools view: view first few alignment of BAM files
- STAR output file: Parameter Description
```
seed_sampleAligned.sortedByCoord.out.bam	Alignment in BAM format (sorted by coordinate)
seed_sampleLog.final.out			Alignment summary statistics such as uniquely mapped reads, percent mapping, number of unmapped reads, etc.
seed_sampleLog.out				Alignment log for commands and parameters (useful in troubleshooting)
seed_sampleLog.progress.out			Alignment progress report (e.g. number of reads processed during particular span of time, mapped and
						unmapped reads, etc.)
seed_sampleSJ.out.tab				Filtered splice junctions found during the mapping stage
```
  - Splice junctions.
```
SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that
STAR defines the junction start/end as intronic bases, while many other software define them as
exonic bases. The columns have the following meaning:
column 1: chromosome
column 2: first base of the intron (1-based)
column 3: last base of the intron (1-based)
column 4: strand (0: undefined, 1: +, 2: -)
column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
AT/AC, 6: GT/AT
column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
column 7: number of uniquely mapping reads crossing the junction
column 8: number of multi-mapping reads crossing the junction
column 9: maximum spliced alignment overhang
The filtering for this output file is controlled by the --outSJfilter* parameters, as described in Section 14.16. Output Filtering: Splice Junctions.
```
```
$ iyoonseok95@yoonseok95:/mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star_mapping$ cat  hw_highlow_umd3.1.1SJ.out.tab|head
chromosome
AC_000158.1	125015	135819	2	2	1	6	0	71
AC_000158.1	136002	137157	2	2	1	12	0	66
AC_000158.1	137265	137833	2	2	1	9	0	66
AC_000158.1	137265	138874	2	2	0	1	0	13
AC_000158.1	137960	138874	2	2	1	5	0	75
AC_000158.1	138985	178432	2	2	1	3	0	70
AC_000158.1	242647	254558	1	1	0	2	0	41
AC_000158.1	242647	254564	1	1	0	1	0	59
AC_000158.1	242647	275622	1	1	0	9	0	58
AC_000158.1	242647	318512	1	1	0	2	0	56
```
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

![alt text](mapped_read.png)

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
  - See 'htseq-count_merge.R' file 

## 5. Normalization(R/FKM, TPM)
  - See 'Normalization_RPKM_TPM.R' file

## 5-1. Pattern between sample with high and low group
  - See 'expression_level_EDA.R' file
  - [best.sample.combination](https://bioinformatics-core-shared-training.github.io/Merged_RNASeq-course/html/02_Preprocessing_Data.nb.html)
  - Libarary size and distribution plots
  [best.sample.combination](https://bioinformatics-core-shared-training.github.io/Merged_RNASeq-course/html/02_Preprocessing_Data.nb.html)
  - To construct best sample combination in total reads through quality control
![Multidimensional scalliing](MDS.png)
  - QC of RNA-seq reads: distance between reads with high and low MS
![MA and volcano plot](volcano_plot.png)

## 6. DEG analysis
  - See 'DEG.R' file

## 7. GO(Gene ontology) and KEGG analysis
  - See 'GO_KEGG_visualization.R' file
![GO result](GO_sample.png)
![KEGG result](KEGG.png)


## Reference
[RNAseq_QC_HISAT2_DEseq2](https://www.nature.com/articles/s41597-022-01149-0)


