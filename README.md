# RNA-seq DEG analysis
  - focusing on the analysing differential expression genes associated with the increasing of intramuscular fat deposition in Korean native cattle

## DEG Analysis pipeline
![workflow](RNA-seq_pipeline.png)

## 1. RNA library construction and NGS

## 2. Data pre-processing
  2-1. Importing data and download reference genome data
```
cd /[directory_path]/fastp_input
gunzip *.fq.gz
```

  2-2. Quality contorl
```
for infile in *_1.fq; do base=$(basename ${infile} _1.fq); echo ${base}_1.fq; done

for infile in *_1.fq; do base=$(basename ${infile} _1.fq); fastp -i ${base}_1.fq -I ${base}_2.fq -o /[directory_path]/${base}_trimmed_1.fastq -O /[directory_path]/${base}_trimmed_2.fastq --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 12 --html /[directory_path]/${base}_trimmed.fastp.html --json /[directory_path]/${base}_trimmed.fastp.json; done

mkdir fastq_results/

mv *.html /[directory_path]/fastp_results

mv *.json /[directory_path]/fastp_results

mkdir fastqc_results/

fastqc -o /[directory_path]/ *.fastq
```  
  - Viewing the FastQC results #check fastqc_html file out

  2-3. 

## 3-1. Exploratory data analysis


## 3-2. Expression analysis

