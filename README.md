# RNA-seq DEG analysis
  - focusing on the analysing differential expression genes associated with the increasing of intramuscular fat deposition in Korean native cattle

## DEG Analysis pipline

## 1. RNA library construction and NGS
  1-1. Importing data
```
cd /home/biolab302/jiyeon/fastp_input
gunzip *.fq.gz
```

  1-2. Quality contorl
```
for infile in *_1.fq; do base=$(basename ${infile} _1.fq); echo ${base}_1.fq; done

for infile in *_1.fq; do base=$(basename ${infile} _1.fq); fastp -i ${base}_1.fq -I ${base}_2.fq -o /home/biolab302/jiyeon/${base}_trimmed_1.fastq -O /home/biolab302/jiyeon/${base}_trimmed_2.fastq --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 12 --html /home/biolab302/jiyeon/${base}_trimmed.fastp.html --json /home/biolab302/jiyeon/${base}_trimmed.fastp.json; done

mkdir fastq_results/

mv *.html /home/biolab302/jiyeon/fastp_results

mv *.json /home/biolab302/jiyeon/fastp_results

mkdir fastqc_results/

fastqc -o /home/biolab302/jiyeon/fastqc_results/ *.fastq
```  
  - Viewing the FastQC results #check fastqc_html file out
  
## 2. Data pre-processing


## 3-1. Exploratory data analysis


## 3-2. Expression analysis

