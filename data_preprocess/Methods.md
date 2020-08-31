# Data preprocess pipeline

[TOC]
## Exome Analysis pipeline

### Samples Features

| Samples                                         | Kit                                          | Target Size |
| ----------------------------------------------- | -------------------------------------------- | ----------- |
| Three local matched trios                       | Illumina TrueSeq Exome kit                   | 45MB        |
| Siventeen local matched trios                   | Agilent SureSelect Human All Exon V6         | 58MB        |
| Twenty matched trios of Oncotarget 2016         | SeqCap EZ Exome v2.0 from Roche Life Science | 44MB        |
| Sixty nine matched trios of Genome Biology 2014 | MSK-IMPACT 230 Gene Panel                    |             |
| Matched trios of Cancer Cell 2018               | MSK-IMPACT 341 or 410 Gene Panel             |             |

### Prepare the target region for exome seq

Sort and merge bed files:

```bash
bed="./S07604624_Padded.bed"  ## Original target region file
sorted_bed="./V6_padded.sorted.bed"  
merged_bed="./V6_padded.merged.bed"  ## Output for next step

sort -k1,1 -k2,2n $bed > $sorted_bed
# The `-sorted` option in `bedtools` can let the program run faster
# but need the sort the bed file before above step
bedtools merge -i $sorted_bed > $merged_bed
```

Merge two bed files:

```bash
bedtools intersect -a target_regionA.bed -b target_regionB.bed -sorted > intersect_output.bed
```

At last, the intersection of _illumina_, _agilent_ and _rocha_ exome capture kit is used.

Before using the intersected target region bed file, check the chromosomes: if the chromosome are in the reference sequnce  used for mapping.

```perl
#!/usr/bin/env perl6

## Get genome list list
## Read bed 
## Sort by fai list

my @fai = qx[awk '{print $1}' ./GATK_bundle/genome_list.txt].words;
my $bed = "./intersect_illumina_agilen_roche_merge.bed";
my $output = 
	"./intersect_illumina_agilen_roche_merge_sortByFai.bed";

qqx[rm $output];

for @fai -> $chr {
    say $chr;
    qqx[grep "$chr\t" $bed >> $output];
}
```

 

## Local sample Analysis

### Mapping & Variation Calling

1. Raw fastq to uBam

   ```bash
   ## fastq to sam
   ${javaCmd} -Xmx${memory}G -jar ${picardJar} FastqToSam \
     FASTQ=${inputFastq1} \
     FASTQ2=${inputFastq2} \
     OUTPUT=${uBam} \
     READ_GROUP_NAME=${sample_name} \
     SAMPLE_NAME=${sample_name} \
     LIBRARY_NAME=${sample_name} \
     PLATFORM_UNIT=illumina \
     PLATFORM=illumina \
     TMP_DIR=/tmp
   ```

   

2. Mark adapter

   ```bash
   ## mark adapter
   ${javaCmd} -Xmx${memory}G -jar ${picardJar} \
     MarkIlluminaAdapters \
     I=${uBam} \
     O=${adaptersMarkedBam} \
     M=${metricsFile} \
     TMP_DIR=/tmp
   ```

   

3. Mapping, Sort and index

   ```bash
   ${javaCmd} -Xmx${memory}G -jar ${picardJar} SamToFastq \
     I=${adaptersMarkedBam} \
     FASTQ=/dev/stdout \
     CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 \
     INTERLEAVE=true NON_PF=true \
     TMP_DIR=/tmp | \
     bwa mem -M -t ${threadBwa} -p ${ref_fasta} /dev/stdin | \
     ${javaCmd} -Xmx${memory}G -jar ${picardJar} \
     MergeBamAlignment \
     ALIGNED_BAM=/dev/stdin \
     UNMAPPED_BAM=${uBam} \
     OUTPUT=${mappedBam} \
     R=${ref_fasta} CREATE_INDEX=true ADD_MATE_CIGAR=true \
     CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
     INCLUDE_SECONDARY_ALIGNMENTS=true \
     MAX_INSERTIONS_OR_DELETIONS=-1 \
     PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
     ATTRIBUTES_TO_RETAIN=XS \
     TMP_DIR=/tmp
   ```

   

4. Indel Realignment

   ```bash
   ${javaCmd} -Xmx${memory}G -Djava.io.tmpdir=${tmpDir} -jar ${gatkJar} -T RealignerTargetCreator \
     -nt ${nt} \
     -R ${ref_fasta} \
     -I ${mappedBam} \
     -known ${know_indels} \
     -o ${realignerIntervals}

   ${javaCmd} -Xmx${memory}G -Djava.io.tmpdir=${tmpDir} \
   -jar ${gatkJar} -T IndelRealigner \
     -R ${ref_fasta} \
     -I ${mappedBam} \
     -known ${know_indels} \
     -targetIntervals ${realignerIntervals} \
     -o ${realignedBam}
   ```

   

5. Basepair Calibration

   ```bash
   ${javaCmd} -Xmx${memory}G -Djava.io.tmpdir=${tmpDir} -jar ${gatkJar} -T BaseRecalibrator\
     -nct ${nct} \
     -R ${ref_fasta} \
     -I ${mappedBam} \
     -knownSites ${dbsnp138} \
     -knownSites ${know_indels} \
     -o ${recalibTable}

   ${javaCmd} -Xmx${memory}G -Djava.io.tmpdir=${tmpDir} -jar ${gatkJar} -T PrintReads\
     -nct ${nct} \
     -R ${ref_fasta} \
     -I ${realignedBam} \
     -BQSR ${recalibTable} \
     -o ${recalibBam}
   ```

   

6. Remove duplication

   ```bash
   java -jar ~/bin/picard.jar \
   MarkDuplicatesWithMateCigar \
         I=$recalBam \
         O=$rmdup \
         M=${rmdup}_mark_dup_matrix.txt \
         REMOVE_DUPLICATES=true
   ```

   

7. UnifiedGenotyper

   ```bash
   ## Variant call
   java -Xmx70G -jar $gatk_jar -T UnifiedGenotyper \
     -nt 10 \
     -R $refseq \
     -I $rmdupN \  ## Normal
     -I $rmdupC \  ## Primary Colorectal cancer
     -I $rmdupM \  ## Liver Metastasis Colorectal cancer
     -o $vcf_raw \
     -glm BOTH \
     --min_base_quality_score 15 \
     --min_indel_count_for_genotyping 2 \
     --min_indel_fraction_per_sample 0.05 \
     --output_mode EMIT_VARIANTS_ONLY \
     --sample_ploidy 2 \
     --standard_min_confidence_threshold_for_calling 10.0 \
     --standard_min_confidence_threshold_for_emitting 10.0 

   vcftools --vcf ${vcf_raw} \
     --min-meanDP 10 \
     --out ${vcf_deepth} \
     --recode Annotations
   ```

   

8. Annotation

   ```bash
   ## target resgion selection
   grep "^#" $vcf > $vcf_exon
   bedtools intersect -a $vcf -b $target_bed >> $vcf_exon

   ## convert vcf to annovar format
   convert2annovar.pl -includeinfo \
     -allsample \
     -withfreq \
     -format vcf4 \
     -outfile ${annovar_input} \
     ${vcf_exon}
   ```


9. Correct Allele Deepth 

   ```bash
   ## RUN
   ## INPUT    1:CHROM 2:LOC 4:REF 5:ALT 18:C 19:M 20:N
   ## OUTPUT   1:CHROM 2:LOC 4:REF 5:ALT 18:C.AD 19:M.AD 20:N.AD
   ## Annotation: C = primary cancer, M = Metastasis cancer, N = Normal tissue
   awk '{
   split($18, array18, /[:]/);
   split($19, array19, /[:]/);
   split($20, array20, /[:]/);
   print \
   $1"\t"$2"\t"$4"\t"$5"\t"array18[2]"\t"array19[2]"\t"array20[2];
   }' $annovar_input | ./lib/multi_allele_depth_correct.pl6 > \
   $AD_correct_genotype_field
   ```



### Annotation gene

```bash
table_annovar.pl --csvout \
    -buildver hg19 \
    -out ${annovar_gene_prefix} \
    -remove \
    -protocol refGene,cytoBand,genomicSuperDups,bed \
    -operation g,r,r,r \
    -bedfile hg19_cosmic70_region.txt \
    -nastring . \
    ${annovar_input} ${humandb}
```


### Remove intron

```bash
Rscript ./lib/rm_FuncInRefGene.R ${rm_intron_input} ${rm_intron_annovar_input} "intergenic" "intronic" "ncRNA_intronic"
```



### Annotation location

```bash
table_annovar.pl --csvout \
    -buildver hg19 \
    -out ${annovar_local_prefix} \
    -remove \
    -protocol esp6500siv2_all,1000g2015aug_all,exac03,cosmic70,dbnsfp33a,clinvar_20170130 \
    -operation f,f,f,f,f,f \
    -nastring . \
    ${rm_intron_annovar_input} ${humandb}
```



## SRA sample Analysis

### Download SRA files

0. Install aspera transfer software.
1. Find the SRA accession number in SRA database.
2. Goto SRA ftp site:
    https://www.ncbi.nlm.nih.gov/public/?/ftp/sra/sra-instant/reads/
    Choose the right SRA dataset and download
3. Install SRA-toolkit:
    https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std
4. Use SRA-toolkit to transform the SRA file into fastq format.


### Transform SRA files to Fastq

```bash
##  sra to fastq
fastq-dump -X 100000 --skip-technical --gzip \
--split-files $SRR_file \
-O ./data/fastq/
```

__The following preparation of SRA samples follows the the same pipeline of the local exome sequence samples__

## Genome Biology dataset Analysis

This datasets does not provide the original sequencing fastq data. It only provides the mutation sites and Variant Allele frequency. Only need to transform its format and re-annotation.

Transform format:

```R
mutTab <- fread("data_mutations_extended.txt", sep = "\t", skip=1)

## Annovar input
# OUTPUT FORMAT(NO HEADER): CHROM	POS REF ALT
outDir <- "../tmp/02_annovar/"
dir.create(outDir, showWarning=F, recursive=T)

tab <- mutTab[, .(
	CHROM        = paste0("chr", Chromosome)
	, START      = Start_Position
	, POS        = Start_Position
	, END        = End_Position
	, REF        = Reference_Allele
	, ALT        = Tumor_Seq_Allele2
	, personID   = Tumor_Sample_Barcode %>% sub("-.$", "", .)
	, T.AD       = paste0(t_ref_count, ",", t_alt_count)
	, N.AD       = paste0(n_ref_count, ",", n_alt_count)
	, sampleType = Tumor_Sample_Barcode %>% sub(".*-", "", .)
)] %>% unique

tabC <- tab[sampleType == "P", .(CHROM, START, END, POS, REF, ALT, personID, C.AD=T.AD, N.AD)]
tabM <- tab[sampleType == "M", .(CHROM, START, END, POS, REF, ALT, personID, M.AD=T.AD, N.AD)]
setkeyv(tabC, c("CHROM", "POS", "REF", "ALT", "personID"))
setkeyv(tabM, c("CHROM", "POS", "REF", "ALT", "personID"))

tabCM <- rbind(tabC[tabM], tabM[tabC]) %>% unique
tabCM[is.na(C.AD), C.AD := "0,0"]
tabCM[is.na(M.AD), M.AD := "0,0"]
tabCM[is.na(N.AD), N.AD := i.N.AD]
tabCM[is.na(START), START := i.START]
tabCM[is.na(END), END := i.END]

persons <- tab$personID %>% unique
for (i in persons) {
	file <- paste0(outDir, i, "_annovar", ".input")
	write_tsv(tabCM[personID == i, .(CHROM, START, END, REF, ALT)], file, col_names=F)
}

## Allele frequency
# OUTPUT FORMAT: CHROM	POS	REF	ALT	C.AD	M.AD	N.AD
outDir <- "../tmp/03_genotype_field/"
dir.create(outDir, showWarning=F, recursive=T)

for (i in persons) {
	file <- paste0(outDir, i, "_genotype_field_AD_correct.tsv")
	write_tsv(tabCM[personID == i, .(CHROM, POS, REF, ALT, C.AD, M.AD, N.AD)] , file)
}
```

Annotation:

```bash
#!/bin/bash

# INPUT: annotation.sh sample_mark humandb prefix


set -o nounset
#################################################
# Parameters
#################################################
sampleMark=$1
humandb=$2
prefix=$3


#################################################
## Setting parameter
#################################################

mkdir -p $prefix

annovar_input=$prefix/${sampleMark}_annovar.input
annovar_gene_prefix=$prefix/${sampleMark}_gene
rm_intron_input=${annovar_gene_prefix}.hg19_multianno.csv
rm_intron_annovar_input=$prefix/${sampleMark}_rm_intron_annovar.input
annovar_local_prefix=$prefix/${sampleMark}_local


#################################################
## Run
## annotation gene 
table_annovar.pl --csvout \
	-buildver hg19 \
	-out ${annovar_gene_prefix} \
	-remove \
	-protocol refGene,cytoBand,genomicSuperDups,bed \
	-operation g,r,r,r \
	-bedfile hg19_cosmic70_region.txt \
	-nastring . \
	${annovar_input} ${humandb}


## extract remove intron
rm_FuncInRefGene.R ${rm_intron_input} ${rm_intron_annovar_input} "intergenic" "intronic" "ncRNA_intronic"

## annotation location
table_annovar.pl --csvout \
	-buildver hg19 \
	-out ${annovar_local_prefix} \
	-remove \
	-protocol esp6500siv2_all,1000g2015aug_all,exac03,cosmic70,dbnsfp33a,clinvar_20170130 \
	-operation f,f,f,f,f,f \
	-nastring . \
	${rm_intron_annovar_input} ${humandb}

```



## Cancer Cell dataset Analysis

Most sample of Cancer Cell 2018 is unpaired. Only 21 patients are paired, and 3 patients have more than 3 sample per person.

_Only use the matched trios sample._

Annotation: Using the same script to Genome Biology 2014.



## Multiregion Biopsis samples

There are mapped bam files available in SRA database.

From the `samtools view -H SRR.bam` output , we can know that the duplicats have been removed. The file is also _Sorted_ and _Indexed_. 

```
@PG     ID:GATK IndelRealigner  VN:2.1-13-g1706365      CL:knownAlleles=[] targetIntervals=CHET9_1D_ATCACG_L007_dupRemoved.intervals LODThresho
ldForCleaning=5.0 consensusDeterminationModel=USE_READS entropyThreshold=0.15 maxReadsInMemory=150000 maxIsizeForMovement=3000 maxPositionalMov
eAllowed=200 maxConsensuses=30 maxReadsForConsensuses=120 maxReadsForRealignment=20000 noOriginalAlignmentTags=false nWayOut=null generate_nWay
Out_md5s=false check_early=false noPGTag=false keepPGTags=false indelsFileForDebugging=null statisticsFileForDebugging=null SNPsFileForDebuggin
g=null
@PG     ID:MarkDuplicates       PN:MarkDuplicates       VN:1.85(1345)   CL:net.sf.picard.sam.MarkDuplicates INPUT=[CHET9_1D_ATCACG_L007_RG.bam]
 OUTPUT=CHET9_1D_ATCACG_L007_dupRemoved.bam METRICS_FILE=CHET9_1D_ATCACG_L007_dupRemoved.met REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILEN
T    PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates ASSUME_SORTED=false MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE
_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPL
ICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
```

So we only need to run GATK pipeline which is described in the local sample process part:

1. Indel Realignment
2. Basepair Calibration
3. UnifiedGenotyper
4. Annotation
5. Correct Allele Deepth

## Muation filtering

### 00 QC

### 01 Mutation Filtering

#### Filter by VAF

cutpointT: tumor mutation VAF cutpoint is 5%, larger than that is considered somatic mutation.

cutpointN: normal tissue contain a mutation VAF larger than 1%, which will be considered as germline mutation.

The tumor variant allele depth should be > 3 and the depth of the normal tissue shoud be > 5.

#### Filter by annotation

```
(!(ExonicFunc.refGene %in% c("synonymous SNV")))
& (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
& (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
& (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
& (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
& ((ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "stopgain"))
        | (SIFT_pred == "D" | Polyphen2_HDIV_pred == "D")
        | (Func.refGene %in% "splicing")
  )
```

Mutation is not synomous mutation.

**AND**

Mutation should be one of the types of "exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic" and "splicing".

**AND**

The mutation rate of site in esp6500siv2all datasets less than 0.01 or none. 

**AND**

The mutation rate of the site in 1000g2015aug_all datasetes is less than 0.01 or none.

**AND**

The mutation rate of the site in ExAC_ALL datasets is less that 0.05 or none.

**AND**

The mutation type is  one of "frameshift deletion", "frameshift insertion" or stopgain" 

​    **OR** 

​    The site effect predict SIFT_pred is "D" or Polyphen2_HDIV_pred is "D"

​    **OR**

​    The mutation site is splicing site.

###  02 More filter

cutpointT: tumor mutation VAF cutpoint is 10%, larger than that is considered somatic mutation.

cutpointN: normal tissue contain a mutation VAF larger than 1%, which will be considered as germline mutation.

### 03 Then the data are really to further analysis







