#!/bin/bash

source activate micro

cd
cd icarus

WT=612281-16_S1
MUT=612279-15_S2
NCPU=16

mkdir trimmed
# TRIM
trimmomatic PE -threads $NCPU -phred33 \
ILM.raw/"$WT"_L001_R1_001.fastq.gz ILM.raw/"$WT"_L001_R1_001.fastq.gz trimmed/trim."$WT".R1.fastq.gz trimmed/trim.discarded."$WT".R1.fastq.gz \
trimmed/trim."$WT".R2.fastq.gz trimmed/trim.discarded."$WT".R2.fastq.gz \
ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 MINLEN:100

rm trimmed/trim.discarded."$WT".R1.fastq.gz trimmed/trim.discarded."$WT".R2.fastq.gz

trimmomatic PE -threads $NCPU -phred33 \
ILM.raw/"$MUT"_L001_R1_001.fastq.gz ILM.raw/"$MUT"_L001_R1_001.fastq.gz trimmed/trim."$MUT".R1.fastq.gz trimmed/trim.discarded."$MUT".R1.fastq.gz \
trimmed/trim."$MUT".R2.fastq.gz trimmed/trim.discarded."$MUT".R2.fastq.gz \
ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 MINLEN:100

rm trimmed/trim.discarded."$MUT".R1.fastq.gz trimmed/trim.discarded."$MUT".R2.fastq.gz


# assemble the genome
mkdir -p unicycler.out
unicycler --threads $NCPU --short1 trimmed/trim."$WT".R1.fastq.gz  --short2 trimmed/trim."$WT".R2.fastq.gz --out unicycler.out/$WT


mkdir -p prokka.out

prokka \
--kingdom Bacteria \
--outdir prokka.out/$WT \
--prefix WT.psaer \
--genus Pseudomonas \
--locustag WT.psaer unicycler.out/$WT/assembly.fasta \
--cpus $NCPU \
--force


# use snippy
snippy --cpus $NCPU --outdir snippy.out --ref prokka.out/612281-16_S1/WT.psaer.gbk --R1 trimmed/trim."$MUT".R1.fastq.gz --R2 trimmed/trim."$MUT".R2.fastq.gz


# snippy correctly identifies the two mutations in oprD and ampR genes, now I estract the two genes for following nanopore reads mapping
cat > extract.genes << EOF
WT.psaer_01028
WT.psaer_04014
EOF

seqtk subseq prokka.out/612281-16_S1/WT.psaer.ffn extract.genes > oprD_ampR_WT.fasta


# map MUT only against these two genes (only for visualization)
mkdir -p bams
bwa index oprD_ampR_WT.fasta
bwa mem -t $NCPU oprD_ampR_WT.fasta trimmed/trim."$MUT".R1.fastq.gz trimmed/trim."$MUT".R2.fastq.gz | samtools sort -o bams/"$MUT".OprdAmpr.bam -@$NCPU
samtools index bams/"$MUT".OprdAmpr.bam


##########
# NOW WITH NANOPORE READS
# R4BC11 seems to be the mutated one
mkdir ONT.map

# fastqs are located in "ONT.raw/" folder
ONTMUT=R4BC11 # MUT sample ID
ONTWT=R4BC12 # WT sample ID
F5FOLDER=/data/admin/icarus/extracted.fast5/ # folder with Nanopore fast5 files

minimap2 -t $NCPU -x ava-ont ONT.raw/$ONTMUT.fastq.gz ONT.raw/$ONTMUT.fastq.gz > ONT.map/$ONTMUT.reads.paf
racon -t $NCPU -f ONT.raw/$ONTMUT.fastq.gz ONT.map/$ONTMUT.reads.paf ONT.raw/$ONTMUT.fastq.gz > ONT.map/polished.$ONTMUT.reads.fasta
minimap2 -ax map-ont -t $NCPU oprD_ampR_WT.fasta ONT.map/polished.$ONTMUT.reads.fasta | samtools sort -o ONT.map/$ONTMUT.OprdAmpr.bam -@ $NCPU
samtools index -@ $NCPU ONT.map/$ONTMUT.OprdAmpr.bam

minimap2 -t $NCPU -x ava-ont ONT.raw/$ONTWT.fastq.gz ONT.raw/$ONTWT.fastq.gz > ONT.map/$ONTWT.reads.paf
racon -t $NCPU -f ONT.raw/$ONTWT.fastq.gz ONT.map/$ONTWT.reads.paf ONT.raw/$ONTWT.fastq.gz > ONT.map/polished.$ONTWT.reads.fasta
minimap2 -ax map-ont -t $NCPU oprD_ampR_WT.fasta ONT.map/polished.$ONTWT.reads.fasta | samtools sort -o ONT.map/$ONTWT.OprdAmpr.bam -@ $NCPU
samtools index -@ $NCPU ONT.map/$ONTWT.OprdAmpr.bam

# assemble WT genome
# assemble the wt with flye and annotate with prokka
mkdir -p flye.out
flye --nano-raw ONT.raw/$ONTWT.fastq.gz --out-dir flye.out/$ONTWT --genome-size 6.5m --threads $NCPU

prokka \
--kingdom Bacteria \
--outdir prokka.out/$ONTWT \
--prefix WT.psaer.nanopore \
--genus Pseudomonas \
--locustag WT.psaer.nanopore flye.out/$ONTWT/assembly.fasta \
--cpus $NCPU \
--force



# NANOPOLISH
nanopolish index -d $F5FOLDER ONT.raw/$ONTMUT.fastq.gz
nanopolish index -d $F5FOLDER ONT.raw/$ONTWT.fastq.gz
mkdir -p nanopolish.out

minimap2 -ax map-ont -t $NCPU oprD_ampR_WT.fasta ONT.raw/$ONTMUT.fastq.gz | samtools sort -o ONT.map/$ONTMUT.raw.OprdAmpr.bam -@ $NCPU
samtools index -@ $NCPU ONT.map/$ONTMUT.raw.OprdAmpr.bam


# CALL MUTATIONS ONLY FOR THE TWO GENES
nanopolish variants \
-t $NCPU \
--reads ONT.raw/$ONTMUT.fastq.gz \
--bam ONT.map/$ONTMUT.raw.OprdAmpr.bam \
--genome oprD_ampR_WT.fasta \
--ploidy 1 \
-d 10 \
--snps \
-w "WT.psaer_01028:1-1332" \
-o nanopolish.out/$ONTMUT.polished.oprD.vcf

nanopolish variants \
-t $NCPU \
--reads ONT.raw/$ONTMUT.fastq.gz \
--bam ONT.map/$ONTMUT.raw.OprdAmpr.bam \
--genome oprD_ampR_WT.fasta \
--ploidy 1 \
-d 10 \
--snps \
-w "WT.psaer_04014:1-891" \
-o nanopolish.out/$ONTMUT.polished.ampR.vcf


minimap2 -ax map-ont -t $NCPU oprD_ampR_WT.fasta ONT.raw/$ONTWT.fastq.gz | samtools sort -o ONT.map/$ONTWT.raw.OprdAmpr.bam -@ $NCPU
samtools index -@ $NCPU ONT.map/$ONTWT.raw.OprdAmpr.bam

nanopolish variants \
-t $NCPU \
--reads ONT.raw/$ONTWT.fastq.gz \
--bam ONT.map/$ONTWT.raw.OprdAmpr.bam \
--genome oprD_ampR_WT.fasta \
--ploidy 1 \
-d 10 \
--snps \
-w "WT.psaer_01028:1-1332" \
-o nanopolish.out/$ONTWT.polished.oprD.vcf

nanopolish variants \
-t $NCPU \
--reads ONT.raw/$ONTWT.fastq.gz \
--bam ONT.map/$ONTWT.raw.OprdAmpr.bam \
--genome oprD_ampR_WT.fasta \
--ploidy 1 \
-d 10 \
--snps \
-w "WT.psaer_04014:1-891" \
-o nanopolish.out/$ONTWT.polished.ampR.vcf


# # call variants full assembly
# nanopolish_makerange.py flye.out/$ONTWT/assembly.fasta > ranges.nanopolish
# python3 nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
#     nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1


# minimap2 -ax map-ont -t $NCPU flye.out/$ONTWT/assembly.fasta ONT.raw/$ONTMUT.fastq.gz | samtools sort -o ONT.map/$ONTMUT.raw.assembly.bam -@ $NCPU
# samtools index -@ $NCPU ONT.map/$ONTMUT.raw.assembly.bam

# parallel "
# nanopolish variants \
# -t 1 \
# --reads ONT.raw/$ONTMUT.fastq.gz \
# --bam ONT.map/$ONTMUT.raw.assembly.bam \
# --genome flye.out/$ONTWT/assembly.fasta \
# --ploidy 1 \
# -d 20 \
# --snps \
# -w {} \
# -o nanopolish.out/$ONTMUT.{}.vcf
# " :::: ranges.nanopolish

# bcftools cat nanopolish.out/$ONTMUT.contig* > nanopolish.out/$ONTMUT.fullgenome.vcf

