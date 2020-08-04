#!/bin/bash


# ANALYSIS OF 12 SAMPLES SEQUENCED ON 1 FLOWCELL (DATA RECEIVED FROM FUNCTIONAL GENOMICS FACILITY)

# FOLDER CONTAINING FAST5 DIVIDED IN fast5_pass/ AND fast5_fail/
F5FILEDIR=/scicore/home/egliadr/GROUP/runQC/nanopore/Icarus/20190610_Icarus_13_24/data/AP20190605-GridION/p3077_5467_62/20190605_1515_GA20000_FAK22571_10c53a7a

# rearrange fast5 directory
cd $F5FILEDIR
mkdir -p fast5
mv fast5_pass/* fast5/
mv fast5_fail/* fast5/

# put each fast5 in a subfolder for easier basecalling
ls fast5/*fast5 | parallel -n1 mkdir fast5/{#}\;mv {} fast5/{#}

# SET VARIABLES
ID=run2
RAWDIR=$F5FILEDIR/fast5
CPUS=24

cd
mkdir -p Icarus2
cd Icarus2
mkdir -p $ID
cd $ID
mkdir -p slurms
mkdir -p logs
mkdir -p 01.basecall

# count number of raw folders
NRS=$(ls $RAWDIR | wc -l)
echo $NRS

# basecalling with guppy flip flop using GPU (faster!)
cat > slurms/all.gpu.basecall.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=$ID.all.basecalling
#SBATCH --gres=gpu:1
#SBATCH --partition=pascal
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=01.basecall/$ID.all.basecalling.log
#SBATCH --array=0-$NRS%1
~/ont-guppy/bin/guppy_basecaller -x "cuda:0" --kit SQK-LSK109 --flowcell FLO-MIN106 --gpu_runners_per_device 1 --disable_pings -s 01.basecall/\$SLURM_ARRAY_TASK_ID -i $RAWDIR/\$SLURM_ARRAY_TASK_ID/ 
EOF
sbatch slurms/all.gpu.basecall.pbs

# # concatenate all reads
# cat 01.basecall/*/*.fastq > all_reads.fastq
# cat 01.basecall.fast/*/*.fastq > all_reads.fast.fastq

# merge sequencing summary
module load R
R --vanilla <<RSCRIPT
library(data.table)
library(parallel)
fls <- list.files(path = "01.basecall", pattern = "summary.txt", full.names = T, recursive = T, include.dirs = T)
x <- rbindlist(mclapply(fls, fread, mc.cores = detectCores()))
fwrite(x, "merged_sequencing_summary.txt", sep="\t")
RSCRIPT

# run QC
MinIONQC.R --combined-only TRUE --smallfigures TRUE -p $CPUS -i merged_sequencing_summary.txt -o 03.QC

# # demux 
# conda activate icarus2
# qcat -t 24 --guppy --trim -o logs/$ID.demux.log -b qcatted -k NBD103/NBD104 -f all_reads.fastq
~/ont-guppy/bin/guppy_barcoder -t 24 -i 01.basecall/ --recursive -s 02.demux --barcode_kits EXP-NBD104 -q 0 

# trim and filter by quality
conda activate icarus # a conda environment with porechop, centrifuge, ktImportTaxonomy, kma, minimap2, racon
parallel -j 6 "
porechop -t 6 -i 02.demux/barcode{}/*.fastq --check_reads 1000 -o 02.demux/BC{}.trim.fastq --discard_middle
cat 02.demux/BC{}.trim.fastq | NanoFilt -q 7 | pigz -p4 > 02.demux/BC{}.clean.fastq.gz
rm 02.demux/BC{}.trim.fastq
" ::: {01..12}

# run centrifuge and remove human reads
mkdir -p 04.centrifuge
mkdir -p 05.bacterial
parallel -j12 "
centrifuge -x ~/tools/p_compressed+h+v -q -U 02.demux/BC{}.clean.fastq.gz -k 1 -p 2 --report-file 04.centrifuge/BC{}.centrifuge_report.tsv -S 04.centrifuge/BC{}.centrifuge_classification.tsv
cat 04.centrifuge/BC{}.centrifuge_classification.tsv | cut -f 1,3 > 04.centrifuge/BC{}.krona_classification.tsv
cat 04.centrifuge/BC{}.centrifuge_classification.tsv | awk '\$3 != \"9606\" {print \$1}' > 04.centrifuge/BC{}.nonhuman.reads.tsv
zcat 02.demux/BC{}.clean.fastq.gz | seqkit grep -f 04.centrifuge/BC{}.nonhuman.reads.tsv | pigz -p 4 > 05.bacterial/BC{}.fastq.gz
" ::: {01..12}

~/miniconda3/bin/ktImportTaxonomy 04.centrifuge/BC*.krona_classification.tsv -o 04.centrifuge/BCs.krona.html -tax ~/tools/taxonomy

# run KMA with latest NCBI AMR database
mkdir -p 06.amr
parallel "
kma -t 2 -i 05.bacterial/BC{}.fastq.gz -o 06.amr/BC{}.raw.kma -t_db ~/amrfinder/data/latest/AMR_CDS -mp 7 -mrs 0.0 -bcNano
minimap2 -t 2 -x map-ont ~/amrfinder/data/latest/AMR_CDS 05.bacterial/BC{}.fastq.gz | cut -f1 | uniq > 06.amr/BC{}.argreads.list 
zcat 05.bacterial/BC{}.fastq.gz | seqkit grep -j 4 -f 06.amr/BC{}.argreads.list  | pigz -p 2  > 06.amr/BC{}.argreads.fastq.gz
minimap2 -t 2 -x ava-ont 06.amr/BC{}.argreads.fastq.gz 06.amr/BC{}.argreads.fastq.gz | pigz -2 > 06.amr/BC{}.args.ava.paf.gz
racon -t 2 -f 06.amr/BC{}.argreads.fastq.gz 06.amr/BC{}.args.ava.paf.gz 06.amr/BC{}.argreads.fastq.gz > 06.amr/BC{}.polished.fasta
kma -t 2 -i 06.amr/BC{}.polished.fasta -o 06.amr/BC{}.polished.kma -t_db ~/amrfinder/data/latest/AMR_CDS -mp 7 -mrs 0.0 -bcNano
" ::: {01..12}

# assemble the metagenomes with flye
conda deactivate
conda activate nano.assembler # simpy a conda env with flye
mkdir -p 07.assembly
parallel -j4 "
flye --nano-raw 05.bacterial/BC{}.fastq.gz --out-dir 07.assembly/BC{} --genome-size 5m --threads 6 --meta
" ::: {01..12}

# run KMA on assemblies
parallel "
kma -t 2 -i 07.assembly/BC{}/assembly.fasta -o 06.amr/BC{}.assembly.kma -t_db ~/amrfinder/data/latest/AMR_CDS
" ::: {01..12}

# # count reads and bases
# for i in {01..12}
# do
# echo $i
# zcat 02.demux/BC$i.clean.fastq.gz | echo $((`wc -l`/4))
# zcat 02.demux/BC$i.clean.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c 
# done

# for i in {01..12}
# do
# echo $i
# zcat 05.bacterial/BC$i.fastq.gz | echo $((`wc -l`/4))
# zcat 05.bacterial/BC$i.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c 
# done



