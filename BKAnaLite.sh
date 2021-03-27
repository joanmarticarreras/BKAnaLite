#!/usr/bin/env bash

# Author:	Joan Martí-Carreras
# E-mail:	joan.marti.carreras@gmail.com
# Info: 	www.joanmarticarreras.com
# BKAnaLite is an automated pipeline meant to assemble polyoma BK virus genomes and performe viral-aware variant calling from Illumina reads
# Released under GPL3

set -e

# OPTIONS
## Variable initialization
r=""
R1=""
R2=""
THREADS=0

## Missing options execption case
if [ $# -eq 0 ]
then
        echo "Missing options!"
        echo "Run $0 -h for help"
        echo ""
        exit 0
fi

OPTIND=1
# Resetting OPTIND is necessary if getopts was used previously in the script.
# It is a good idea to make OPTIND local if you process options in a function.

## Option assignation
while getopts r:1:2:s:t:x:o:h OPTIONS
do
	case ${OPTIONS} in

                r)
			r=${OPTARG};;
                1)
			R1=${OPTARG};;
                2)
			R2=${OPTARG};;
		s)
			S=${OPTARG};;

                t)
			THREADS=${OPTARG};;
                x)
			BASENAME=${OPTARG};;
                o)	OUTPUT_DIR=${OPTARG};;
		h)
                        echo "Usage: BKAnaLite.sh [OPTIONS] "
                        echo ""
                        echo "   -r     Path to reference genome in FASTA format"
                        echo "   -1     Path to R1 Illunmina PE reads"
                        echo "   -2     Path to R2 Illumina PE reads"
			echo "   -s     Path to Sanger reads in FASTQ format"
                        echo "   -t     Number of threads to use when multithreading is possible"
			echo "   -x     Basename for the output files"
			echo "   -o     Path to output directory"
                        echo "   -h     help (this output)"
                        exit 0
                        ;;

        esac
done

shift "$((OPTIND-1))"   # Discard the options and sentinel --

# INITIALIZE SAMPLE
## Create directory structure
mkdir -p $OUTPUT_DIR/$BASENAME/{raw,trim,ref_assembly/{map,lofreq,annotation},novo_assembly/{pilon,map,lofreq,annotation},comparison_reference/lofreq,results}
## Copy raw reads
cp $R1 $OUTPUT_DIR/$BASENAME/raw/$BASENAME.R1.fastq.gz
cp $R2 $OUTPUT_DIR/$BASENAME/raw/$BASENAME.R2.fastq.gz
## Copy reference
cp ${r} $OUTPUT_DIR/$BASENAME/reference.fa

# TRIMMING Illumina reads for ASSEMBLY  Q5 + Nextera adapters + minlength 30
bbduk.sh in=$R1 in2=$R2 \
	ref=<PATH-TO-ADAPTER-REF>/BKAnaLite/references/Nextera-polyoma.fasta \
	out=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R1.fastq.gz \
	out2=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R2.fastq.gz \
	threads=$THREADS minlength=150 trimq=30 qtrim=w minavgquality=20 overwrite=true # minlength=30

bbnorm.sh in=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R1.fastq.gz in2=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R2.fastq.gz  \
	out=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R1.fastq.gz out2=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R2.fastq.gz \
	threads=$THREADS target=99999999 min=500 passes=1

bbnorm.sh in=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R1.fastq.gz in2=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R2.fastq.gz \
	out=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R1.fastq.gz out2=$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R2.fastq.gz \
	threads=$THREADS target=1000 ecc=t keepall passes=1 bits=16 prefilter

# ASSEMBLY
## De novo assembly
iva -f $OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R1.fastq.gz -r $OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.R2.fastq.gz $OUTPUT_DIR/$BASENAME/novo_assembly/iva
cp $OUTPUT_DIR/$BASENAME/novo_assembly/iva/contigs.fasta $OUTPUT_DIR/$BASENAME/novo_assembly/contigs.fasta

#### Filter contigs that map against the reference
blastn -query $OUTPUT_DIR/$BASENAME/novo_assembly/contigs.fasta -subject ${r} -outfmt 6 > $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.blast6
awk '{print $1}' $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.blast6  | sort | uniq > $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.list
if [[ $(awk 'NR==1{print $9}' $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.blast6) -lt  $(awk 'NR==1{print $10}' $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.blast6) ]]
then
	seqtk subseq $OUTPUT_DIR/$BASENAME/novo_assembly/contigs.fasta $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.list > $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta
	echo "+"
else
	seqtk subseq $OUTPUT_DIR/$BASENAME/novo_assembly/contigs.fasta $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.list | seqtk seq -r -a -  > $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta
	echo "-"
fi

### Change name of fasta header
mv $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta  $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs_multiple.fasta
reformat.sh in=$OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs_multiple.fasta out=$OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta reads=1

#Circularization
circules.py -k 120 -f $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta
circules.py -f $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta -l $(python /home/luna.kuleuven.be/u0117540/BK/bin/circules.py -k 120 -f $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs.fasta | tail -2 - | head -n 1 - | sed -E "s/.*\-\ (.*)\ \(.*/\1/g" -) -p $OUTPUT_DIR/$BASENAME/novo_assembly/circules
mv $OUTPUT_DIR/$BASENAME/novo_assembly/circules.circular* $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs_circular.fasta
sed -i "s/^>.*/>$BASENAME denovo-assembly/gi" $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs_circular.fasta

### Translocate beguinning of the sequence to end, as Dunlop reference coordiantes
cat $OUTPUT_DIR/$BASENAME/novo_assembly/polyoma_contigs_circular.fasta | \
	sed 's/TTTTGCAAAA/@@TTTTGCAAAA##/' | sed 's/^\(.*@@\)\(TTTTGCAAAA\)\(##.*\)$/\2\3\1/' | \
	sed 's/@@//' | sed 's/##//' | sed 's/TTTTGCAAA/@@@TTTTGCAAA/2' | sed 's/\(.*\)\(@@@TTTTGCAAAA.*\)/\1/' | sed 's/@@//' > $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta

### Polishing
bwa index $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta
bwa mem -t $THREADS  $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R1.fastq.gz \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.R2.fastq.gz > $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.sam

samtools view -@ $THREADS -b $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.sam | samtools sort -l 9 -@ $THREADS -T $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
	-o $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.sort.bam
samtools index $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.sort.bam

pilon-1.23.jar --genome $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta --bam $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.sort.bam \
	--outdir $OUTPUT_DIR/$BASENAME/novo_assembly/pilon --output $BASENAME.denovo.genome.coord.corrected --changes --fix all --threads $THREADS --verbose

seqtk seq -A $OUTPUT_DIR/$BASENAME/novo_assembly/pilon/$BASENAME.denovo.genome.coord.corrected.fasta | sed 's/_pilon//' >  $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.pilon.fasta

### Move to results
cp $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.pilon.fasta $OUTPUT_DIR/$BASENAME/results
cp $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.pilon.fasta $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta

### Calculate genome length
NOVO_LENGTH=$(awk -v FS="" '!/^>/ {print NF}'  $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.pilon.fasta)

## Reference-based assembly
echo "Reference-based assembly"

### Map reads against the reference
bwa index ${r}
bwa mem -t $THREADS  ~/BK/ref/polyoma_dunlop_reference.fasta \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R1.fastq.gz \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R2.fastq.gz > $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sam

### SAM to BAM to sort BAM
samtools view -@ $THREADS -b $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sam | \
	samtools sort -l 9 -@ $THREADS \
	-T ${r} \
	-o $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.bam

### Index reference
echo "Index reference by lofreq"
lofreq faidx ${r}

### Correct alignments in BAM
lofreq  viterbi \
	-f ${r} \
	$OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.bam |\
	samtools sort - > $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.bam

### Add SNV (including indel) qualities to BAM and re-sort BAM
lofreq alnqual -b -r \
	$OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.bam \
	${r} | samtools sort - > $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.alnqual.bam

### Index corrected and quality-annotated BAM
lofreq  index $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.alnqual.bam

### Run low-frquency variant calling in parallel
lofreq call-parallel --pp-threads $THREADS \
	-f ${r} \
	-o $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.vcf \
	--call-indels -s -a 0.01 --use-orphan $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.alnqual.bam
bcftools mpileup --max-depth 10000 -Ou -f ${r} $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.sort.vit.alnqual.bam | bcftools call -mv -Ov -o $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.vcf

### Compress the called variant file
bgzip -f $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.vcf

### Index the compress variant file
tabix -f $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.vcf.gz

### Compute the consensus by most frequent allel per position (sites are independent)
bcftools consensus -f ${r} \
	$OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.vsDunlop.vcf.gz -o $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.Dunlop.fasta
seqtk seq -A $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.Dunlop.fasta  > $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta

### Calculate genome length
REF_LENGTH=$(awk -v FS="" '!/^>/ {print NF}' $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta)

### Change the header to the appropiate for the sample
sed -i "s/>.*/>$BASENAME ref-assembly/" $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta

### Copy reference-based assembly consensus sequence at results folder
cp $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta $OUTPUT_DIR/$BASENAME/results

# ANNOTATION
## Transfer annotation from BK Dunlop reference to BK de novo assembly
(pushd $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/ && <PATH-TO-RATT>/start.ratt.sh <PATH-TO-EMBL-REF-FOLDER> ../$BASENAME.denovo.genome.coord.pilon.fasta $BASENAME Species.Repetitive; popd) ##

### Fix EMBL
sed "1s/.*/ID   GBK_ACC; 1; circular; DNA; HTG; VRL; $NOVO_LENGTH BP./" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.embl > $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.mod.embl
sed -i -e "s/\.5.../@@@/" -e "s/@@@/\.$NOVO_LENGTH/" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.mod.embl

### Convert EMBL to Genbank
embl2gbk.pl $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.mod.embl $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gbk

### Extract CDS from the newly annotated BK genome
genbank2fasta.py -i $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.mod.embl -m embl -s nt -f CDS -q note -o $BASENAME.gene.nt.fasta

### Extrat proteins from the newly annotated BK genome
genbank2fasta.py -i $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.$BASENAME.final.mod.embl -m embl -s aa -f CDS -q note -o $BASENAME.gene.aa.fasta

### Modify header names
sed -i "s/ /_/g" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.nt.fasta
sed -i "s/>/>${BASENAME}_/g" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.nt.fasta
sed -i "s/ /_/g" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.aa.fasta
sed -i "s/>/>${BASENAME}_/g" $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.aa.fasta

## Transfer annotation from BK Dunlop reference to BK de novo assembly
(pushd $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/ && /home/luna.kuleuven.be/u0117540/BK/bin/ratt-code/start.ratt.sh /home/luna.kuleuven.be/u0117540/BK/ref/embl/ ../$BASENAME.ref.genome.vsDunlop.fasta $BASENAME Strain; popd) ##

### Fix EMBL
sed -i "1s/.*/ID   GBK_ACC; 1; circular; DNA; HTG; VRL; $NOVO_LENGTH BP./" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.$BASENAME.final.embl

### Convert EMBL to Genbank
sed -i -e "s/\.5.../@@@/" -e "s/@@@/\.$REF_LENGTH/" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.$BASENAME.final.embl
embl2gbk.pl $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.$BASENAME.final.embl $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gbk

### Extract CDS from the newly annotated BK genome
genbank2fasta.py -i $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.$BASENAME.final.embl -m embl -s nt -f CDS -q note -o $BASENAME.gene.nt.fasta

### Extrat proteins from the newly annotated BK genome
genbank2fasta.py -i $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.$BASENAME.final.embl -m embl -s aa -f CDS -q note -o $BASENAME.gene.aa.fasta

### Modify header names
sed -i "s/ /_/g" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.nt.fasta
sed -i "s/>/>${BASENAME}_/g" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.nt.fasta
sed -i "s/ /_/g" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.aa.fasta
sed -i "s/>/>${BASENAME}_/g" $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.aa.fasta

### Copy annotation files and gene fastas to results
cp $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gbk $OUTPUT_DIR/$BASENAME/results/$BASENAME.denovo.genome.coord.gbk
cp $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gbk $OUTPUT_DIR/$BASENAME/results/$BASENAME.ref.genome.vsDunlop.gbk
cp $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.nt.fasta $OUTPUT_DIR/$BASENAME/results/$BASENAME.denovo.genome.coord.gene.nt.fasta
cp $OUTPUT_DIR/$BASENAME/novo_assembly/annotation/$BASENAME.gene.aa.fasta $OUTPUT_DIR/$BASENAME/results/$BASENAME.denovo.genome.coord.gene.aa.fasta
cp $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.nt.fasta $OUTPUT_DIR/$BASENAME/results/$BASENAME.ref.genome.vsDunlop.gene.nt.fasta
cp $OUTPUT_DIR/$BASENAME/ref_assembly/annotation/$BASENAME.gene.aa.fasta $OUTPUT_DIR/$BASENAME/results/$BASENAME.ref.genome.vsDunlop.gene.aa.fasta

# VARIANT CALLING
## VC for de novo assembled BK genome
### Index the de novo BK genome
bwa index $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta

### Map high-quality reads to the indexed de novo BK genome
bwa mem -t $THREADS $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R1.fastq.gz \
	$OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R2.fastq.gz > $OUTPUT_DIR/$BASENAME/novo_assembly/map/$BASENAME.vsdenovo.sam

### SAM to BAM to sort BAM
samtools view -@ $THREADS -b $OUTPUT_DIR/$BASENAME/novo_assembly/map/$BASENAME.vsdenovo.sam | \
       samtools sort -l 9 -@ $THREADS \
       -T $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
       -o $OUTPUT_DIR/$BASENAME/novo_assembly/map/$BASENAME.vsdenovo.sort.bam

### Index reference
lofreq faidx $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta

### Correct alignments in BAM
lofreq  viterbi \
       -f $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
       $OUTPUT_DIR/$BASENAME/novo_assembly/map/$BASENAME.vsdenovo.sort.bam | \
       samtools sort - > $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.sort.vit.bam

### Add SNV (including indel) qualities to BAM and re-sort BAM
lofreq  alnqual -b -r \
       $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.sort.vit.bam \
       $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta | samtools sort - > $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.sort.vit.alnqual.bam

### Index corrected and quality-annotated BAM
lofreq  index $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.sort.vit.alnqual.bam

### Run low-frquency variant calling in parallel
lofreq  call-parallel --pp-threads $THREADS \
       -f $OUTPUT_DIR/$BASENAME/novo_assembly/$BASENAME.denovo.genome.coord.fasta \
       -o $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.vcf \
       --call-indels -b dynamic -s -a 0.001 --use-orphan $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.sort.vit.alnqual.bam

### Copy VC file to results folder
cp $OUTPUT_DIR/$BASENAME/novo_assembly/lofreq/$BASENAME.vsdenovo.vcf $OUTPUT_DIR/$BASENAME/results/

## VC for reference assembled BK genome
### Map high-quality reads to the indexed reference assembled BK genome
bwa index $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta
bwa mem -t $THREADS $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta \
        $OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R1.fastq.gz \
        $OUTPUT_DIR/$BASENAME/trim/$BASENAME.trim.min500.target1000.R2.fastq.gz > $OUTPUT_DIR/$BASENAME/ref_assembly/map/$BASENAME.vsref.sam

### SAM to BAM to sort BAM
samtools view -@ $THREADS -b $OUTPUT_DIR/$BASENAME/ref_assembly/map/$BASENAME.vsref.sam | \
       samtools sort -l 9 -@ $THREADS \
       -T $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta \
       -o $OUTPUT_DIR/$BASENAME/ref_assembly/map/$BASENAME.vsref.sort.bam

### Index reference
lofreq faidx $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta

### Correct alignments in BAM
lofreq  viterbi \
       -f $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta \
       $OUTPUT_DIR/$BASENAME/ref_assembly/map/$BASENAME.vsref.sort.bam | \
       samtools sort - > $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.sort.vit.bam

### Add SNV (including indel) qualities to BAM and re-sort BAM
lofreq  alnqual -b -r \
	$OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.sort.vit.bam \
	$OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta| samtools sort - > $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.sort.vit.alnqual.bam

### Index corrected and quality-annotated BAM
lofreq  index $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.sort.vit.alnqual.bam

### Run low-frquency variant calling in parallel
lofreq call-parallel --pp-threads $THREADS \
       -f $OUTPUT_DIR/$BASENAME/ref_assembly/$BASENAME.ref.genome.vsDunlop.fasta \
       -o $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.vcf \
       --call-indels -b dynamic -s -a 0.001 --use-orphan $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.sort.vit.alnqual.bam

### Copy VC file to results folder
cp $OUTPUT_DIR/$BASENAME/ref_assembly/lofreq/$BASENAME.vsref.vcf $OUTPUT_DIR/$BASENAME/results/
