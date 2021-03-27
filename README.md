# BKAnaLite
Genome assembly and Variant Calling for Polyoma BK Virus (PyBKV).

## Installation
To install the pipeline just:

```
git clone https://github.com/joanmarticarreras/BKAnaLite
cd BKAnaLite
./BKAnaLite.sh --help
```

BKAnaLite depends on a series of dependencies, which executables should be in `$PATH`:

```
BBMap
IVA
BLAST
SeqTK
Python2.7
Python3
BWA
SAMTools
LoFreq
bgzip
tabix
bcftools
RATT
```
Additionally, `BKAnaLite` depends on the script [circules.py](https://github.com/chrishah/MITObim/blob/master/misc_scripts/circules.py) from the MITObim package. Accessory scripts from RATT should be in `$PATH` too (NOTE: Depend on `Python2.7`).

Absolute paths to `RATT` and the annotated reference genome in EMBL format should be added in L209. Path to a file containing the sequencing adapters and PCR primers (if used) should be added in L80.

## Use

```
Usage: BKAnaLite.sh [OPTIONS]
-r     Path to reference genome in FASTA format                                                                                                                                                        
-1     Path to R1 Illunmina PE reads                                                                                                                                                                   
-2     Path to R2 Illumina PE reads                                                                                                                                                                    
-s     Path to Sanger reads in FASTQ format                                                                                                                                                            
-t     Number of threads to use when multithreading is possible                                                                                                                                        
-x     Basename for the output files                                                                                                                                                                   
-o     Path to output directory                                                                                                                                                                        
-h     help (this output)      
```
