# phasex
Written by Erick C. Castelli, erick.castelli@unesp.br

Current version: 0.8.1

Phasex is a software written in C++ to automate and compare multiple haplotyping runs using Shapeit4 and/or Beagle 4.

It has been used mainly for haplotyping of HLA and KIR alleles in different studies. Still, it can be used for haplotyping of other genes. It is suitable for datasets of thousands of samples but a limited number of variants (e.g., 5000 samples, 3000 variants). Phasex will automate parallel runs and compare the results, fixing the haplotypes with concordance rates over a threshold to subsequent runs.

## How to compile and install phasex:
First, you need Boost to compile phasex. On Ubuntu, we recommend using sudo-apt-get install libboost-all-dev. On Macos, we recommend use homebrew: brew install boost

Second, we use cmake to generate the make file. Thus, be sure that you have cmake installed.


Download this git repository, and follow these instructions:
- enter folder /build
- type "cmake .."
- if everything worked, type "make"
- the binary will be placed in the folder /build

Finally, you need a working copy of Shapeit4 (https://odelaneau.github.io/shapeit4/) and Beagle 4.1 (https://faculty.washington.edu/browning/beagle/b4_1.html)

## Manuscripts using phasex:
- Castelli et al. Immunogenetics of resistance to SARS-CoV-2 infection in discordant couples. MedRxiv (doi 10.1101/2021.04.21.21255872)
- Naslavsky et al. Whole-genome sequencing of 1,171 elderly admixed individuals from the largest Latin
American metropolis (São Paulo, Brazil). BioRxiv (doi 10.1101/2020.09.15.298026)
- Sonon et al. Peripheral spectrum neurological disorder after arbovirus infection is associated with HLA-F variants among Northeastern Brazilians. Infect Genet Evol. 2021 Apr 8;92:104855. doi: 10.1016/j.meegid.2021.104855
- Sonon et al. Human leukocyte antigen (HLA)-F and -G gene polymorphisms and haplotypes are associated with malaria susceptibility in the Beninese Toffin children. Infect Genet Evol. 2021 Mar 27;92:104828. doi: 10.1016/j.meegid.2021.104828
- Weiss et al. KIR2DL4 genetic diversity in a Brazilian population sample: implications for transcription regulation and protein diversity in samples with different ancestry backgrounds. Immunogenetics. 2021 Jun;73(3):227-241. doi: 10.1007/s00251-021-01206-9
- Souza et al. Hla-C genetic diversity and evolutionary insights in two samples from Brazil and Benin. HLA. 2020 Oct;96(4):468-486. doi: 10.1111/tan.13996
- Ramos et al. A large familial cluster and sporadic cases of frontal fibrosing alopecia in Brazil reinforce known human leucocyte antigen (HLA) associations and indicate new HLA susceptibility haplotypes. J Eur Acad Dermatol Venereol. 2020 Oct;34(10):2409-2413. doi: 10.1111/jdv.16629
- and many others...

## How to use phasex:
PHASEX used shapeit4 to phase bi-allelic variants (considering the PS field) and BEAGLE 4.1 to phase multiallelic variants considering the scaffold inferred by Shapeit4.

Assuming that you have moved the binary to /usr/local/bin or another directory in the PATH, follow these instructions:
> phasex

You will see the main functions
                                                                     
     phase-ps: to phase variants considering Phase Sets (PS)                      
     recreate: to recreate the final VCF in case you have edited the final results                                          
     hp-ps: to recode GATK ReadBackedPhasing format to PS format
     
> phasex hp-ps

You will see the options for the hp-ps function. This method converts GATK ReadBackedPhasing phased VCF to the PS format to be compatible with phasex. If you use WhatsHap, this step is not necessary if you used "--tag=PS" when running WhatsHap.

To use this function, you must type:
> phasex hp-ps vcf=THE_VCF_FILE output=THE_NEW_VCF_FILE

Do not use spaces before and after "=". The --quiet mode forces the program not to output any comment.

> phasex phase-ps

You will see all the options for this function. This is a typical phasex run:

> phasex phase-ps vcf=VCF_FILE_IN_PS_FORMAT iterations=10 replicates=20 shapeit=SHAPEIT4_BINARY beagle=BEAGLE4_JAR

This configuration informs PHASEX to run 20 parallel haplotyping runs (replicates), fixing concordant haplotypes using the threshold value (95% of the runs), and performing these steps 10 times (iteractions). The output folder will be placed next to the input VCF unless modified with "output=". 

The threshold for fixing a haplotype is 95% (the default), i.e., a haplotype is fixed as true if 19 runs (20 replicates * 0.95) indicate the same hapotype for a sample. You can modify this using "threshold=".

Phasex uses half of the number of cores of the system unless modified by "threads=".

Phasex will perform 10 iterations, i.e., 10 steps of 20 parallel runs and haplotype comparison. You can modify this using "iterations=" and "replicates=".

After the final iteration (in this case, the 10th interation), phasex will output the final haplotypes, considering only samples in which the same haplotype was inferred in at least 70% of the replicates in the final run. You can modify this using "select=".

Option scheme is used by Shapeit4. By default, this scheme is 15b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,15m. you can modify this using scheme="10b,1p,1b,1p,1b,1p,1b,1p,10m".

Option shapeit_others is used to indicate other shapeit4 parameters.

Option map is used to indicate a genetic map for Shapeit4 (not mandatory). Please download these maps at the Shapeit4 website.

Flag --quiet forces PHASEX not to output any progress or comment.

Flag --biallelic forces PHASEX to deal only with biallelic variants, using only Shapeit4.


## The phasex outputs:

The output structure is as follows:

phasex.log: Record all the parameters and some quality-control information

results.vcf: This is the final PHASED VCF file. Only the samples passing the select threshold are included in this file (by default: 70% of the runs presenting the same haplotype in the final run).

results.freq: The haplotypes, their global count, and frequency

sample_list.txt: The list of samples that passed the SELECT threshold.

/shapeit : the shapeit results for each iteration, and the final results in "results.txt"

/shapeit/results.txt: the final results when using shapeit4. This file presents the following format:

Sample	h1	h2	Freq(1)	Info(1)	Freq(2)	Info(2)	Freq(n)	Info(n)	Status
 - Sample: the sample id
 - h1: first haplotype
 - h2: second haplotype
 - Freq(n): proportion of parallel runs indicating this pair of haplotypes in iteration N
 - info(n): "-" if under the threshold, "def" if fixed for the next iteration
 - Status: "-" if not this haplotype pair is under the SELECT threshold, "pass" if it is above the SELECT threhold. Only the samples with "pass" are included in the final VCF.
																					
/beagle : the beagle results for each iteration, and the final results in "results.txt"

/beagle/results.txt: the final results when using Beagle. Same format as for Shapeit.

/source : the files used for the haplotyping procedure.





