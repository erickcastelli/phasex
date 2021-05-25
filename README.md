# phasex
Written by Erick C. Castelli, erick.castelli@unesp.br

Current version: 0.8.1

Phasex is piace of software written in C++ to automate and compare multiple haplotyping runs using Shapeit4 and/or Beagle 4.

It has been used mainly for haplotyping of HLA and KIR alleles in different studies, but it can be used for haplotyping of other genes. It is suitable for datasets of thosands of samples but limited number os variants (e.g., 5000 samples, 3000 variants). Phasex will automate parallel runs and compare the results, fixing the haplotypes with concordance rates over a threshold to subsequent runs.

## How to compile and install phasex:
First, you need Boost to compile phasex. On Ubuntu, we recommend using sudo-apt-get install libboost-all-dev. On Macos, we recommend use homebrew: brew install boost

Second, we use cmake to generate the make file. Thus, be sure that you have cmake installed.


Download this git repository, and follow these instrocutions:
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
Assuming that you have moved the binary to /usr/local/bin or onother directory in the PATH, follow these instructions:
> phasex

You will see the main functions
                                                                     
     phase-ps: to phase variants considering Phase Sets (PS)                      
     recreate: to recreate the final VCF in case you have edited the final results                                          
     hp-ps: to recode GATK ReadBackedPhasing format to PS format
     
> phasex hp-ps

You will see the options for the hp-ps function. This method converts GATK ReadBackedPhasing phased VCF to the PS format to be compatible with phasex. If you use WhatsHap, this step is not necessary if you used "--tag=PS" when running WhatsHap.

To use this function, you must type:
> phasex hp-ps vcf=THE_VCF_FILE output=THE_NEW_VCF_FILE

Do not use spaces before and after "=". The --quiet mode forces the program to not output any comment.




