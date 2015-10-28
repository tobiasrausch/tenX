SV Re-genotyping Methods for 10X Genomics
=========================================

Installing
----------

`git clone --recursive https://github.com/tobiasrausch/tenX.git`

`cd tenX/`

`make all`

Example
-------

Download the NA12878 10X Genomics phased bam file from http://10xgenomics.com/ and the 1000 Genomes structural variant (SV) release from the EBI ftp site:

`wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz'`

Run the haplotype-aware SV re-genotyping of the SV site list:

`./src/genoDEL -v ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz NA12878.10X.bam > del.tsv`

`python compareGeno.py -t del.tsv > deletion.summary`
