
Overview      {#mainpage}
===========

_align2bed_ is a light-weight tool that extracts single nucleotide polymorphisms (SNPs) from FASTA alignments and saves them to a [plink](http://zzz.bwh.harvard.edu/plink/index.shtml) binary-format BED file. It is fast and can process whole-genome data on a laptop in minutes. The FASTA alignments should be in the format provided by the _Drosophila_ Genome Nexus (<http://www.johnpool.net/genomes.html>), i.e. the data for each line (or individual) are in a separate file with nucleotides in one line and without the customary FASTA header. Each chromosome arm has its own set of files.

The sofware uses control files that list paths to each FASTA file to be processed. An example data set is included to illustrate the necessary features of the data. One of the individuals must be marked as the reference, and it is also assumed that this reference is the outgroup. The reference genotypes are not included in the output. BED format files require the SNPs to be biallelic, so SNPs that do not meet this criterion are not included as output. However, if a SNP is biallelic within the population (non-reference) sample, but both alleles are different from the outgroup, it is included. Names of such SNPs are marked with "d" at the end (in the accompanying .bim file) to enable downstream filtering. In addition, SNPs with outgroup missing are included but their names marked with "m" at the end. SNP names are "s" plus position, then "m" or "d" if applicable, then underscore ("_"), then chromosome arm name. Note that the ancestral nucleotide (if available) is listed last in the .bim file.

While _align2bed_ is tailored for the _Drosophila_ Genome Nexus data, there are three ways it can be extended to similar data sets from other species. Data can be arranged to mimic the _Drosophila_ set by treating chromosomes in groups of five. Alignment length can vary indefinitely. Furthermore, I wrote the program using a class that has wider applicability. Taking the _align2bed_ source code as an exmaple, and reading the provided interface documentation, someone with even very limited experience in C++ can write software that applies to different data sets and hardware configurations. Finally, anyone who would like to extend functionality even further is welcome to modify the class implementation to suit their needs.

Neither _align2bed_ itself nor the class used to implement it has any dependencies outside of the C++ STL. Only a compiler capable of recognizing the C++11 standard is required (I successfully compiled with LLVM and GCC). The implementation is multithreaded, with chromosome arms processed in parallel. Each thread allocates a 2Gb buffer to read the FASTA files. The _align2bed_ source can be easily modified to change the threading and memory allocation parameters (see included class documentation for details).

To compile, make sure you are in the directory with the source code files and run

	g++ align2bed.cpp sequence.cpp -o align2bed -O3 -march=native -std=c++11

then copy the binary where you need it. Run by typing `./align2bed` in the directory with the binary and the data, or move into an appropriate /bin folder for global access.

The example data set includes control files for each autosome and 20 kb of alignments extracted from 284 _Drosophila_ lines (283 _D. melanogaster_ and a _D. simulans_ outgroup). Each chromosome needs a separate control file, which simply lists the paths to FASTA files (which can include directories), one file per line. The file containing the outgroup sequence should be marked with "r:".
