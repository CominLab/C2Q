# C2Q

 Alignment-free genome comparison based on Sequencing Reads and Quality Values

Sequencing technologies are generating enormous amounts of read data, however the assembly of genomes and metagenomes, from fragments, is still one of the most challenging task. If a reference genome is not available the basic task of reads mapping cannot be carried out. For metagenomic communities the problem is even more challenging because we can not know in advance which genomes are present in the metagenomic sample. In this paper we study the comparison of genomes and metagenomes based on read data, without assembly, using alignment-free sequence statistics.
Quality scores produced by sequencing platforms are fundamental for various analysis, like reads mapping and error detection. Moreover future-generation sequencing platforms, e.g. PacBio and MinION, will produce long reads but with an error rate around 15\%. In this context it will be fundamental to exploit quality value information within the framework of alignment-free statistics.
In this paper we present a family of alignment-free measures, called $d^q$-type, that are based on $k$-mer counts and quality values. These statistics can be used to compare genomes and metagenomes based on their read sets, without the need to assembly. Results on simulated genomes show that the evolutionary relationship of genomes can be reconstructed based on the direct comparison of sets of reads. The use of quality values on average improves the classification accuracy, and its contribution increases when the reads are more noisy. Preliminary experiments on microbial communities show that similar metagenomes can be detected just by processing their read data.


Software

The file contains source code (file c2q.cpp) and example reads files (reads subdirectory).


To compile 'make' can be used or the line (the former includes -O3 optimization flag)
c++ c2q.cpp -o c2q

The program c2q requires 4 mandatory arguments

c2q [k] [list] [type] [M]

k k-mer length

list a file containing a list of input sequences file in the form [fast file] [seq name]

type either Q for quality value based output or A for non quality output

M length of the reads (needed to compute averages)

To run the provided example with k=5 use the command

c2q 5 reads.txt Q 300

The software outputs (along with several support files) the 4 matrices d2,d2*,d2q and d2q* in 4 files with extension .out


Licence

The software is freely available for academic use.
For questions about the tool, please contact Matteo Comin.

Reference

Please cite the following paper:

M.Comin, M.Schimd
"Fast Comparison of Genomic and Meta-Genomic Reads with Alignment-Free Measures based on Quality Values",
BMC Medical Genomics 2016, 9(Suppl 1):36 Open Access
