### novoBreak v1.1.3rc

Author: Zechen Chong

Email: chongzechen@gmail.com or kchen3@mdanderson.org

Draft date: Aug. 15, 2018

## Description

novoBreak is a tool used in cancer genomic studies to discover SV (both somatic and germline)
breakpoints. It can report accurate breakpoints of Deletions (DEL), Duplications (DUP), Inversions (INV) and
Translocations (TRA) (you should consider some of them are mobile elements insertions or templated
insertions). For novel insertions, we may only report the breakpoints but not the inserted
sequence. Please forget about novel insertions at the moment. We will work on that later. novoBreak
was designed for Illumina paired-end data.

Please cite novoBreak:
> Chong, Zechen, et al. "novoBreak: local assembly for breakpoint detection in cancer genomes." Nature Methods (2016) [PMID: 27892959](https://www.nature.com/nmeth/journal/v14/n1/full/nmeth.4084.html).	


## System requirements and dependency

novoBreak runs on a x86_64 Linux system with a ~40GB physical memory. It depends on SSAKE
for local assembly, bwa-mem for contig mapping and samtools and picard(SamToFastq) to extract
reads. I have already put them in the release. If you cannot run the dependencies, please download
the lastest versions of these tools.

## Installation

```
git clone https://github.com/czc/nb_distribution.git
```
Then, please also add this directory to your PATH:
```
export PATH=$PWD/nb_distribution/:$PATH
```
Type
```
novoBreak
```
You should get the following output:

```
novoBreak - a tool for discovering somatic sv breakpoints
Auther: Zechen Chong <zchong@mdanderson.org>
Version: 1.1 (r20151007)
Usage:
  novoBreak -i <tumorbam> -c <normalbam> -r <reference> -o <output.kmer> [options]
Options:
  -h             This help
  -i <string>    Tumor bam file
  -c <string>    Normal bam file
  -r <string>    Reference file in fasta format
  -k <int>       Kmer size, <=31 [31]
  -o <string>    Output kmer
  -g <int>       Output germline events [0]
  -m <int>       Minimum kmer count regarded as novo kmers [3].
```

Otherwise, please install novoBreak first
```
git clone https://github.com/czc/novobreak_src.git
cd novobreak_src && make
```
**Replace novoBreak in nb_distribution with the compiled novoBreak***


## Usage

Preprocessing: For virus integration analysis, please add the virus genome(s) to reference and realign the reads to the 
              new reference. The integration should be a TRA event.
Input files: a reference sequence file and two bam files (tumor and normal) of mapping results of
			 the paired-end reads
Output file: a filtered high confident VCF file (novoBreak.pass.flt.vcf). I left the intermedium
             files. You can delete them or inspect them as you need.

run novoBreak:
*The shell script for easy (hopefully) run of novoBreak is in the release directory. You can tune the
parameters as you wish.*
```
bash <A_PATH>/novoBreak/run_novoBreak.sh <novoBreak_exe_dir> <ref> <tumor_bam> <normal_bam> <n_CPUs:INT> [outputdir:-PWD]
```

## About the default filter

To increase sensitivity, novoBreak tries to infer as many SVs as possible based on the local assembly results. But many of the inferred
SVs may be false positives due to misassembly or lack of enough evidence. So we provided a default filter to get a relatively stringent
filtered callset based on real data experience. We empirically defined the minimum SV size as 100 bp and no upper limit. Users can change
the filter and cutoffs based on the utility and the knowledge as needed. An empirical filter can be made based on the column 6 of novoBreak's
output. A higher value of column 6 indicates a more reliable event. But sometimes it can be real event when breakpoints fall in the repetitive
regions but with a low value.

## NEWS

* 20150814: Updated cluster module. Improved speed a lot. 
* 20150924: Introduced multiple-core for the shell script. Further improved speed. 
* 20151009: novoBreak can directly read from bam files. No raw fastq files required. Save lots of space. Should have done earlier!
* 20151015: Added an output redirection option
* 20160126: Removed Picard's SamToFastq and changed to samtools 1.3. Added breakpoint consensus information.
* 20160412: Reduced memory consumption. No short reads realignments required. Fixed a few bugs.
* 20160915: Added header to the final output and provided an alternative filter "filter_sv2.pl".
* 20160815: Fixed header by adding the '##END' line and fixed 'CT' values in VCF output. Thanks Xiaotong Yao!
