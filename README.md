\# LTR\_FINDER\_parallel: LTR Discovery, Density, and LAI (Wheat Genomes)



This guide shows how to install and run \*\*LTR\_FINDER\_parallel\*\* to discover LTR retrotransposons, compute \*\*LTR density\*\* in 1 Mb bins, and estimate \*\*LAI-like summaries\*\* from GFF3 outputs on HPC systems (e.g., Atlas/CCAST).



---



\## 1) Install \& Set PATH



```bash

\# Clone

git clone https://github.com/oushujun/LTR\_FINDER\_parallel.git

cd LTR\_FINDER\_parallel/



\# (Optional) locate the bundled ltr\_finder binary

find . -name "ltr\_finder" -type f



\# Test the binary directly

./bin/LTR\_FINDER.x86\_64-1.0.7/ltr\_finder -h



\# Add to PATH for this session

export PATH=$PATH:/directory/where/this/saved/LTR\_FINDER\_parallel/bin/LTR\_FINDER.x86\_64-1.0.7

ltr\_finder -h



\# Make PATH permanent

echo 'export PATH=$PATH:/directory/where/this/saved/LTR\_FINDER\_parallel/bin/LTR\_FINDER.x86\_64-1.0.7' >> ~/.bashrc

source ~/.bashrc

```



\# Run LTR\_FINDER\_parallel (interactive)

```bash

ml perl/5.38.0



cd /directory/where/this/saved/LTR\_FINDER\_parallel/



\# Example: Glenn

perl LTR\_FINDER\_parallel \\

&nbsp; -seq /directory/where/this/saved/wheat.fasta \\

&nbsp; -size 10000000 \\

&nbsp; -threads 48 \\

&nbsp; -harvest\_out \\

&nbsp; -time 300

```

\# SLURM Batch Jobs (Atlas)

```bash

\#!/bin/bash

\#SBATCH -N 1

\#SBATCH --cpus-per-task=10

\#SBATCH -p atlas

\#SBATCH --mem=300GB

\#SBATCH -J ltr\_fparallel\_wheat

\#SBATCH -A genolabswheatphg



ml perl/5.38.0

export PATH=$PATH:/directory/where/this/saved/LTR\_FINDER\_parallel/bin/LTR\_FINDER.x86\_64-1.0.7



cd /directory/where/this/saved/LTR\_FINDER\_parallel/



perl LTR\_FINDER\_parallel \\

&nbsp; -seq /directory/where/this/saved/wheat.fasta \\

&nbsp; -size 10000000 -threads 10 -harvest\_out -time 300

```



\# Compute LTR Density (1 Mb windows)

```bash

\# 1 Mb window density (bp per bin) for long\_terminal\_repeat features

awk -v window=1000000 '

$3 == "long\_terminal\_repeat" {

&nbsp;   chr = $1

&nbsp;   start = $4

&nbsp;   end = $5

&nbsp;   bin = int(start / window) \* window

&nbsp;   density\[chr, bin] += (end - start + 1)

}

END {

&nbsp;   for (key in density) {

&nbsp;       split(key, keys, SUBSEP)

&nbsp;       chr = keys\[1]

&nbsp;       bin\_start = keys\[2]

&nbsp;       bin\_end = bin\_start + window

&nbsp;       print chr, bin\_start, bin\_end, density\[key]

&nbsp;   }

}' wheat\_ltr\_retrotransposon.gff3 > wheat\_ltr\_density\_circos.txt



\# Rename chr\* to ta\* (e.g., chr1A -> ta1A)

sed -E 's/chr(\[1-7]\[ABD])/ta\\1/' wheat\_ltr\_density\_circos.txt > wheat\_ltr\_density\_circos\_renamed.txt



\# Get min/max density

awk 'BEGIN {max=0; min=1e18} {if ($4>max) max=$4; if ($4<min) min=$4} END {print "Max Density:", max; print "Min Density:", min}' wheat\_ltr\_density\_circos\_renamed.txt

```



\#  LTR-Length Summaries (LAI-like ratio)

```bash

cd /directory/where/this/saved/LTR\_FINDER\_parallel



\# Extract only LTR\_retrotransposon entries

grep -w "LTR\_retrotransposon" wheat.fasta.finder.combine.gff3 > ltr\_lai\_wheat.gff3



\# Total LTR length (sum of feature lengths)

awk '{sum += ($5 - $4)} END {print sum}' ltr\_lai\_sumai3.gff3



\# Per-chromosome LTR length

awk '{len\[$1] += ($5 - $4)} END {for (c in len) print c, len\[c]}' ltr\_lai\_sumai3.gff3 | sort -k1,1V

```



Maintainer:



Ruby Mijan

