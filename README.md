# LTR_FINDER_parallel: LTR Discovery, Density, and LAI (Wheat Genomes)

This guide shows how to install and run **LTR_FINDER_parallel** to discover LTR retrotransposons, compute **LTR density** in 1 Mb bins, and estimate **LAI-like summaries** from GFF3 outputs on HPC systems (e.g., Atlas/CCAST).

---

## 1) Install & Set PATH
```bash
# Clone
git clone https://github.com/oushujun/LTR_FINDER_parallel.git
cd LTR_FINDER_parallel/

# (Optional) locate the bundled ltr_finder binary
find . -name "ltr_finder" -type f

# Test the binary directly
./bin/LTR_FINDER.x86_64-1.0.7/ltr_finder -h

# Add to PATH for this session
export PATH=$PATH:/directory/where/this/saved/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7
ltr_finder -h

# Make PATH permanent
echo 'export PATH=$PATH:/directory/where/this/saved/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7' >> ~/.bashrc
source ~/.bashrc
```

## 2) Run LTR_FINDER_parallel (interactive)
```bash
ml perl/5.38.0
cd /directory/where/this/saved/LTR_FINDER_parallel/

# Example: Glenn
perl LTR_FINDER_parallel \
  -seq /directory/where/this/saved/wheat.fasta \
  -size 10000000 \
  -threads 48 \
  -harvest_out \
  -time 300
```

## 3) SLURM Batch Jobs (Atlas)
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH -p atlas
#SBATCH --mem=300GB
#SBATCH -J ltr_fparallel_wheat
#SBATCH -A genolabswheatphg

ml perl/5.38.0
export PATH=$PATH:/directory/where/this/saved/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7
cd /directory/where/this/saved/LTR_FINDER_parallel/

perl LTR_FINDER_parallel \
  -seq /directory/where/this/saved/wheat.fasta \
  -size 10000000 -threads 10 -harvest_out -time 300
```

## 4) Compute LTR Density (1 Mb windows)
```bash
# 1 Mb window density (bp per bin) for long_terminal_repeat features
awk -v window=1000000 '
$3 == "long_terminal_repeat" {
    chr = $1
    start = $4
    end = $5
    bin = int(start / window) * window
    density[chr, bin] += (end - start + 1)
}
END {
    for (key in density) {
        split(key, keys, SUBSEP)
        chr = keys[1]
        bin_start = keys[2]
        bin_end = bin_start + window
        print chr, bin_start, bin_end, density[key]
    }
}' wheat_ltr_retrotransposon.gff3 > wheat_ltr_density_circos.txt

# Rename chr* to ta* (e.g., chr1A -> ta1A)
sed -E 's/chr([1-7][ABD])/ta\1/' wheat_ltr_density_circos.txt > wheat_ltr_density_circos_renamed.txt

# Get min/max density
awk 'BEGIN {max=0; min=1e18} {if ($4>max) max=$4; if ($4<min) min=$4} END {print "Max Density:", max; print "Min Density:", min}' wheat_ltr_density_circos_renamed.txt
```

## 5) LTR-Length Summaries (LAI-like ratio)
```bash
cd /directory/where/this/saved/LTR_FINDER_parallel

# Extract only LTR_retrotransposon entries
grep -w "LTR_retrotransposon" wheat.fasta.finder.combine.gff3 > ltr_lai_wheat.gff3

# Total LTR length (sum of feature lengths)
awk '{sum += ($5 - $4)} END {print sum}' ltr_lai_wheat.gff3

# Per-chromosome LTR length
awk '{len[$1] += ($5 - $4)} END {for (c in len) print c, len[c]}' ltr_lai_wheat.gff3 | sort -k1,1V
```

Maintainer:

Ruby Mijan
