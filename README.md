# 🎮 Atari

**A**lterna**T**e **A**llele **R**ead v**I**sualizer

![Version](https://img.shields.io/badge/version-0.1.1-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

**Atari** is a command-line tool to count reference and alternate alleles from BAM files across genomic regions and visualize allele depths.  
It uses `pysam` for BAM parsing, `matplotlib` for plotting, and supports flexible input with BED files or manual region specification.

---

## ✨ Features

- Count reference and alternate reads per position
- Aggregate counts into non-overlapping windows
- Apply quality filters (mapping, base quality, depth)
- Output clean TSV tables
- Generate high-quality plots
- Annotate visualizations with genes from GTF files
- Multi-sample comparisons
- Detailed per-base alternate allele breakdown (A, C, G, T, DEL, INS, N)

---

## 📦 Installation

```bash
git clone https://github.com/yourusername/atari.git
cd atari
pip install -r requirements.txt
```

Requirements:

- Python ≥ 3.7
- pysam
- matplotlib
- numpy
- pandas
- tqdm
- rich

---

## 🚀 Quick Start

> Ensure BAM and FASTA files are indexed (.bai, .fai).

### 1. Basic Usage

```bash
python atari.py \
  --bam sample1.bam sample2.bam \
  --reference genome.fa \
  --bed regions.bed \
  --window 100 \
  --output results.tsv
```
### 2. Full Pipeline with Visualization and Gene Annotation

```bash
python atari.py \
  --bam sample1.bam sample2.bam \
  --reference genome.fa \
  --bed regions.bed \
  --window 100 \
  --output results.tsv \
  --plot \
  --plot-output coverage_plot.png \
  --gtf annotations.gtf \
  --verbose
```
### 3. Analyze Specific Chromosomal Region

```bash
python atari.py \
  --bam sample.bam \
  --reference genome.fa \
  --chromosome chr20 \
  --start 1000000 \
  --end 1100000 \
  --window 100 \
  --output chr20_results.tsv \
  --plot \
  --plot-output chr20_coverage.png
```

### 4. Detailed Base-Level Counts
```bash
python atari.py \
  --bam sample.bam \
  --reference genome.fa \
  --bed regions.bed \
  --window 100 \
  --output detailed_results.tsv \
  --detailed \
  --verbose
```

| Argument | Description | Required |
| --- | --- | --- |
| --bam | BAM file(s) (indexed) | ✅ |
| --reference | Reference FASTA (indexed) | ✅ |
| --bed | BED file with regions | ➡️ (or --chromosome) |
| --chromosome | Chromosome name (e.g., chr1) | ➡️ (or --bed) |
| --start, --end | Region coordinates (if using --chromosome) | ✅ |
| --window | Window size in bp for summarization | ✅ |
| --output | Output file name (TSV) | ✅ |
| --plot | Generate plot | Optional |
| --plot-output | Plot output file | Optional |
| --gtf | Gene annotations (GTF format) | Optional |
| --verbose | Verbose mode | Optional |
| ...and many more! (see --help) |  |  |


### 🐛 Troubleshooting

1. Ensure your BAM and FASTA files are indexed (.bai and .fai present).
2. Use `--verbose` to print detailed processing steps.
3. For debugging errors, set ATARI_DEBUG=1 environment variable.

### 👨‍💻 Contributing
PRs welcome! Feel free to open issues for bugs, feature requests, or improvements.

### 📫 Contact
For questions, contact samarth8392@gmail.com.