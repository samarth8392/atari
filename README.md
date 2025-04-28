# ATARI - AlternaTe Allele Read vIsualizer

ATARI is a tool for analyzing genomic data from BAM files. It counts reference and alternate alleles at genomic positions and visualizes the results.

## Installation

### From PyPI

```bash
# From PyPI
pip install atari

# From Conda
conda install -c bioconda atari

# Clone the repository
git clone https://github.com/yourusername/atari.git
cd atari

# Option 1: Install with pip
pip install .

# Option 2: Create conda environment from environment.yml
conda env create -f environment.yml
conda activate atari

# usage
# Basic usage with BED file
atari \
  --bam sample1.bam sample2.bam \
  --reference genome.fa \
  --bed regions.bed \
  --window 100 \
  --output results.tsv

# Full pipeline with visualization and gene annotations
atari \
  --bam sample1.bam sample2.bam \
  --reference genome.fa \
  --bed regions.bed \
  --window 100 \
  --output results.tsv \
  --plot \
  --plot-output coverage_plot.png \
  --gtf annotations.gtf \
  --verbose