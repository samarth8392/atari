"""ATARI - AlternaTe Allele Read vIsualizer."""

from .core import count_alleles, process_regions, process_windows
from .visualize import plot_allele_counts, plot_allele_counts_with_annotations, plot_gene_annotations
from .utils import get_total_regions_size, parse_bed_file, parse_gtf_file

__version__ = '0.1.1'