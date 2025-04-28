#!/usr/bin/env python3
import argparse
import sys
import os
from rich.console import Console

# Internal imports from the package
from .core import process_regions, process_windows
from .utils import parse_bed_file, parse_gtf_file
from .visualize import plot_allele_counts, plot_allele_counts_with_annotations
from . import __version__  # Import version from __init__.py

# Set up environment variables
os.environ['XDG_RUNTIME_DIR'] = '/tmp/runtime-dir'
os.makedirs('/tmp/runtime-dir', exist_ok=True)

console = Console()
ascii_banner = f"""
[cyan]
░█████╗░████████╗░█████╗░██████╗░██╗
██╔══██╗╚══██╔══╝██╔══██╗██╔══██╗██║
███████║░░░██║░░░███████║██████╔╝██║
██╔══██║░░░██║░░░██╔══██║██╔══██╗██║
██║░░██║░░░██║░░░██║░░██║██║░░██║██║
╚═╝░░╚═╝░░░╚═╝░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝
[/cyan]

[bold]A[/bold]lterna[bold]T[/bold]e [bold]A[/bold]llele [bold]R[/bold]ead v[bold]I[/bold]sualizer 

[magenta] [bold] Version:  {__version__} [/bold] [/magenta] 
"""

description=f"""
A tool to count reference and alternate alleles at genomic positions using samtools via pysam.
This script analyzes BAM files and reports the number of reads supporting reference and alternate
alleles at each position within specified regions, then visualizes the results.
"""

usage_epilog=f"""
Examples:
  # Basic usage with BED file
  atari \\
    --bam sample1.bam sample2.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output results.tsv

  # Full pipeline with visualization and gene annotations
  atari \\
    --bam sample1.bam sample2.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output results.tsv \\
    --plot \\
    --plot-output coverage_plot.png \\
    --gtf annotations.gtf \\
    --verbose
    
  # Generate detailed base-level counts
  atari \\
    --bam sample1.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output detailed_results.tsv \\
    --detailed \\
    --verbose

  # Analyze a specific region with advanced visualization options
  atari \\
    --bam sample.bam \\
    --reference genome.fa \\
    --chromosome chr20 \\
    --start 1000000 \\
    --end 1100000 \\
    --window 100 \\
    --output chr20_results.tsv \\
    --plot \\
    --plot-output chr20_coverage.png \\
    --gtf annotations.gtf \\
    --figsize 14,10 \\
    --dpi 600 \\
    --max-genes 15

  # Multi-sample comparison with quality filters
  atari \\
    --bam tumor.bam normal.bam \\
    --reference genome.fa \\
    --bed hotspot_regions.bed \\
    --window 50 \\
    --min-mapping-quality 30 \\
    --min-base-quality 30 \\
    --min-depth 20 \\
    --output comparison_results.tsv \\
    --plot \\
    --plot-output sample_comparison.png
    """

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Count reference and alternate alleles from BAM files and visualize results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=usage_epilog)
    
    # Required arguments
    parser.add_argument('-b', '--bam', required=True, nargs='+',
                        help='Input BAM file(s). Must be indexed.')
    parser.add_argument('-r', '--reference', required=True,
                        help='Reference genome in FASTA format. Must be indexed.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file name (TSV format).')
    parser.add_argument('-w', '--window', required=True, type=int, 
                    help='Window size (bp) for calculating mean ref/alt depths in non-overlapping windows.')
    
    # Region specification (mutually exclusive)
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument('-B', '--bed',
                            help='BED file with regions to analyze.')
    region_group.add_argument('-c', '--chromosome',
                            help='Chromosome to analyze (e.g., chr1).')
    
    # Optional parameters for manual region specification
    parser.add_argument('-s', '--start', type=int,
                        help='Start position (1-based, inclusive). Required if --chromosome is used.')
    parser.add_argument('-e', '--end', type=int,
                        help='End position (1-based, inclusive). Required if --chromosome is used.')
    
    # Optional filters
    parser.add_argument('-q', '--min-mapping-quality', type=int, default=20,
                        help='Minimum mapping quality score (default: 20).')
    parser.add_argument('-Q', '--min-base-quality', type=int, default=20,
                        help='Minimum base quality score (default: 20).')
    parser.add_argument('-d', '--min-depth', type=int, default=0,
                        help='Minimum read depth to report a position (default: 0).')
    
    # Output options
    parser.add_argument('-D', '--detailed', action='store_true',
                        help='Include detailed counts for each alternate base (A,C,G,T,DEL,INS).')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress information.')
    parser.add_argument('--no-progress-bar', action='store_true',
                        help='Disable progress bar.')
    
    # Visualization options
    visualization_group = parser.add_argument_group('Visualization Options')
    visualization_group.add_argument('--plot', action='store_true',
                        help='Generate a visualization plot of the results.')
    visualization_group.add_argument('--plot-output',
                        help='Output file for visualization (e.g., plot.png). Required if --plot is specified.')
    visualization_group.add_argument('--dpi', type=int, default=300,
                        help='DPI (resolution) for output image (default: 300).')
    visualization_group.add_argument('--figsize', type=str, default="12,8",
                        help='Figure size in inches, format "width,height" (default: "12,8").')
    visualization_group.add_argument('--no-show-plot', action='store_true',
                        help='Do not display plot (only save if --plot-output is specified).')
    visualization_group.add_argument('--gtf',
                        help='GTF file with gene annotations for visualization.')
    visualization_group.add_argument('--max-genes', type=int, default=10,
                        help='Maximum number of genes to display in each region (default: 10).')
    
    args = parser.parse_args()

    # Validate arguments
    if args.chromosome and (args.start is None or args.end is None):
        parser.error("--start and --end are required when --chromosome is specified.")
    
    if args.start is not None and args.end is not None and args.start > args.end:
        parser.error("Start position must be less than or equal to end position.")
    
    for bam_file in args.bam:
        if not os.path.exists(bam_file):
            parser.error(f"BAM file not found: {bam_file}")
        if not os.path.exists(f"{bam_file}.bai") and not os.path.exists(f"{os.path.splitext(bam_file)[0]}.bai"):
            parser.error(f"Index not found for BAM file: {bam_file}")
    
    if not os.path.exists(args.reference):
        parser.error(f"Reference file not found: {args.reference}")
    if not os.path.exists(f"{args.reference}.fai"):
        parser.error(f"Index not found for reference file: {args.reference}")
    
    if args.bed and not os.path.exists(args.bed):
        parser.error(f"BED file not found: {args.bed}")
    
    if args.plot and not args.plot_output and args.no_show_plot:
        parser.error("--plot-output is required when --plot and --no-show-plot are both specified.")
    
    if args.gtf and not os.path.exists(args.gtf):
        parser.error(f"GTF file not found: {args.gtf}")
    
    return args

def run_main():
    """Run the main processing functions with the parsed arguments."""
    args = parse_arguments()
    
    # Parse regions
    if args.bed:
        regions = parse_bed_file(args.bed)
    else:
        regions = [{
            'chrom': args.chromosome,
            'start': args.start,
            'end': args.end
        }]
    
    if args.verbose:
        console.print(f"Reference genome: [bold]{args.reference}[/bold]", style="green")
        console.print(f"Number of BAM files: [bold]{len(args.bam)}[/bold]", style="green")
        console.print(f"Number of regions: [bold]{len(regions)}[/bold]", style="green")
        console.print(f"Total positions to analyze: [bold]{get_total_regions_size(regions) * len(args.bam):,}[/bold]", style="green")
    
    # Process regions
    results_df = process_regions(
        args.bam,
        args.reference,
        regions,
        args.min_mapping_quality,
        args.min_base_quality,
        args.min_depth,
        args.verbose,
        not args.no_progress_bar
    )
    
    if results_df.empty:
        console.print("No positions found with the specified criteria.", style="bold red")
        sys.exit(1)

    # Process windows if requested
    if args.window:
        if args.verbose:
            console.print(f"Grouping results into windows of size [bold]{args.window}[/bold] bp", style="blue")
        
        # Create windowed summary data (this replaces the position-level data)
        window_df = process_windows(results_df, args.window)
        
        if window_df.empty:
            console.print("No windows found with the specified criteria.", style="bold red")
            sys.exit(1)
        
        # Define window-specific columns
        if args.detailed:
            columns = [
                'Bam_file', 'CHROM', 'Window_start', 'Window_end', 'Positions_count',
                'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads',
                'Mean_A_count', 'Mean_C_count', 'Mean_G_count', 'Mean_T_count', 
                'Mean_DEL_count', 'Mean_INS_count', 'Mean_N_count'
            ]
        else:
            columns = [
                'Bam_file', 'CHROM', 'Window_start', 'Window_end', 'Positions_count',
                'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads'
            ]
            
        results_df = window_df
    else:
        # Position-level output (original behavior)
        if args.detailed:
            columns = [
                'Bam_file', 'CHROM', 'POS', 'REF', 'DP', 'Ref_reads', 'Alt_reads',
                'A_count', 'C_count', 'G_count', 'T_count', 'DEL_count', 'INS_count', 'N_count'
            ]
        else:
            columns = ['Bam_file', 'CHROM', 'POS', 'REF', 'DP', 'Ref_reads', 'Alt_reads']

    # Write output TSV
    results_df[columns].to_csv(args.output, sep='\t', index=False)

    if args.verbose:
        console.print(f"[bold green]Results written to {args.output}[/bold green]")
        if args.window:
            console.print(f"Total windows processed: [bold]{len(results_df)}[/bold]", style="blue")
        else:
            console.print(f"Total positions processed: [bold]{len(results_df)}[/bold]", style="blue")
    
    # Generate visualization if requested
    if args.plot:
        if args.verbose:
            console.print("Generating visualization...", style="magenta")
        
        # Parse GTF file if provided
        gtf_features = None
        if args.gtf:
            if args.verbose:
                console.print(f"Parsing gene annotations from [bold]{args.gtf}[/bold]", style="magenta")
            gtf_features = parse_gtf_file(args.gtf)
        
        # Parse figure size
        try:
            figsize = tuple(float(x) for x in args.figsize.split(','))
            if len(figsize) != 2:
                figsize = (12, 8)
        except:
            figsize = (12, 8)
            
        # Plot with window data and gene annotations
        if gtf_features:
            plot_allele_counts_with_annotations(
                results_df,
                gtf_features=gtf_features,
                output_file=args.plot_output,
                dpi=args.dpi,
                figsize=figsize,
                show_plot=not args.no_show_plot,
                max_genes=args.max_genes
            )
        else:
            # Fall back to original plotting function if no GTF
            plot_allele_counts(
                results_df,
                output_file=args.plot_output,
                dpi=args.dpi,
                figsize=figsize,
                show_plot=not args.no_show_plot
            )
        
        if args.verbose and args.plot_output:
            console.print(f"[bold green]Visualization saved to {args.plot_output}[/bold green]")
    
    return 0

def main():
    """Entry point for the command line interface."""
    # Print the banner
    console.print(ascii_banner)
    console.print(description)
    
    # Parse arguments
    args = parse_arguments()
    
    # Run the main function
    try:
        return run_main(args)
    except KeyboardInterrupt:
        console.print("\n[bold red]Process interrupted by user[/bold red]")
        return 1
    except Exception as e:
        console.print(f"[bold red]Error: {str(e)}[/bold red]")
        if os.environ.get('ATARI_DEBUG'):
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())