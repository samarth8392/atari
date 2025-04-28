import os
import sys
import pandas as pd
import pysam
import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from tqdm import tqdm
from rich.console import Console

from .utils import get_total_regions_size, parse_bed_file, parse_gtf_file

console = Console()

def count_alleles(bam_file: str, reference_fasta: str, chromosome: str, start_pos: int, end_pos: int, 
                  min_mapping_quality: int = 20, min_base_quality: int = 20, 
                  progress_callback=None) -> List[Dict[str, Any]]:
    """
    Count reference and alternate alleles at each position in a region.
    
    Args:
        bam_file: Path to the BAM file
        reference_fasta: Path to the reference genome
        chromosome: Chromosome name
        start_pos: Start position (1-based)
        end_pos: End position (1-based)
        min_mapping_quality: Minimum mapping quality
        min_base_quality: Minimum base quality
        progress_callback: Function to call for progress updates with position as argument
    """
    # Open files
    samfile = pysam.AlignmentFile(bam_file, "rb")
    fastafile = pysam.FastaFile(reference_fasta)
    results = []
    
    # Check if chromosome exists in both BAM and reference
    try:
        test_fetch = samfile.fetch(chromosome, start_pos-1, start_pos)
        next(test_fetch, None)
    except (ValueError, StopIteration):
        return results
    
    try:
        fastafile.fetch(chromosome, start_pos-1, start_pos)
    except (ValueError, KeyError):
        return results
    
    # Process each position
    for pos in range(start_pos, end_pos + 1):
        # Call progress callback for each position
        if progress_callback:
            progress_callback()
            
        # Get reference base
        try:
            ref_base = fastafile.fetch(chromosome, pos-1, pos).upper()
        except (ValueError, KeyError):
            continue
        
        # Initialize counters
        base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'DEL': 0, 'INS': 0}
        coverage = 0
        
        # Pileup at the current position
        for pileupcolumn in samfile.pileup(chromosome, pos-1, pos, truncate=True, min_base_quality=min_base_quality):
            if pileupcolumn.pos != pos-1:
                continue
                
            coverage = pileupcolumn.n
            
            # Process each read
            for pileupread in pileupcolumn.pileups:
                # Skip low quality reads
                if pileupread.alignment.mapping_quality < min_mapping_quality:
                    continue
                
                # Count deletions
                if pileupread.is_del:
                    base_counts['DEL'] += 1
                    continue
                
                # Count insertions
                if pileupread.indel > 0:
                    base_counts['INS'] += 1
                
                # Get base
                if pileupread.query_position is None:
                    continue
                    
                base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                
                # Count base
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1
        
        # Calculate reference and alternate allele counts
        ref_count = base_counts.get(ref_base, 0)
        total_alt_count = sum(count for base, count in base_counts.items() 
                             if base != ref_base and base != 'N' and count > 0)
        
        # Add result
        result = {
            'Bam_file': os.path.basename(bam_file),
            'CHROM': chromosome,
            'POS': pos,
            'REF': ref_base,
            'DP': coverage,
            'Ref_reads': ref_count,
            'Alt_reads': total_alt_count,
            'A_count': base_counts['A'],
            'C_count': base_counts['C'],
            'G_count': base_counts['G'],
            'T_count': base_counts['T'],
            'DEL_count': base_counts['DEL'],
            'INS_count': base_counts['INS'],
            'N_count': base_counts['N'],
        }
        
        results.append(result)
    
    # Close files
    samfile.close()
    fastafile.close()
    
    return results


def process_regions(bam_files: List[str], reference_fasta: str, regions: List[Dict[str, Any]], 
                    min_mapping_quality: int, min_base_quality: int, min_depth: int, 
                    verbose: bool, show_progress: bool) -> pd.DataFrame:
    """Process all specified regions across all BAM files."""
    all_results = []
    
    # Calculate the total work - number of positions to process
    total_positions = get_total_regions_size(regions) * len(bam_files)
    positions_processed = 0
    
    # Create progress bar based on total positions
    progress_bar = None
    if show_progress:
        progress_bar = tqdm(
            total=total_positions,
            desc="Processing genomic positions",
            unit="pos",
            bar_format="{desc}: {percentage:3.1f}% |{bar}| {n_fmt}/{total_fmt} positions"
        )
    
    # Process each BAM file and region
    for i, bam_file in enumerate(bam_files, 1):
        if verbose:
            console.print(f"Processing BAM file {i}/{len(bam_files)}: {bam_file}", style="blue")
        
        for j, region in enumerate(regions, 1):
            if verbose and not show_progress:
                console.print(f"  Region {j}/{len(regions)}: {region['chrom']}:{region['start']}-{region['end']}", style="blue")
            
            # Update progress function
            def update_progress():
                nonlocal positions_processed
                positions_processed += 1
                if progress_bar:
                    progress_bar.update(1)
            
            # Process the region
            results = count_alleles(
                bam_file, 
                reference_fasta, 
                region['chrom'], 
                region['start'], 
                region['end'], 
                min_mapping_quality, 
                min_base_quality,
                update_progress
            )
            
            all_results.extend(results)
    
    # Close progress bar
    if progress_bar:
        progress_bar.close()
    
    # Convert results to DataFrame
    if not all_results:
        return pd.DataFrame()
        
    df = pd.DataFrame(all_results)
    
    # Filter by minimum depth
    if min_depth > 0:
        df = df[df['DP'] >= min_depth]
    
    return df


def process_windows(df: pd.DataFrame, window_size: int) -> pd.DataFrame:
    """
    Group results into non-overlapping windows and calculate mean values.
    
    Args:
        df: DataFrame with position-level results
        window_size: Size of the windows in base pairs
    
    Returns:
        DataFrame with window-level summaries
    """
    if df.empty:
        return df
    
    # Add a window ID column based on position
    df['Window'] = (df['POS'] - 1) // window_size
    
    # Group by BAM file, chromosome, and window
    window_groups = df.groupby(['Bam_file', 'CHROM', 'Window'])
    
    # Calculate window summaries
    window_results = []
    
    for (bam, chrom, window), group in window_groups:
        # Calculate window bounds
        start_pos = window * window_size + 1
        end_pos = (window + 1) * window_size
        
        # Calculate means
        mean_ref_reads = group['Ref_reads'].mean()
        mean_alt_reads = group['Alt_reads'].mean()
        mean_depth = group['DP'].mean()
        
        # Add detailed base counts if they exist
        result = {
            'Bam_file': bam,
            'CHROM': chrom,
            'Window_start': start_pos,
            'Window_end': end_pos,
            'Mean_depth': mean_depth,
            'Mean_ref_reads': mean_ref_reads,
            'Mean_alt_reads': mean_alt_reads,
            'Positions_count': len(group)
        }
        
        # Include mean base counts if detailed output is requested
        if 'A_count' in group.columns:
            result.update({
                'Mean_A_count': group['A_count'].mean(),
                'Mean_C_count': group['C_count'].mean(),
                'Mean_G_count': group['G_count'].mean(),
                'Mean_T_count': group['T_count'].mean(),
                'Mean_DEL_count': group['DEL_count'].mean(),
                'Mean_INS_count': group['INS_count'].mean(),
                'Mean_N_count': group['N_count'].mean(),
            })
        
        window_results.append(result)
    
    return pd.DataFrame(window_results)

