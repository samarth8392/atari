import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
from typing import Dict, List, Optional, Tuple, Any
from matplotlib.ticker import FuncFormatter

def plot_gene_annotations(ax, features, start_pos, end_pos, max_genes=10):
    """Plot gene annotations on the given axes.
    
    Args:
        ax: Matplotlib axes to plot on
        features: List of gene features for the current chromosome
        start_pos: Start position for the current view
        end_pos: End position for the current view
        max_genes: Maximum number of genes to display
        
    Returns:
        Boolean indicating if a legend was added
    """
    # Filter features for the visible region
    visible_genes = {}
    visible_transcripts = {}
    visible_exons = {}
    
    for feature in features:
        # Skip features outside the visible region
        if feature['end'] < start_pos or feature['start'] > end_pos:
            continue
            
        if feature['type'] == 'gene':
            gene_id = feature['gene_id']
            visible_genes[gene_id] = feature
            
        elif feature['type'] == 'transcript':
            transcript_id = feature['transcript_id']
            gene_id = feature['gene_id']
            
            if gene_id not in visible_transcripts:
                visible_transcripts[gene_id] = []
                
            visible_transcripts[gene_id].append(feature)
            
        elif feature['type'] in ['exon', 'CDS']:
            transcript_id = feature['transcript_id']
            gene_id = feature['gene_id']
            
            if gene_id not in visible_exons:
                visible_exons[gene_id] = []
                
            visible_exons[gene_id].append(feature)
    
    # Limit the number of genes if there are too many
    gene_ids = list(visible_genes.keys())
    if len(gene_ids) > max_genes:
        # Sort genes by length (show longest genes)
        gene_ids.sort(key=lambda g: visible_genes[g]['end'] - visible_genes[g]['start'], reverse=True)
        gene_ids = gene_ids[:max_genes]
        
    # No genes to display
    if not gene_ids:
        ax.text(0.5, 0.5, 'No gene annotations in view', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_yticks([])
        return False
        
    # Plot genes
    y_pos = 0
    y_positions = {}
    colors = plt.cm.tab10.colors
    
    # For legend - collect items to add to the legend
    legend_elements = []
    
    for i, gene_id in enumerate(gene_ids):
        gene = visible_genes[gene_id]
        color = colors[i % len(colors)]
        y_pos = i
        y_positions[gene_id] = y_pos
        
        # Plot gene as a thick line
        gene_start = max(gene['start'], start_pos)
        gene_end = min(gene['end'], end_pos)
        
        # Gene body
        gene_line = ax.plot([gene_start, gene_end], [y_pos, y_pos], 
                color=color, linewidth=2, solid_capstyle='butt')[0]
        
        # Direction arrow
        if gene['strand'] == '+':
            ax.arrow(gene_end - (gene_end - gene_start) * 0.05, y_pos, 
                    (gene_end - gene_start) * 0.03, 0, 
                    head_width=0.2, head_length=(gene_end - gene_start) * 0.02, 
                    fc=color, ec=color)
        else:
            ax.arrow(gene_start + (gene_end - gene_start) * 0.05, y_pos, 
                    -(gene_end - gene_start) * 0.03, 0, 
                    head_width=0.2, head_length=(gene_end - gene_start) * 0.02, 
                    fc=color, ec=color)
        
        # Add to legend items for this gene
        if i == 0:  # Only add one gene example to the legend
            legend_elements.append(plt.Line2D([0], [0], color=color, lw=2, 
                                            label='Gene body'))
            # Add arrow style to legend 
            legend_elements.append(plt.Line2D([0], [0], marker='>',
                                         color='w', markerfacecolor=color, 
                                         markersize=8, label='Direction'))
        
        exon_legend_added = False
        cds_legend_added = False
        
        # Plot exons if available
        if gene_id in visible_exons:
            exons = visible_exons[gene_id]
            for exon in exons:
                exon_start = max(exon['start'], start_pos)
                exon_end = min(exon['end'], end_pos)
                
                if exon['type'] == 'CDS':
                    # CDS as taller box
                    cds_patch = plt.Rectangle(
                        (exon_start, y_pos - 0.3), 
                        exon_end - exon_start, 
                        0.6, 
                        facecolor=color, 
                        edgecolor='none', 
                        alpha=0.8)
                    ax.add_patch(cds_patch)
                    
                    # Add to legend once
                    if i == 0 and not cds_legend_added:
                        legend_elements.append(plt.Rectangle((0, 0), 1, 1, 
                                              facecolor=color, alpha=0.8, 
                                              label='CDS (coding)'))
                        cds_legend_added = True
                else:
                    # Regular exon as box
                    exon_patch = plt.Rectangle(
                        (exon_start, y_pos - 0.2), 
                        exon_end - exon_start, 
                        0.4, 
                        facecolor=color, 
                        edgecolor='none', 
                        alpha=0.6)
                    ax.add_patch(exon_patch)
                    
                    # Add to legend once
                    if i == 0 and not exon_legend_added:
                        legend_elements.append(plt.Rectangle((0, 0), 1, 1, 
                                              facecolor=color, alpha=0.6, 
                                              label='Exon'))
                        exon_legend_added = True
        
        # Add gene name
        ax.text(gene_start, y_pos + 0.4, gene['gene_name'], 
               fontsize=8, ha='left', va='bottom', 
               bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0))
    
    # Set y-axis limits and labels
    ax.set_ylim(-1, len(gene_ids))
    ax.set_yticks([])
    
    # Format x-axis in Mbp
    def format_kbp(x, pos):
        return f'{x/1_000:.2f}'
    
    ax.xaxis.set_major_formatter(FuncFormatter(format_kbp))
    ax.set_xlabel('Position (Kbp)')
    
    # Add legend for gene features
    if legend_elements:
        ax.legend(handles=legend_elements, 
                 loc='upper right', 
                 fontsize='small',
                 framealpha=0.7,
                 title="Gene Features")
        return True
    
    return False

def plot_allele_counts_with_annotations(df, gtf_features=None, output_file=None, dpi=300, 
                                      figsize=(12, 8), show_plot=True, max_genes=10,
                                      min_coverage_threshold=None, highlight_regions=None):
    """
    Generate an improved visualization of allele counts with gene annotations.
    
    Args:
        df: DataFrame with allele count data
        gtf_features: Dictionary of gene features by chromosome from parse_gtf_file
        output_file: Path to save the plot (if None, will display plot)
        dpi: Resolution for the output image
        figsize: Figure size as a tuple (width, height) in inches
        show_plot: Whether to display the plot (in addition to saving it)
        max_genes: Maximum number of genes to display per region
        min_coverage_threshold: Optional minimum coverage threshold to display as horizontal line
        highlight_regions: Optional list of dict with 'start', 'end', and 'label' keys for regions to highlight
    """
    # Check if the required columns exist
    required_cols = ['CHROM', 'Window_start', 'Window_end', 'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads']
    if not all(col in df.columns for col in required_cols):
        print(f"DataFrame must contain columns: {', '.join(required_cols)}", file=sys.stderr)
        return
    
    # Get unique BAM files and chromosomes
    bam_files = df['Bam_file'].unique()
    chromosomes = df['CHROM'].unique()
    
    # Set up the plot
    num_bams = len(bam_files)
    num_chroms = len(chromosomes)
    
    # Create a custom color palette with more distinct colors
    palette = {
        'depth': '#1f77b4',  # Blue
        'ref': '#2ca02c',    # Green
        'alt': '#d62728',    # Red
        'ratio': '#9467bd'   # Purple
    }
    
    # Adjust layout to include ONE gene annotation panel per chromosome
    # Modified to have a single gene annotation panel per chromosome
    if num_bams == 1 and num_chroms == 1:
        fig = plt.figure(figsize=figsize)
        gs = plt.GridSpec(2 + num_bams, 1, height_ratios=[3] * num_bams + [1, 1], figure=fig)
        
        # Coverage plot
        ax_cov = fig.add_subplot(gs[0])  
        # Ratio plot
        ax_ratio = fig.add_subplot(gs[1], sharex=ax_cov)
        # Single gene annotation plot
        ax_genes = fig.add_subplot(gs[2], sharex=ax_cov)
        
        axes_cov = [[ax_cov]]
        axes_ratio = [[ax_ratio]]
        axes_genes = [[ax_genes]]  # One genes panel
    elif num_bams >= 1 and num_chroms == 1:
        # One chromosome, multiple BAMs
        fig = plt.figure(figsize=(figsize[0], figsize[1] * (num_bams + 0.5)))
        # Each BAM gets a coverage and ratio plot, plus one shared gene plot at the bottom
        gs = plt.GridSpec(2 * num_bams + 1, 1, 
                         height_ratios=[3, 1] * num_bams + [1], 
                         figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            # Coverage plot for this BAM
            ax_cov = fig.add_subplot(gs[2*i])
            # Ratio plot for this BAM
            ax_ratio = fig.add_subplot(gs[2*i+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
        
        # Single gene annotation at the bottom for this chromosome
        ax_genes = fig.add_subplot(gs[2*num_bams], sharex=axes_cov[0][0])
        # Store in a 2D array for consistency
        axes_genes = [[ax_genes]]
    elif num_bams == 1 and num_chroms > 1:
        # Multiple chromosomes, one BAM
        fig = plt.figure(figsize=(figsize[0], figsize[1] * num_chroms))
        gs = plt.GridSpec(3 * num_chroms, 1, height_ratios=[3, 1, 1] * num_chroms, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        axes_genes = []
        
        for j in range(num_chroms):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[3*j])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[3*j+1], sharex=ax_cov)
            # Gene annotation for this chromosome
            ax_genes = fig.add_subplot(gs[3*j+2], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
            axes_genes.append([ax_genes])
    else:
        # Multiple BAMs and chromosomes
        fig = plt.figure(figsize=(figsize[0] * num_chroms, figsize[1] * (num_bams + 0.5)))
        
        # Calculate heights: each BAM gets a coverage and ratio plot, 
        # plus one gene plot per chromosome
        gs = plt.GridSpec(2 * num_bams + 1, num_chroms, 
                         height_ratios=[3, 1] * num_bams + [1], 
                         figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            axes_cov_row = []
            axes_ratio_row = []
            
            for j in range(num_chroms):
                # Coverage plot
                ax_cov = fig.add_subplot(gs[2*i, j])
                # Ratio plot
                ax_ratio = fig.add_subplot(gs[2*i+1, j], sharex=ax_cov)
                
                axes_cov_row.append(ax_cov)
                axes_ratio_row.append(ax_ratio)
            
            axes_cov.append(axes_cov_row)
            axes_ratio.append(axes_ratio_row)
        
        # Create one gene annotation row per chromosome at the bottom
        axes_genes = [[fig.add_subplot(gs[2*num_bams, j], sharex=axes_cov[0][j])] 
                      for j in range(num_chroms)]
    
    # Format function for x-axis (millions of base pairs)
    def format_kbp(x, pos):
        return f'{x/1_000:.2f}'
    
    # Process each chromosome to determine plot ranges
    chrom_data = {}
    for j, chrom in enumerate(chromosomes):
        # Get all data for this chromosome across all BAMs
        chrom_subset = df[df['CHROM'] == chrom]
        if chrom_subset.empty:
            continue
            
        x_min = chrom_subset['Window_start'].min()
        x_max = chrom_subset['Window_end'].max()
        chrom_data[chrom] = {
            'x_min': x_min,
            'x_max': x_max
        }
    
    # Create each subplot
    for i, bam in enumerate(bam_files):
        for j, chrom in enumerate(chromosomes):
            subset = df[(df['Bam_file'] == bam) & (df['CHROM'] == chrom)].copy()
            
            if subset.empty:
                continue
                
            # Calculate window midpoints for x-axis
            subset['window_mid'] = (subset['Window_start'] + subset['Window_end']) / 2
            
            # Calculate alt/total ratio
            subset['alt_ratio'] = (subset['Mean_alt_reads'] / (subset['Mean_depth'].replace(0, np.nan))) * 100  # Convert to percentage
            subset['alt_ratio'] = subset['alt_ratio'].fillna(0)  # Replace NaNs with 0
            
            # Get current axes
            ax_cov = axes_cov[i][j] if num_bams > 1 or num_chroms > 1 else axes_cov[0][0]
            ax_ratio = axes_ratio[i][j] if num_bams > 1 or num_chroms > 1 else axes_ratio[0][0]
            # For genes, use the chromosome index only since we share gene plots across BAMs
            ax_genes = axes_genes[j][0] if num_chroms > 1 else axes_genes[0][0]
            
            # Get axis boundaries for highlighting
            x_min = subset['Window_start'].min()
            x_max = subset['Window_end'].max()
            
            # Plot coverage data with enhanced styling
            # Add shaded background for alt reads area
            ax_cov.set_facecolor('#f9f9f9')  # Light background
            ax_cov.axvspan(x_min, x_max, color='#fff3f3', alpha=0.5)  # Light red shading for alt reads area
            
            # Add grid for better readability
            ax_cov.grid(True, linestyle='--', alpha=0.3)
            
            # Draw coverage threshold if specified
            if min_coverage_threshold is not None:
                ax_cov.axhline(y=min_coverage_threshold, color='#aaaaaa', linestyle='--', alpha=0.7, 
                              label=f'Min. Coverage ({min_coverage_threshold}x)')
            
            # Plot main coverage data with both lines and markers
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], '-', 
                       color=palette['depth'], label='Mean Depth', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], 'o', 
                       color=palette['depth'], markersize=4, alpha=0.5)
            
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], '-', 
                       color=palette['ref'], label='Ref Reads', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], 's', 
                       color=palette['ref'], markersize=4, alpha=0.5)
            
            # Calculate reasonable y-max (max depth + 10%)
            y_max = max(subset['Mean_depth'].max(), subset['Mean_ref_reads'].max()) * 1.1
            
            # Add alternate allele subplot with different scale and clearer styling
            alt_color = palette['alt']
            ax2 = ax_cov.twinx()
            # Add light red background to make it clear this is a different axis
            ax2.patch.set_alpha(0.0)  # Make background transparent
            
            # Add lines and markers for alt reads
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '-', 
                    color=alt_color, linewidth=1.5, label='Alt Reads')
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '^', 
                    color=alt_color, markersize=4, alpha=0.5)
            
            # Set alt reads scale
            alt_max = max(subset['Mean_alt_reads'].max() * 1.5, 5)  # Ensure non-zero max
            ax2.set_ylim(0, alt_max)
            
            # Styling for alt reads axis - Fix for Issue 2: No box around the label
            alt_label = 'Alternate Reads (Right Axis)'
            ax2.set_ylabel(alt_label, color=alt_color, fontweight='bold')
            ax2.tick_params(axis='y', labelcolor=alt_color)
            
            # Add annotation to clarify alt reads is on right axis
            ax_cov.annotate('Alt Reads →', xy=(0.98, 0.5), xycoords='axes fraction',
                          ha='right', va='center', color=alt_color, fontsize=9)
            
            # Plot ratio data
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], '-', 
                        color=palette['ratio'], label='Alt/Total Ratio', linewidth=1.5)
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], 'o', 
                        color=palette['ratio'], markersize=4, alpha=0.5)
            ax_ratio.set_ylim(0, min(100.0, subset['alt_ratio'].max() * 1.5 or 10.0))
            ax_ratio.set_ylabel('Alt/Total Ratio (%)')
            ax_ratio.grid(True, linestyle='--', alpha=0.3)
            
            # Highlight regions of interest if specified
            if highlight_regions:
                for region in highlight_regions:
                    if region['start'] <= x_max and region['end'] >= x_min:  # Only if overlaps visible region
                        # Add transparent highlighting to all subplots
                        for ax in [ax_cov, ax_ratio]:
                            ax.axvspan(region['start'], region['end'], 
                                     color='yellow', alpha=0.2, label=region.get('label', 'ROI'))
                        
                        # Add label above the highlighted region
                        if 'label' in region:
                            ax_cov.annotate(region['label'], 
                                          xy=((region['start'] + region['end'])/2, 0.95),
                                          xycoords=('data', 'axes fraction'),
                                          ha='center', va='center',
                                          bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3))
            
            # Set titles and labels for coverage plot
            ax_cov.set_title(f"{os.path.basename(bam)} - {chrom}", fontweight='bold')
            ax_cov.set_ylabel('Mean Read Count (Left Axis)', color='black', fontweight='bold')
            ax_cov.set_ylim(0, y_max)
            
            # Remove x-label from coverage and ratio plots
            ax_cov.set_xlabel('')
            ax_cov.tick_params(axis='x', labelbottom=False)
            ax_ratio.set_xlabel('')
            ax_ratio.tick_params(axis='x', labelbottom=False)
            
            # Add legend to first plot only with better positioning
            if i == 0 and j == 0:
                lines1, labels1 = ax_cov.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines3, labels3 = ax_ratio.get_legend_handles_labels()
                # Use a more visible location for the legend
                leg = ax_cov.legend(lines1 + lines2 + lines3, 
                                  labels1 + labels2 + labels3, 
                                  loc='upper right', 
                                  framealpha=0.9, 
                                  fancybox=True,
                                  shadow=True)
                leg.get_frame().set_edgecolor('#888888')
    
    # Plot gene annotations at the bottom, once per chromosome, after all data plots
    for j, chrom in enumerate(chromosomes):
        if chrom not in chrom_data:
            continue
            
        # Get gene annotation axis - one per chromosome
        ax_genes = axes_genes[j][0] if num_chroms > 1 else axes_genes[0][0]
        
        # Use chromosome bounds for annotations
        min_pos = chrom_data[chrom]['x_min']
        max_pos = chrom_data[chrom]['x_max']
        
        # Add gene annotations if available
        if gtf_features and chrom in gtf_features:
            gene_legend_added = plot_gene_annotations(ax_genes, gtf_features[chrom], min_pos, max_pos, max_genes)
            # Format x-axis in Mbp
            ax_genes.xaxis.set_major_formatter(FuncFormatter(format_kbp))
            ax_genes.set_xlabel('Position (kbp)')
        else:
            ax_genes.text(0.5, 0.5, 'No gene annotations available', 
                         ha='center', va='center', transform=ax_genes.transAxes)
            ax_genes.set_yticks([])
            ax_genes.set_xlabel('Position (bp)')
    
    plt.tight_layout()
    
    # Save the plot if output file is specified
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to {output_file}", file=sys.stderr)
    
    # Display the plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_allele_counts(df, output_file=None, dpi=300, figsize=(12, 8), show_plot=True,
                              min_coverage_threshold=None, highlight_regions=None):
    """
    Generate an improved visualization of allele counts with better visual clarity.
    
    Args:
        df: DataFrame with allele count data
        output_file: Path to save the plot (if None, will display plot)
        dpi: Resolution for the output image
        figsize: Figure size as a tuple (width, height) in inches
        show_plot: Whether to display the plot (in addition to saving it) 
        min_coverage_threshold: Optional minimum coverage threshold to display as horizontal line
        highlight_regions: Optional list of dict with 'start', 'end', and 'label' keys for regions to highlight
    """
    # Check if the required columns exist
    required_cols = ['CHROM', 'Window_start', 'Window_end', 'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads']
    if not all(col in df.columns for col in required_cols):
        print(f"DataFrame must contain columns: {', '.join(required_cols)}", file=sys.stderr)
        return
    
    # Get unique BAM files and chromosomes
    bam_files = df['Bam_file'].unique()
    chromosomes = df['CHROM'].unique()
    
    # Set up the plot
    num_bams = len(bam_files)
    num_chroms = len(chromosomes)
    
    # Create a custom color palette with more distinct colors
    palette = {
        'depth': '#1f77b4',  # Blue
        'ref': '#2ca02c',    # Green
        'alt': '#d62728',    # Red
        'ratio': '#9467bd'   # Purple
    }
    
    # Adjust layout to include ratio subplot
    if num_bams == 1 and num_chroms == 1:
        fig = plt.figure(figsize=figsize)
        gs = plt.GridSpec(2, 1, height_ratios=[3, 1], figure=fig)
        
        # Coverage plot
        ax_cov = fig.add_subplot(gs[0])
        # Ratio plot
        ax_ratio = fig.add_subplot(gs[1], sharex=ax_cov)
        
        axes_cov = [[ax_cov]]
        axes_ratio = [[ax_ratio]]
    elif num_bams == 1:
        fig = plt.figure(figsize=(figsize[0], figsize[1]*num_chroms))
        gs = plt.GridSpec(2*num_chroms, 1, height_ratios=[3, 1] * num_chroms, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for j in range(num_chroms):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[2*j])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[2*j+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
    elif num_chroms == 1:
        fig = plt.figure(figsize=(figsize[0], figsize[1]*num_bams))
        gs = plt.GridSpec(2*num_bams, 1, height_ratios=[3, 1] * num_bams, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[2*i])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[2*i+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
    else:
        fig = plt.figure(figsize=(figsize[0]*num_chroms/2, figsize[1]*num_bams))
        gs = plt.GridSpec(2*num_bams, num_chroms, height_ratios=[3, 1] * num_bams, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            axes_cov_row = []
            axes_ratio_row = []
            
            for j in range(num_chroms):
                # Coverage plot
                ax_cov = fig.add_subplot(gs[2*i, j])
                # Ratio plot
                ax_ratio = fig.add_subplot(gs[2*i+1, j], sharex=ax_cov)
                
                axes_cov_row.append(ax_cov)
                axes_ratio_row.append(ax_ratio)
            
            axes_cov.append(axes_cov_row)
            axes_ratio.append(axes_ratio_row)
    
    # Format function for x-axis (millions of base pairs)
    def format_kbp(x, pos):
        return f'{x/1_000:.2f}'
    
    # Create each subplot
    for i, bam in enumerate(bam_files):
        for j, chrom in enumerate(chromosomes):
            subset = df[(df['Bam_file'] == bam) & (df['CHROM'] == chrom)]
            
            if subset.empty:
                continue
                
            # Calculate window midpoints for x-axis
            subset['window_mid'] = (subset['Window_start'] + subset['Window_end']) / 2
            
            # Calculate alt/total ratio
            subset['alt_ratio'] = (subset['Mean_alt_reads'] / (subset['Mean_depth'].replace(0, np.nan))) * 100
            subset['alt_ratio'] = subset['alt_ratio'].fillna(0)  # Replace NaNs with 0
            
            # Get current axes
            ax_cov = axes_cov[i][j] if num_bams > 1 or num_chroms > 1 else axes_cov[0][0]
            ax_ratio = axes_ratio[i][j] if num_bams > 1 or num_chroms > 1 else axes_ratio[0][0]
            
            # Get axis boundaries for highlighting
            x_min = subset['Window_start'].min()
            x_max = subset['Window_end'].max()
            
            # Plot coverage data with enhanced styling
            # Add shaded background for alt reads area
            ax_cov.set_facecolor('#f9f9f9')  # Light background
            ax_cov.axvspan(x_min, x_max, color='#fff3f3', alpha=0.5)  # Light red shading for alt reads area
            
            # Add grid for better readability
            ax_cov.grid(True, linestyle='--', alpha=0.3)
            
            # Draw coverage threshold if specified
            if min_coverage_threshold is not None:
                ax_cov.axhline(y=min_coverage_threshold, color='#aaaaaa', linestyle='--', alpha=0.7, 
                              label=f'Min. Coverage ({min_coverage_threshold}x)')
            
            # Plot main coverage data with both lines and markers
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], '-', 
                       color=palette['depth'], label='Mean Depth', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], 'o', 
                       color=palette['depth'], markersize=4, alpha=0.5)
            
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], '-', 
                       color=palette['ref'], label='Ref Reads', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], 's', 
                       color=palette['ref'], markersize=4, alpha=0.5)
            
            # Calculate reasonable y-max (max depth + 10%)
            y_max = max(subset['Mean_depth'].max(), subset['Mean_ref_reads'].max()) * 1.1
            
            # Add alternate allele subplot with different scale and clearer styling
            alt_color = palette['alt']
            ax2 = ax_cov.twinx()
            # Add light red background to make it clear this is a different axis
            ax2.patch.set_alpha(0.0)  # Make background transparent
            
            # Add lines and markers for alt reads
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '-', 
                    color=alt_color, linewidth=1.5, label='Alt Reads')
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '^', 
                    color=alt_color, markersize=4, alpha=0.5)
            
            # Set alt reads scale
            alt_max = max(subset['Mean_alt_reads'].max() * 1.5, 5)  # Ensure non-zero max
            ax2.set_ylim(0, alt_max)
            
            # Styling for alt reads axis
            alt_label = 'Alternate Reads (Right Axis)'
            ax2.set_ylabel(alt_label, color=alt_color, fontweight='bold')
            ax2.tick_params(axis='y', labelcolor=alt_color)
            
            # Draw a border around the right axis label to make it stand out
            ax2.yaxis.label.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor=alt_color, pad=5, boxstyle='round'))
            
            # Add annotation to clarify alt reads is on right axis
            ax_cov.annotate('Alt Reads →', xy=(0.98, 0.5), xycoords='axes fraction',
                          ha='right', va='center', color=alt_color, fontsize=9,
                          bbox=dict(facecolor='white', alpha=0.8, edgecolor=alt_color))
            
            # Plot ratio data
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], '-', 
                        color=palette['ratio'], label='Alt/Total Ratio', linewidth=1.5)
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], 'o', 
                        color=palette['ratio'], markersize=4, alpha=0.5)
            ax_ratio.set_ylim(0, min(100.0, subset['alt_ratio'].max() * 1.5 or 10.0))
            ax_ratio.set_ylabel('Alt/Total Ratio')
            ax_ratio.grid(True, linestyle='--', alpha=0.3)
            
            # Highlight regions of interest if specified
            if highlight_regions:
                for region in highlight_regions:
                    if region['start'] <= x_max and region['end'] >= x_min:  # Only if overlaps visible region
                        # Add transparent highlighting to all subplots
                        for ax in [ax_cov, ax_ratio]:
                            ax.axvspan(region['start'], region['end'], 
                                     color='yellow', alpha=0.2, label=region.get('label', 'ROI'))
                        
                        # Add label above the highlighted region
                        if 'label' in region:
                            ax_cov.annotate(region['label'], 
                                          xy=((region['start'] + region['end'])/2, 0.95),
                                          xycoords=('data', 'axes fraction'),
                                          ha='center', va='center',
                                          bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3))
            
            # Set titles and labels for coverage plot
            ax_cov.set_title(f"{os.path.basename(bam)} - {chrom}", fontweight='bold')
            ax_cov.set_ylabel('Mean Read Count (Left Axis)', color='black', fontweight='bold')
            ax_cov.set_ylim(0, y_max)
            
            # Remove x-label from coverage plot
            ax_cov.set_xlabel('')
            ax_cov.tick_params(axis='x', labelbottom=False)
            
            # Format x-axis in Mbp for ratio plot
            ax_ratio.xaxis.set_major_formatter(FuncFormatter(format_kbp))
            ax_ratio.set_xlabel('Position (kbp)')
            
            # Add legend to first plot only with better positioning
            if i == 0 and j == 0:
                lines1, labels1 = ax_cov.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines3, labels3 = ax_ratio.get_legend_handles_labels()
                # Use a more visible location for the legend
                leg = ax_cov.legend(lines1 + lines2 + lines3, 
                                  labels1 + labels2 + labels3, 
                                  loc='upper right', 
                                  framealpha=0.9, 
                                  fancybox=True,
                                  shadow=True)
                leg.get_frame().set_edgecolor('#888888')
    
    plt.tight_layout()
    
    # Save the plot if output file is specified
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to {output_file}", file=sys.stderr)
    
    # Display the plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()