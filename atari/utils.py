from typing import List, Dict, Any
import os

def parse_bed_file(bed_file: str) -> List[Dict[str, any]]:
    """Parse regions from a BED file."""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
                
            regions.append({
                'chrom': fields[0],
                'start': int(fields[1]) + 1,  # Convert from 0-based to 1-based
                'end': int(fields[2])
            })
    
    return regions


def get_total_regions_size(regions: List[Dict[str, any]]) -> int:
    """Calculate total size of all regions."""
    total_size = 0
    for region in regions:
        total_size += region['end'] - region['start'] + 1
    return total_size


def parse_gtf_file(gtf_file: str) -> Dict[str, List[Dict]]:
    """Parse gene and exon information from a GTF file.
    
    Args:
        gtf_file: Path to the GTF file
        
    Returns:
        Dictionary mapping chromosome names to lists of features
    """
    features_by_chrom = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Parse GTF line
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])  # GTF is 1-based
            end = int(fields[4])
            strand = fields[6]
            
            # Only process genes and exons
            if feature_type not in ['gene', 'exon', 'transcript', 'CDS']:
                continue
            
            # Parse attributes
            attr_str = fields[8]
            attrs = {}
            
            # Extract key attributes
            for attr_pair in attr_str.split(';'):
                attr_pair = attr_pair.strip()
                if not attr_pair:
                    continue
                
                # Handle different GTF formats (Ensembl vs UCSC)
                if '=' in attr_pair:
                    key, value = attr_pair.split('=', 1)
                    value = value.strip('"')
                elif ' ' in attr_pair:
                    key, value = attr_pair.split(' ', 1)
                    value = value.strip('"')
                else:
                    continue
                
                attrs[key] = value
            
            # Get gene name and ID
            gene_name = attrs.get('gene_name', attrs.get('gene_id', 'Unknown'))
            gene_id = attrs.get('gene_id', 'Unknown')
            
            # Store feature information
            feature = {
                'type': feature_type,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_name': gene_name,
                'gene_id': gene_id
            }
            
            # Add transcript and exon IDs if available
            if feature_type in ['transcript', 'exon', 'CDS']:
                feature['transcript_id'] = attrs.get('transcript_id', 'Unknown')
                
            # Add to chromosome dictionary
            if chrom not in features_by_chrom:
                features_by_chrom[chrom] = []
            
            features_by_chrom[chrom].append(feature)
    
    return features_by_chrom