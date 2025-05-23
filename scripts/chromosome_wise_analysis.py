#!/usr/bin/env python3
"""
Chromosome-wise OpenCRAVAT Analysis
==================================
Analyzes all 1.38M variants by processing each chromosome separately
Creates organized output folders for each chromosome
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import time
import os
import gc
from pathlib import Path

warnings.filterwarnings('ignore')
plt.ioff()

def create_output_structure():
    """Create organized folder structure for outputs"""
    base_dir = Path("opencravat_analysis_results")
    base_dir.mkdir(exist_ok=True)
    
    # Create main folders
    folders = {
        'summary': base_dir / "00_overall_summary",
        'chromosomes': base_dir / "chromosomes"
    }
    
    for folder in folders.values():
        folder.mkdir(exist_ok=True)
    
    return folders

def load_data_by_chromosome(filename):
    """Load all data and group by chromosome"""
    print("Loading complete OpenCRAVAT dataset...")
    print("This will process all 1.38 million variants")
    
    start_time = time.time()
    
    try:
        # Load in chunks and group by chromosome
        print("Reading data in chunks...")
        chromosome_data = {}
        chunk_count = 0
        total_variants = 0
        
        chunk_size = 50000  # Smaller chunks for better memory management
        
        for chunk in pd.read_csv(filename, comment='#', chunksize=chunk_size, low_memory=False):
            chunk_count += 1
            total_variants += len(chunk)
            print(f"  Processing chunk {chunk_count}: {len(chunk):,} variants (Total: {total_variants:,})")
            
            # Group this chunk by chromosome
            if 'chrom' in chunk.columns:
                for chrom in chunk['chrom'].unique():
                    chrom_variants = chunk[chunk['chrom'] == chrom].copy()
                    
                    if chrom not in chromosome_data:
                        chromosome_data[chrom] = []
                    
                    chromosome_data[chrom].append(chrom_variants)
            
            # Clear memory
            del chunk
            gc.collect()
        
        print(f"\nData loading completed!")
        print(f"Total variants processed: {total_variants:,}")
        print(f"Loading time: {time.time() - start_time:.1f} seconds")
        print(f"Chromosomes found: {list(chromosome_data.keys())}")
        
        return chromosome_data
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def analyze_single_chromosome(chrom_name, chrom_chunks, output_dir):
    """Analyze a single chromosome's data"""
    print(f"\n{'='*60}")
    print(f"ANALYZING CHROMOSOME {chrom_name}")
    print(f"{'='*60}")
    
    # Combine all chunks for this chromosome
    print(f"Combining {len(chrom_chunks)} chunks for chromosome {chrom_name}...")
    df = pd.concat(chrom_chunks, ignore_index=True)
    
    # Clear chunks from memory
    del chrom_chunks
    gc.collect()
    
    print(f"Chromosome {chrom_name}: {len(df):,} variants")
    
    # Create chromosome-specific output directory
    chrom_dir = output_dir / f"chr_{chrom_name.replace('chr', '')}"
    chrom_dir.mkdir(exist_ok=True)
    
    # Basic statistics
    stats = {
        'chromosome': chrom_name,
        'total_variants': len(df),
        'unique_genes': 0,
        'clinvar_annotated': 0,
        'pathogenic_variants': 0,
        'high_impact_cadd': 0,
        'ultra_rare_variants': 0
    }
    
    # Gene analysis
    if 'hugo' in df.columns:
        gene_data = df['hugo'].dropna()
        stats['unique_genes'] = gene_data.nunique()
        
        if len(gene_data) > 0:
            gene_counts = gene_data.value_counts()
            
            # Save top genes
            top_genes = gene_counts.head(20)
            
            plt.figure(figsize=(12, 8))
            top_genes.plot(kind='barh', color='lightgreen', edgecolor='black')
            plt.title(f'Top 20 Genes - Chromosome {chrom_name}', fontsize=16)
            plt.xlabel('Number of Variants')
            plt.tight_layout()
            plt.savefig(chrom_dir / f'top_genes_chr_{chrom_name.replace("chr", "")}.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save gene counts to CSV
            gene_counts.to_csv(chrom_dir / f'gene_counts_chr_{chrom_name.replace("chr", "")}.csv')
            
            print(f"  Unique genes: {stats['unique_genes']:,}")
            print(f"  Top gene: {top_genes.index[0]} ({top_genes.iloc[0]:,} variants)")
    
    # Clinical significance analysis
    if 'clinvar.sig' in df.columns:
        clin_data = df['clinvar.sig'].dropna()
        stats['clinvar_annotated'] = len(clin_data)
        
        if len(clin_data) > 0:
            clin_counts = clin_data.value_counts()
            
            # Count pathogenic variants
            pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
            stats['pathogenic_variants'] = sum(clin_counts.get(term, 0) for term in pathogenic_terms)
            
            # Save clinical significance plot
            plt.figure(figsize=(12, 6))
            top_clin = clin_counts.head(10)
            top_clin.plot(kind='bar', color='orange', edgecolor='black')
            plt.title(f'Clinical Significance - Chromosome {chrom_name}', fontsize=16)
            plt.xlabel('Clinical Significance')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(chrom_dir / f'clinical_significance_chr_{chrom_name.replace("chr", "")}.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save clinical data to CSV
            clin_counts.to_csv(chrom_dir / f'clinical_significance_chr_{chrom_name.replace("chr", "")}.csv')
            
            print(f"  ClinVar annotated: {stats['clinvar_annotated']:,}")
            print(f"  Pathogenic variants: {stats['pathogenic_variants']:,}")
    
    # Frequency analysis
    if 'gnomad.af' in df.columns:
        freq_data = pd.to_numeric(df['gnomad.af'], errors='coerce').dropna()
        
        if len(freq_data) > 0:
            stats['ultra_rare_variants'] = (freq_data < 0.0001).sum()
            
            # Save frequency distribution plot
            plt.figure(figsize=(12, 6))
            plt.hist(freq_data, bins=50, alpha=0.7, color='purple', edgecolor='black')
            plt.title(f'Allele Frequency Distribution - Chromosome {chrom_name}', fontsize=16)
            plt.xlabel('Allele Frequency')
            plt.ylabel('Count')
            plt.yscale('log')
            plt.tight_layout()
            plt.savefig(chrom_dir / f'allele_frequency_chr_{chrom_name.replace("chr", "")}.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Variants with frequency data: {len(freq_data):,}")
            print(f"  Ultra-rare variants: {stats['ultra_rare_variants']:,}")
    
    # CADD scores analysis
    if 'cadd.phred' in df.columns:
        cadd_data = pd.to_numeric(df['cadd.phred'], errors='coerce').dropna()
        
        if len(cadd_data) > 0:
            stats['high_impact_cadd'] = (cadd_data > 20).sum()
            
            # Save CADD distribution plot
            plt.figure(figsize=(12, 6))
            plt.hist(cadd_data, bins=50, alpha=0.7, color='gold', edgecolor='black')
            plt.title(f'CADD Score Distribution - Chromosome {chrom_name}', fontsize=16)
            plt.xlabel('CADD Score')
            plt.ylabel('Count')
            plt.axvline(x=20, color='red', linestyle='--', label='High impact (>20)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(chrom_dir / f'cadd_scores_chr_{chrom_name.replace("chr", "")}.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  High impact CADD variants: {stats['high_impact_cadd']:,}")
    
    # Variant types analysis
    if 'so' in df.columns:
        so_data = df['so'].dropna()
        
        if len(so_data) > 0:
            so_counts = so_data.value_counts()
            
            # Save variant types plot
            plt.figure(figsize=(12, 8))
            top_so = so_counts.head(15)
            top_so.plot(kind='barh', color='lightcoral', edgecolor='black')
            plt.title(f'Top 15 Variant Types - Chromosome {chrom_name}', fontsize=16)
            plt.xlabel('Count')
            plt.tight_layout()
            plt.savefig(chrom_dir / f'variant_types_chr_{chrom_name.replace("chr", "")}.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save variant types to CSV
            so_counts.to_csv(chrom_dir / f'variant_types_chr_{chrom_name.replace("chr", "")}.csv')
            
            print(f"  Most common variant type: {so_counts.index[0]} ({so_counts.iloc[0]:,})")
    
    # Create chromosome summary dashboard
    create_chromosome_dashboard(df, chrom_name, chrom_dir)
    
    # Save chromosome statistics
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(chrom_dir / f'summary_stats_chr_{chrom_name.replace("chr", "")}.csv', index=False)
    
    # Clear memory
    del df
    gc.collect()
    
    return stats

def create_chromosome_dashboard(df, chrom_name, output_dir):
    """Create a comprehensive dashboard for a single chromosome"""
    fig = plt.figure(figsize=(20, 12))
    
    # 1. Gene distribution
    if 'hugo' in df.columns:
        ax1 = plt.subplot(2, 3, 1)
        gene_counts = df['hugo'].dropna().value_counts().head(10)
        if len(gene_counts) > 0:
            gene_counts.plot(kind='barh', ax=ax1, color='lightgreen')
            ax1.set_title(f'Top 10 Genes - Chr {chrom_name}')
    
    # 2. Clinical significance
    if 'clinvar.sig' in df.columns:
        ax2 = plt.subplot(2, 3, 2)
        clin_data = df['clinvar.sig'].dropna().value_counts().head(8)
        if len(clin_data) > 0:
            clin_data.plot(kind='bar', ax=ax2, color='orange')
            ax2.set_title(f'Clinical Significance - Chr {chrom_name}')
            ax2.tick_params(axis='x', rotation=45)
    
    # 3. Allele frequency
    if 'gnomad.af' in df.columns:
        ax3 = plt.subplot(2, 3, 3)
        freq_data = pd.to_numeric(df['gnomad.af'], errors='coerce').dropna()
        if len(freq_data) > 0:
            ax3.hist(freq_data, bins=30, alpha=0.7, color='purple')
            ax3.set_title(f'Allele Frequency - Chr {chrom_name}')
            ax3.set_yscale('log')
    
    # 4. Variant types
    if 'so' in df.columns:
        ax4 = plt.subplot(2, 3, 4)
        so_counts = df['so'].dropna().value_counts().head(8)
        if len(so_counts) > 0:
            so_counts.plot(kind='bar', ax=ax4, color='lightcoral')
            ax4.set_title(f'Variant Types - Chr {chrom_name}')
            ax4.tick_params(axis='x', rotation=45)
    
    # 5. CADD scores
    if 'cadd.phred' in df.columns:
        ax5 = plt.subplot(2, 3, 5)
        cadd_data = pd.to_numeric(df['cadd.phred'], errors='coerce').dropna()
        if len(cadd_data) > 0:
            ax5.hist(cadd_data, bins=30, alpha=0.7, color='gold')
            ax5.set_title(f'CADD Scores - Chr {chrom_name}')
            ax5.axvline(x=20, color='red', linestyle='--')
    
    # 6. Summary stats
    ax6 = plt.subplot(2, 3, 6)
    summary_text = f"CHROMOSOME {chrom_name} SUMMARY\n\n"
    summary_text += f"Total Variants: {len(df):,}\n"
    
    if 'hugo' in df.columns:
        unique_genes = df['hugo'].dropna().nunique()
        summary_text += f"Unique Genes: {unique_genes:,}\n"
    
    if 'clinvar.sig' in df.columns:
        clin_count = df['clinvar.sig'].dropna().shape[0]
        summary_text += f"ClinVar Annotated: {clin_count:,}\n"
    
    if 'cadd.phred' in df.columns:
        cadd_data = pd.to_numeric(df['cadd.phred'], errors='coerce').dropna()
        high_impact = (cadd_data > 20).sum()
        summary_text += f"High Impact: {high_impact:,}\n"
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, 
             fontsize=12, verticalalignment='top', fontfamily='monospace')
    ax6.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / f'dashboard_chr_{chrom_name.replace("chr", "")}.png', 
               dpi=300, bbox_inches='tight')
    plt.close()

def create_overall_summary(all_stats, output_dir):
    """Create overall summary across all chromosomes"""
    print(f"\n{'='*60}")
    print("CREATING OVERALL SUMMARY")
    print(f"{'='*60}")
    
    # Convert stats to DataFrame
    summary_df = pd.DataFrame(all_stats)
    
    # Sort chromosomes properly
    def sort_chrom(chrom_str):
        try: return int(chrom_str.replace('chr', ''))
        except: return 23 if 'X' in chrom_str else 24 if 'Y' in chrom_str else 25
    
    summary_df['sort_key'] = summary_df['chromosome'].apply(sort_chrom)
    summary_df = summary_df.sort_values('sort_key').drop('sort_key', axis=1)
    
    # Save summary statistics
    summary_df.to_csv(output_dir / 'overall_summary_statistics.csv', index=False)
    
    # Create overall summary plots
    fig = plt.figure(figsize=(20, 12))
    
    # 1. Variants per chromosome
    ax1 = plt.subplot(2, 3, 1)
    summary_df.set_index('chromosome')['total_variants'].plot(kind='bar', ax=ax1, color='skyblue')
    ax1.set_title('Variants per Chromosome')
    ax1.set_ylabel('Number of Variants')
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. Genes per chromosome
    ax2 = plt.subplot(2, 3, 2)
    summary_df.set_index('chromosome')['unique_genes'].plot(kind='bar', ax=ax2, color='lightgreen')
    ax2.set_title('Unique Genes per Chromosome')
    ax2.set_ylabel('Number of Genes')
    ax2.tick_params(axis='x', rotation=45)
    
    # 3. Clinical annotations
    ax3 = plt.subplot(2, 3, 3)
    summary_df.set_index('chromosome')['clinvar_annotated'].plot(kind='bar', ax=ax3, color='orange')
    ax3.set_title('ClinVar Annotations per Chromosome')
    ax3.set_ylabel('Number of Annotated Variants')
    ax3.tick_params(axis='x', rotation=45)
    
    # 4. Pathogenic variants
    ax4 = plt.subplot(2, 3, 4)
    summary_df.set_index('chromosome')['pathogenic_variants'].plot(kind='bar', ax=ax4, color='red')
    ax4.set_title('Pathogenic Variants per Chromosome')
    ax4.set_ylabel('Number of Pathogenic Variants')
    ax4.tick_params(axis='x', rotation=45)
    
    # 5. High impact variants
    ax5 = plt.subplot(2, 3, 5)
    summary_df.set_index('chromosome')['high_impact_cadd'].plot(kind='bar', ax=ax5, color='gold')
    ax5.set_title('High Impact Variants per Chromosome (CADD>20)')
    ax5.set_ylabel('Number of High Impact Variants')
    ax5.tick_params(axis='x', rotation=45)
    
    # 6. Overall statistics
    ax6 = plt.subplot(2, 3, 6)
    total_variants = summary_df['total_variants'].sum()
    total_genes = summary_df['unique_genes'].sum()
    total_pathogenic = summary_df['pathogenic_variants'].sum()
    total_high_impact = summary_df['high_impact_cadd'].sum()
    
    overall_text = f"OVERALL SUMMARY\n\n"
    overall_text += f"Total Variants: {total_variants:,}\n"
    overall_text += f"Total Unique Genes: {total_genes:,}\n"
    overall_text += f"Total Pathogenic: {total_pathogenic:,}\n"
    overall_text += f"Total High Impact: {total_high_impact:,}\n"
    overall_text += f"Chromosomes Analyzed: {len(summary_df)}\n"
    
    ax6.text(0.1, 0.9, overall_text, transform=ax6.transAxes, 
             fontsize=14, verticalalignment='top', fontfamily='monospace')
    ax6.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'overall_summary_dashboard.png', 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print summary
    print(f"FINAL SUMMARY:")
    print(f"- Total variants analyzed: {total_variants:,}")
    print(f"- Chromosomes processed: {len(summary_df)}")
    print(f"- Total pathogenic variants: {total_pathogenic:,}")
    print(f"- Total high impact variants: {total_high_impact:,}")
    
    return summary_df

def main():
    """Main analysis function"""
    filename = "your_opencravat_data.csv"  # Replace with your OpenCRAVAT annotated file
    
    print("Chromosome-wise Analysis")
    print("=" * 60)
    print("This will analyze ALL 1.38 million variants")
    print("Results will be organized by chromosome")
    
    start_time = time.time()
    
    # Create output structure
    folders = create_output_structure()
    print(f"Output will be saved to: {folders['summary'].parent}")
    
    # Load data grouped by chromosome
    chromosome_data = load_data_by_chromosome(filename)
    
    if chromosome_data is None:
        print("Failed to load data")
        return
    
    # Analyze each chromosome
    all_stats = []
    
    for chrom_name in sorted(chromosome_data.keys(), 
                           key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() 
                           else 23 if 'X' in x else 24 if 'Y' in x else 25):
        
        chrom_chunks = chromosome_data[chrom_name]
        stats = analyze_single_chromosome(chrom_name, chrom_chunks, folders['chromosomes'])
        all_stats.append(stats)
        
        print(f"✓ Completed analysis for chromosome {chrom_name}")
    
    # Create overall summary
    summary_df = create_overall_summary(all_stats, folders['summary'])
    
    # Final timing
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print("CHROMOSOME-WISE ANALYSIS COMPLETED!")
    print(f"{'='*60}")
    print(f"Total analysis time: {total_time/60:.1f} minutes")
    print(f"Results saved in: {folders['summary'].parent}")
    
    # Clean up
    del chromosome_data
    gc.collect()

if __name__ == "__main__":
    main() 