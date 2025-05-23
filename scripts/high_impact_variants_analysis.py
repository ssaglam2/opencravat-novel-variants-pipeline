#!/usr/bin/env python3
"""
High Impact Variants Analysis
============================
Extracts detailed information about high impact variants including:
- MAF (Minor Allele Frequency)
- Deleteriousness scores (CADD, REVEL, PolyPhen2, SIFT)
- Gene names and variant details
- Clinical significance
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import time

warnings.filterwarnings('ignore')

def load_and_filter_high_impact_variants(filename):
    """Load data and filter for high impact variants"""
    print("Loading OpenCRAVAT data to extract high impact variants...")
    print("This will process all variants and filter for high impact ones")
    
    start_time = time.time()
    
    # Load data in chunks to manage memory
    high_impact_variants = []
    chunk_count = 0
    total_variants = 0
    high_impact_count = 0
    
    chunk_size = 100000
    
    for chunk in pd.read_csv(filename, comment='#', chunksize=chunk_size, low_memory=False):
        chunk_count += 1
        total_variants += len(chunk)
        print(f"Processing chunk {chunk_count}: {len(chunk):,} variants (Total: {total_variants:,})")
        
        # Filter for high impact variants based on multiple criteria
        high_impact_chunk = filter_high_impact_chunk(chunk)
        
        if len(high_impact_chunk) > 0:
            high_impact_variants.append(high_impact_chunk)
            high_impact_count += len(high_impact_chunk)
            print(f"  Found {len(high_impact_chunk):,} high impact variants in this chunk")
    
    if high_impact_variants:
        # Combine all high impact variants
        df_high_impact = pd.concat(high_impact_variants, ignore_index=True)
        print(f"\nTotal high impact variants found: {len(df_high_impact):,}")
        print(f"Loading time: {time.time() - start_time:.1f} seconds")
        
        return df_high_impact
    else:
        print("No high impact variants found")
        return pd.DataFrame()

def filter_high_impact_chunk(chunk):
    """Filter chunk for high impact variants"""
    high_impact_conditions = []
    
    # 1. High CADD score (>20)
    if 'cadd.phred' in chunk.columns:
        cadd_data = pd.to_numeric(chunk['cadd.phred'], errors='coerce')
        high_impact_conditions.append(cadd_data > 20)
    
    # 2. Pathogenic/Likely pathogenic clinical significance
    if 'clinvar.sig' in chunk.columns:
        pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
        pathogenic_condition = chunk['clinvar.sig'].isin(pathogenic_terms)
        high_impact_conditions.append(pathogenic_condition)
    
    # 3. High REVEL score (>0.5)
    if 'revel.score' in chunk.columns:
        revel_data = pd.to_numeric(chunk['revel.score'], errors='coerce')
        high_impact_conditions.append(revel_data > 0.5)
    
    # 4. Loss of function variants
    if 'so' in chunk.columns:
        lof_terms = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 
                     'splice_donor_variant', 'start_lost', 'stop_lost']
        lof_condition = chunk['so'].isin(lof_terms)
        high_impact_conditions.append(lof_condition)
    
    # Combine conditions (any one of them qualifies as high impact)
    if high_impact_conditions:
        combined_condition = high_impact_conditions[0]
        for condition in high_impact_conditions[1:]:
            combined_condition = combined_condition | condition
        
        return chunk[combined_condition].copy()
    else:
        return pd.DataFrame()

def analyze_high_impact_variants(df):
    """Analyze high impact variants in detail"""
    print(f"\n{'='*80}")
    print("HIGH IMPACT VARIANTS DETAILED ANALYSIS")
    print(f"{'='*80}")
    
    print(f"Total high impact variants: {len(df):,}")
    
    # Select key columns for analysis
    key_columns = [
        'chrom', 'pos', 'ref', 'alt', 'hugo', 'so', 'hgvs.p', 'hgvs.c',
        'cadd.phred', 'revel.score', 'polyphen2.hdiv_pred', 'sift.pred',
        'clinvar.sig', 'gnomad.af', 'thousandgenomes.af', 'exac.af'
    ]
    
    # Keep only columns that exist in the data
    available_columns = [col for col in key_columns if col in df.columns]
    df_analysis = df[available_columns].copy()
    
    print(f"Available columns for analysis: {len(available_columns)}")
    
    # 1. Chromosome distribution
    print(f"\n{'='*60}")
    print("CHROMOSOME DISTRIBUTION OF HIGH IMPACT VARIANTS")
    print(f"{'='*60}")
    
    if 'chrom' in df_analysis.columns:
        chrom_counts = df_analysis['chrom'].value_counts()
        print("Top 20 chromosomes with high impact variants:")
        for i, (chrom, count) in enumerate(chrom_counts.head(20).items(), 1):
            percentage = (count / len(df_analysis)) * 100
            print(f"{i:2d}. {chrom}: {count:,} ({percentage:.1f}%)")
    
    # 2. Gene analysis
    print(f"\n{'='*60}")
    print("GENES WITH HIGH IMPACT VARIANTS")
    print(f"{'='*60}")
    
    if 'hugo' in df_analysis.columns:
        gene_data = df_analysis['hugo'].dropna()
        if len(gene_data) > 0:
            gene_counts = gene_data.value_counts()
            print(f"Genes with high impact variants: {gene_counts.nunique():,}")
            print("\nTop 20 genes with most high impact variants:")
            for i, (gene, count) in enumerate(gene_counts.head(20).items(), 1):
                print(f"{i:2d}. {gene}: {count:,} variants")
    
    # 3. Variant types
    print(f"\n{'='*60}")
    print("VARIANT TYPES (SEQUENCE ONTOLOGY)")
    print(f"{'='*60}")
    
    if 'so' in df_analysis.columns:
        so_data = df_analysis['so'].dropna()
        if len(so_data) > 0:
            so_counts = so_data.value_counts()
            print("Distribution of high impact variant types:")
            for i, (so_term, count) in enumerate(so_counts.head(15).items(), 1):
                percentage = (count / len(so_data)) * 100
                print(f"{i:2d}. {so_term}: {count:,} ({percentage:.1f}%)")
    
    # 4. Clinical significance
    print(f"\n{'='*60}")
    print("CLINICAL SIGNIFICANCE")
    print(f"{'='*60}")
    
    if 'clinvar.sig' in df_analysis.columns:
        clin_data = df_analysis['clinvar.sig'].dropna()
        if len(clin_data) > 0:
            clin_counts = clin_data.value_counts()
            print(f"Variants with clinical annotations: {len(clin_data):,}")
            print("Clinical significance distribution:")
            for sig, count in clin_counts.items():
                percentage = (count / len(clin_data)) * 100
                print(f"  {sig}: {count:,} ({percentage:.1f}%)")
    
    # 5. Minor Allele Frequency Analysis
    print(f"\n{'='*60}")
    print("MINOR ALLELE FREQUENCY (MAF) ANALYSIS")
    print(f"{'='*60}")
    
    freq_columns = ['gnomad.af', 'thousandgenomes.af', 'exac.af']
    available_freq_cols = [col for col in freq_columns if col in df_analysis.columns]
    
    for freq_col in available_freq_cols:
        print(f"\n{freq_col.upper()} Analysis:")
        freq_data = pd.to_numeric(df_analysis[freq_col], errors='coerce').dropna()
        
        if len(freq_data) > 0:
            print(f"  Variants with {freq_col} data: {len(freq_data):,}")
            print(f"  Mean frequency: {freq_data.mean():.6f}")
            print(f"  Median frequency: {freq_data.median():.6f}")
            
            # Frequency categories
            ultra_rare = (freq_data < 0.0001).sum()
            rare = ((freq_data >= 0.0001) & (freq_data < 0.01)).sum()
            common = (freq_data >= 0.01).sum()
            
            print(f"  Ultra-rare (AF < 0.01%): {ultra_rare:,} ({ultra_rare/len(freq_data)*100:.1f}%)")
            print(f"  Rare (0.01% ≤ AF < 1%): {rare:,} ({rare/len(freq_data)*100:.1f}%)")
            print(f"  Common (AF ≥ 1%): {common:,} ({common/len(freq_data)*100:.1f}%)")
    
    # 6. Deleteriousness Scores
    print(f"\n{'='*60}")
    print("DELETERIOUSNESS SCORES")
    print(f"{'='*60}")
    
    # CADD scores
    if 'cadd.phred' in df_analysis.columns:
        cadd_data = pd.to_numeric(df_analysis['cadd.phred'], errors='coerce').dropna()
        if len(cadd_data) > 0:
            print(f"CADD Scores:")
            print(f"  Variants with CADD scores: {len(cadd_data):,}")
            print(f"  Mean CADD: {cadd_data.mean():.2f}")
            print(f"  Median CADD: {cadd_data.median():.2f}")
            print(f"  High impact (>20): {(cadd_data > 20).sum():,} ({(cadd_data > 20).sum()/len(cadd_data)*100:.1f}%)")
            print(f"  Ultra-high impact (>30): {(cadd_data > 30).sum():,} ({(cadd_data > 30).sum()/len(cadd_data)*100:.1f}%)")
    
    # REVEL scores
    if 'revel.score' in df_analysis.columns:
        revel_data = pd.to_numeric(df_analysis['revel.score'], errors='coerce').dropna()
        if len(revel_data) > 0:
            print(f"\nREVEL Scores:")
            print(f"  Variants with REVEL scores: {len(revel_data):,}")
            print(f"  Mean REVEL: {revel_data.mean():.3f}")
            print(f"  Median REVEL: {revel_data.median():.3f}")
            print(f"  Pathogenic (>0.5): {(revel_data > 0.5).sum():,} ({(revel_data > 0.5).sum()/len(revel_data)*100:.1f}%)")
            print(f"  Highly pathogenic (>0.75): {(revel_data > 0.75).sum():,} ({(revel_data > 0.75).sum()/len(revel_data)*100:.1f}%)")
    
    # PolyPhen2 predictions
    if 'polyphen2.hdiv_pred' in df_analysis.columns:
        polyphen_data = df_analysis['polyphen2.hdiv_pred'].dropna()
        if len(polyphen_data) > 0:
            print(f"\nPolyPhen2 Predictions:")
            polyphen_counts = polyphen_data.value_counts()
            for pred, count in polyphen_counts.items():
                percentage = (count / len(polyphen_data)) * 100
                print(f"  {pred}: {count:,} ({percentage:.1f}%)")
    
    # SIFT predictions
    if 'sift.pred' in df_analysis.columns:
        sift_data = df_analysis['sift.pred'].dropna()
        if len(sift_data) > 0:
            print(f"\nSIFT Predictions:")
            sift_counts = sift_data.value_counts()
            for pred, count in sift_counts.items():
                percentage = (count / len(sift_data)) * 100
                print(f"  {pred}: {count:,} ({percentage:.1f}%)")
    
    return df_analysis

def create_high_impact_summary_tables(df_analysis):
    """Create summary tables for high impact variants"""
    print(f"\n{'='*60}")
    print("CREATING SUMMARY TABLES")
    print(f"{'='*60}")
    
    # 1. Most critical variants (pathogenic + high CADD + ultra-rare)
    critical_variants = df_analysis.copy()
    
    # Add criticality criteria
    critical_variants['is_pathogenic'] = False
    critical_variants['is_high_cadd'] = False
    critical_variants['is_ultra_rare'] = False
    
    if 'clinvar.sig' in critical_variants.columns:
        pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
        critical_variants['is_pathogenic'] = critical_variants['clinvar.sig'].isin(pathogenic_terms)
    
    if 'cadd.phred' in critical_variants.columns:
        cadd_scores = pd.to_numeric(critical_variants['cadd.phred'], errors='coerce')
        critical_variants['is_high_cadd'] = cadd_scores > 20
    
    if 'gnomad.af' in critical_variants.columns:
        gnomad_freq = pd.to_numeric(critical_variants['gnomad.af'], errors='coerce')
        critical_variants['is_ultra_rare'] = gnomad_freq < 0.0001
    
    # Calculate criticality score
    critical_variants['criticality_score'] = (
        critical_variants['is_pathogenic'].astype(int) +
        critical_variants['is_high_cadd'].astype(int) +
        critical_variants['is_ultra_rare'].astype(int)
    )
    
    # Sort by criticality
    critical_variants_sorted = critical_variants.sort_values(['criticality_score', 'cadd.phred'], 
                                                           ascending=[False, False])
    
    # Save top critical variants
    top_critical = critical_variants_sorted.head(100)
    top_critical.to_csv('top_100_critical_variants.csv', index=False)
    print("✓ Saved top 100 critical variants to: top_100_critical_variants.csv")
    
    # 2. Pathogenic variants summary
    if 'clinvar.sig' in df_analysis.columns:
        pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
        pathogenic_variants = df_analysis[df_analysis['clinvar.sig'].isin(pathogenic_terms)].copy()
        
        if len(pathogenic_variants) > 0:
            pathogenic_variants.to_csv('all_pathogenic_variants.csv', index=False)
            print(f"✓ Saved {len(pathogenic_variants):,} pathogenic variants to: all_pathogenic_variants.csv")
    
    # 3. Ultra-high impact variants (CADD > 30)
    if 'cadd.phred' in df_analysis.columns:
        cadd_scores = pd.to_numeric(df_analysis['cadd.phred'], errors='coerce')
        ultra_high_cadd = df_analysis[cadd_scores > 30].copy()
        
        if len(ultra_high_cadd) > 0:
            ultra_high_cadd.to_csv('ultra_high_cadd_variants.csv', index=False)
            print(f"✓ Saved {len(ultra_high_cadd):,} ultra-high CADD variants to: ultra_high_cadd_variants.csv")
    
    # 4. Loss of function variants
    if 'so' in df_analysis.columns:
        lof_terms = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 
                     'splice_donor_variant', 'start_lost', 'stop_lost']
        lof_variants = df_analysis[df_analysis['so'].isin(lof_terms)].copy()
        
        if len(lof_variants) > 0:
            lof_variants.to_csv('loss_of_function_variants.csv', index=False)
            print(f"✓ Saved {len(lof_variants):,} loss-of-function variants to: loss_of_function_variants.csv")
    
    # 5. Gene-based summary
    if 'hugo' in df_analysis.columns:
        gene_summary = []
        
        for gene in df_analysis['hugo'].dropna().unique():
            gene_variants = df_analysis[df_analysis['hugo'] == gene]
            
            gene_info = {
                'gene': gene,
                'total_variants': len(gene_variants),
                'pathogenic_variants': 0,
                'high_cadd_variants': 0,
                'ultra_rare_variants': 0,
                'max_cadd_score': np.nan,
                'min_gnomad_af': np.nan
            }
            
            # Count pathogenic
            if 'clinvar.sig' in gene_variants.columns:
                pathogenic_terms = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
                gene_info['pathogenic_variants'] = gene_variants['clinvar.sig'].isin(pathogenic_terms).sum()
            
            # Count high CADD
            if 'cadd.phred' in gene_variants.columns:
                cadd_scores = pd.to_numeric(gene_variants['cadd.phred'], errors='coerce')
                gene_info['high_cadd_variants'] = (cadd_scores > 20).sum()
                gene_info['max_cadd_score'] = cadd_scores.max()
            
            # Count ultra-rare
            if 'gnomad.af' in gene_variants.columns:
                gnomad_freq = pd.to_numeric(gene_variants['gnomad.af'], errors='coerce')
                gene_info['ultra_rare_variants'] = (gnomad_freq < 0.0001).sum()
                gene_info['min_gnomad_af'] = gnomad_freq.min()
            
            gene_summary.append(gene_info)
        
        gene_summary_df = pd.DataFrame(gene_summary)
        gene_summary_df = gene_summary_df.sort_values(['pathogenic_variants', 'high_cadd_variants'], 
                                                     ascending=[False, False])
        gene_summary_df.to_csv('gene_based_high_impact_summary.csv', index=False)
        print(f"✓ Saved gene-based summary to: gene_based_high_impact_summary.csv")
    
    return critical_variants_sorted

def create_visualizations(df_analysis):
    """Create visualizations for high impact variants"""
    print(f"\n{'='*60}")
    print("CREATING VISUALIZATIONS")
    print(f"{'='*60}")
    
    # 1. CADD score distribution
    if 'cadd.phred' in df_analysis.columns:
        plt.figure(figsize=(12, 6))
        cadd_data = pd.to_numeric(df_analysis['cadd.phred'], errors='coerce').dropna()
        if len(cadd_data) > 0:
            plt.hist(cadd_data, bins=50, alpha=0.7, color='red', edgecolor='black')
            plt.axvline(x=20, color='blue', linestyle='--', label='High impact threshold')
            plt.axvline(x=30, color='darkred', linestyle='--', label='Ultra-high impact')
            plt.title('CADD Score Distribution - High Impact Variants', fontsize=16)
            plt.xlabel('CADD Score')
            plt.ylabel('Count')
            plt.legend()
            plt.tight_layout()
            plt.savefig('high_impact_cadd_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ CADD distribution plot saved")
    
    # 2. Frequency distribution
    if 'gnomad.af' in df_analysis.columns:
        plt.figure(figsize=(12, 6))
        freq_data = pd.to_numeric(df_analysis['gnomad.af'], errors='coerce').dropna()
        if len(freq_data) > 0:
            # Remove zero frequencies for log scale
            freq_data_nonzero = freq_data[freq_data > 0]
            if len(freq_data_nonzero) > 0:
                plt.hist(freq_data_nonzero, bins=50, alpha=0.7, color='purple', edgecolor='black')
                plt.axvline(x=0.0001, color='red', linestyle='--', label='Ultra-rare threshold')
                plt.axvline(x=0.01, color='orange', linestyle='--', label='Rare threshold')
                plt.title('gnomAD Allele Frequency Distribution - High Impact Variants', fontsize=16)
                plt.xlabel('Allele Frequency')
                plt.ylabel('Count')
                plt.xscale('log')
                plt.legend()
                plt.tight_layout()
                plt.savefig('high_impact_frequency_distribution.png', dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ Frequency distribution plot saved")
    
    # 3. Top genes with high impact variants
    if 'hugo' in df_analysis.columns:
        gene_counts = df_analysis['hugo'].dropna().value_counts().head(20)
        if len(gene_counts) > 0:
            plt.figure(figsize=(12, 10))
            gene_counts.plot(kind='barh', color='lightblue', edgecolor='black')
            plt.title('Top 20 Genes with High Impact Variants', fontsize=16)
            plt.xlabel('Number of High Impact Variants')
            plt.tight_layout()
            plt.savefig('top_genes_high_impact.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Top genes plot saved")

def main():
    """Main analysis function"""
    filename = "your_opencravat_data.csv"  # Replace with your OpenCRAVAT annotated file
    
    print("High Impact Variants Analysis")
    print("=" * 80)
    print("Extracting detailed information about high impact variants...")
    
    start_time = time.time()
    
    # Load and filter high impact variants
    df_high_impact = load_and_filter_high_impact_variants(filename)
    
    if df_high_impact.empty:
        print("No high impact variants found")
        return
    
    # Analyze variants
    df_analysis = analyze_high_impact_variants(df_high_impact)
    
    # Create summary tables
    critical_variants = create_high_impact_summary_tables(df_analysis)
    
    # Create visualizations
    create_visualizations(df_analysis)
    
    # Save complete high impact dataset
    df_analysis.to_csv('all_high_impact_variants.csv', index=False)
    
    # Final summary
    total_time = time.time() - start_time
    print(f"\n{'='*80}")
    print("HIGH IMPACT VARIANTS ANALYSIS COMPLETED!")
    print(f"{'='*80}")
    print(f"Analysis time: {total_time/60:.1f} minutes")
    print(f"Total high impact variants: {len(df_analysis):,}")
    
    print(f"\nFiles created:")
    print(f"✓ all_high_impact_variants.csv - Complete dataset")
    print(f"✓ top_100_critical_variants.csv - Most critical variants")
    print(f"✓ all_pathogenic_variants.csv - ClinVar pathogenic variants")
    print(f"✓ ultra_high_cadd_variants.csv - CADD > 30 variants")
    print(f"✓ loss_of_function_variants.csv - LoF variants")
    print(f"✓ gene_based_high_impact_summary.csv - Gene-level summary")
    print(f"✓ Visualization plots generated")

if __name__ == "__main__":
    main() 