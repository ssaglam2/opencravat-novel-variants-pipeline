#!/usr/bin/env python3
"""
European Population-Specific Novel Variants Discovery
====================================================
Enhanced novel variant discovery tailored for European populations:
- Uses European-specific frequency data (gnomAD NFE, 1000G European)
- Cross-population validation (rare in Europeans, common elsewhere)
- Population-specific pathogenic variant detection
- Enhanced scoring system for European genetic background
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import time
import os
from pathlib import Path

warnings.filterwarnings('ignore')
plt.style.use('default')

def create_output_structure():
    """Create organized folder structure for European novel variant analysis"""
    base_dir = Path("european_novel_variants_analysis")
    base_dir.mkdir(exist_ok=True)
    
    folders = {
        'data': base_dir / "data_tables",
        'plots': base_dir / "visualizations", 
        'reports': base_dir / "reports",
        'population_analysis': base_dir / "population_analysis"
    }
    
    for folder in folders.values():
        folder.mkdir(exist_ok=True)
    
    return folders

def get_european_frequency_columns():
    """Define European-specific frequency columns in order of priority"""
    return {
        'primary_european': [
            'gnomad.af_nfe',  # Non-Finnish European (most representative)
            'gnomad.af_fin'   # Finnish European
        ],
        'thousand_genomes_european': [
            'thousandgenomes_european.gbr_af',  # British
            'thousandgenomes_european.ibs_af',  # Spanish  
            'thousandgenomes_european.tsi_af',  # Italian
            'thousandgenomes_european.ceu_af'   # Utah/CEPH
        ],
        'other_populations': [
            'gnomad.af_afr',  # African
            'gnomad.af_eas',  # East Asian
            'gnomad.af_sas',  # South Asian
            'gnomad.af_amr'   # Latino/American
        ]
    }

def calculate_european_frequency_stats(row):
    """Calculate comprehensive European frequency statistics"""
    freq_cols = get_european_frequency_columns()
    
    # European frequencies
    european_freqs = []
    for col in freq_cols['primary_european'] + freq_cols['thousand_genomes_european']:
        if col in row.index and pd.notna(row[col]):
            european_freqs.append(row[col])
    
    # Non-European frequencies  
    other_freqs = []
    for col in freq_cols['other_populations']:
        if col in row.index and pd.notna(row[col]):
            other_freqs.append(row[col])
    
    stats = {
        'max_european_freq': max(european_freqs) if european_freqs else np.nan,
        'min_european_freq': min(european_freqs) if european_freqs else np.nan,
        'mean_european_freq': np.mean(european_freqs) if european_freqs else np.nan,
        'max_other_freq': max(other_freqs) if other_freqs else np.nan,
        'european_freq_count': len(european_freqs),
        'other_freq_count': len(other_freqs)
    }
    
    return stats

def classify_european_variant(row):
    """Classify variants based on European population frequencies"""
    freq_stats = calculate_european_frequency_stats(row)
    max_euro_freq = freq_stats['max_european_freq']
    max_other_freq = freq_stats['max_other_freq']
    
    if pd.isna(max_euro_freq):
        return 'unknown_frequency'
    
    # European frequency categories
    if max_euro_freq < 0.00001:    # < 0.001%
        euro_category = 'ultra_rare_european'
    elif max_euro_freq < 0.0001:   # 0.001-0.01%
        euro_category = 'very_rare_european'  
    elif max_euro_freq < 0.01:     # 0.01-1%
        euro_category = 'rare_european'
    else:
        euro_category = 'common_european'
    
    # Population specificity
    if not pd.isna(max_other_freq) and max_euro_freq < 0.0001 and max_other_freq > 0.01:
        return f'{euro_category}_population_specific'
    
    return euro_category

def calculate_population_specificity_score(row):
    """Calculate bonus points for population-specific variants"""
    freq_stats = calculate_european_frequency_stats(row)
    max_euro_freq = freq_stats['max_european_freq']
    max_other_freq = freq_stats['max_other_freq']
    
    if pd.isna(max_euro_freq) or pd.isna(max_other_freq):
        return 0
    
    # High score if rare in Europeans but common elsewhere
    if max_euro_freq < 0.0001:  # Rare in Europeans
        if max_other_freq > 0.01:  # But common elsewhere
            # More bonus for greater frequency difference
            ratio = max_other_freq / (max_euro_freq + 1e-8)
            return min(np.log10(ratio), 10)  # Cap at 10 points
    
    return 0

def calculate_european_novelty_score(row):
    """Calculate European population-specific novelty score"""
    score = 0.0
    
    # 1. CADD score component (0-35 points)
    cadd = pd.to_numeric(row.get('cadd.phred', 0), errors='coerce')
    if not pd.isna(cadd):
        cadd_normalized = np.clip(cadd / 40.0, 0, 1)
        score += cadd_normalized * 35
    
    # 2. REVEL score component (0-25 points)
    revel = pd.to_numeric(row.get('revel.score', 0), errors='coerce') 
    if not pd.isna(revel):
        score += revel * 25
    
    # 3. European frequency component (0-20 points)
    freq_stats = calculate_european_frequency_stats(row)
    max_euro_freq = freq_stats['max_european_freq']
    if not pd.isna(max_euro_freq):
        freq_score = np.clip(-np.log10(max_euro_freq + 1e-8) / 8, 0, 1)
        score += freq_score * 20
    
    # 4. Loss of function bonus (0-10 points)
    so = row.get('so', '')
    lof_terms = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 
                 'splice_donor_variant', 'start_lost', 'stop_lost']
    if so in lof_terms:
        score += 10
    
    # 5. Multiple predictor consensus (0-8 points)
    consensus_score = 0
    polyphen = row.get('polyphen2.hdiv_pred', '')
    if polyphen in ['D', 'P']:  # Damaging or Probably damaging
        consensus_score += 3
    
    revel_consensus = pd.to_numeric(row.get('revel.score', 0), errors='coerce')
    if not pd.isna(revel_consensus) and revel_consensus > 0.5:
        consensus_score += 5
    
    score += consensus_score
    
    # 6. Population specificity bonus (0-10 points)
    pop_spec_score = calculate_population_specificity_score(row)
    score += pop_spec_score
    
    return score

def load_and_filter_european_novel_variants(filename):
    """Load data and filter for European population-specific novel variants"""
    print("Loading OpenCRAVAT data for European novel variant discovery...")
    print("Focus: European population-specific pathogenic variants")
    
    start_time = time.time()
    
    novel_variants = []
    chunk_count = 0
    total_variants = 0
    
    chunk_size = 100000
    
    for chunk in pd.read_csv(filename, comment='#', chunksize=chunk_size, low_memory=False):
        chunk_count += 1
        total_variants += len(chunk)
        print(f"Processing chunk {chunk_count}: {len(chunk):,} variants (Total: {total_variants:,})")
        
        # Filter for European population-specific novel variants
        novel_chunk = filter_european_novel_variants_chunk(chunk)
        
        if len(novel_chunk) > 0:
            novel_variants.append(novel_chunk)
            print(f"  Found {len(novel_chunk):,} potentially novel European variants")
    
    if novel_variants:
        df_novel = pd.concat(novel_variants, ignore_index=True)
        print(f"\nTotal European novel variants: {len(df_novel):,}")
        print(f"Loading time: {time.time() - start_time:.1f} seconds")
        return df_novel
    else:
        print("No European novel variants found")
        return pd.DataFrame()

def filter_european_novel_variants_chunk(chunk):
    """Filter chunk for European population-specific novel variants"""
    conditions = []
    freq_cols = get_european_frequency_columns()
    
    # 1. High CADD score
    if 'cadd.phred' in chunk.columns:
        cadd_data = pd.to_numeric(chunk['cadd.phred'], errors='coerce')
        conditions.append(cadd_data > 25)
    
    # 2. High REVEL score  
    if 'revel.score' in chunk.columns:
        revel_data = pd.to_numeric(chunk['revel.score'], errors='coerce')
        conditions.append(revel_data > 0.7)
    
    # 3. Loss of function variants
    if 'so' in chunk.columns:
        lof_terms = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 
                     'splice_donor_variant', 'start_lost', 'stop_lost']
        conditions.append(chunk['so'].isin(lof_terms))
    
    # 4. Damaging missense with PolyPhen2 and CADD consensus
    if all(col in chunk.columns for col in ['so', 'polyphen2.hdiv_pred', 'cadd.phred']):
        cadd_data = pd.to_numeric(chunk['cadd.phred'], errors='coerce')
        is_missense = chunk['so'] == 'missense_variant'
        is_damaging_polyphen = chunk['polyphen2.hdiv_pred'].isin(['D', 'P'])
        is_high_cadd = cadd_data > 20
        damaging_missense = is_missense & is_damaging_polyphen & is_high_cadd
        conditions.append(damaging_missense)
    
    # 5. Ultra-rare in Europeans with high CADD
    euro_freq_condition = None
    for freq_col in freq_cols['primary_european']:
        if freq_col in chunk.columns:
            euro_freq = pd.to_numeric(chunk[freq_col], errors='coerce')
            cadd_scores = pd.to_numeric(chunk['cadd.phred'], errors='coerce')
            ultra_rare_euro = (euro_freq < 0.00001) & (cadd_scores > 20)
            if euro_freq_condition is None:
                euro_freq_condition = ultra_rare_euro
            else:
                euro_freq_condition = euro_freq_condition | ultra_rare_euro
    
    if euro_freq_condition is not None:
        conditions.append(euro_freq_condition)
    
    # 6. Population-specific variants (rare in Europeans, common elsewhere)
    pop_specific_condition = None
    if any(col in chunk.columns for col in freq_cols['primary_european']) and \
       any(col in chunk.columns for col in freq_cols['other_populations']):
        
        # Find European frequency
        euro_freq = None
        for col in freq_cols['primary_european']:
            if col in chunk.columns:
                euro_freq = pd.to_numeric(chunk[col], errors='coerce')
                break
        
        # Find max other population frequency
        other_freq = None
        for col in freq_cols['other_populations']:
            if col in chunk.columns:
                freq = pd.to_numeric(chunk[col], errors='coerce')
                if other_freq is None:
                    other_freq = freq
                else:
                    other_freq = pd.concat([other_freq, freq], axis=1).max(axis=1)
        
        if euro_freq is not None and other_freq is not None:
            cadd_scores = pd.to_numeric(chunk['cadd.phred'], errors='coerce')
            pop_specific_condition = (euro_freq < 0.0001) & (other_freq > 0.01) & (cadd_scores > 20)
            conditions.append(pop_specific_condition)
    
    # Combine conditions
    if conditions:
        combined_condition = conditions[0]
        for condition in conditions[1:]:
            combined_condition = combined_condition | condition
        
        # Exclude known ClinVar pathogenic
        if 'clinvar.sig' in chunk.columns:
            known_pathogenic = chunk['clinvar.sig'].isin(['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'])
            combined_condition = combined_condition & ~known_pathogenic
        
        return chunk[combined_condition].copy()
    else:
        return pd.DataFrame()

def analyze_european_novel_variants(df, folders):
    """Comprehensive analysis of European novel variants"""
    print(f"\n{'='*80}")
    print("EUROPEAN POPULATION-SPECIFIC NOVEL VARIANTS ANALYSIS")
    print(f"{'='*80}")
    
    print(f"Total European novel variants: {len(df):,}")
    
    # Add European-specific analysis columns
    print("Calculating European population-specific scores...")
    
    # Apply European frequency classification
    european_classifications = []
    novelty_scores = []
    population_scores = []
    
    for idx, row in df.iterrows():
        classification = classify_european_variant(row)
        european_classifications.append(classification)
        
        novelty_score = calculate_european_novelty_score(row)
        novelty_scores.append(novelty_score)
        
        pop_score = calculate_population_specificity_score(row)
        population_scores.append(pop_score)
    
    df['european_classification'] = european_classifications
    df['european_novelty_score'] = novelty_scores
    df['population_specificity_score'] = population_scores
    
    # Sort by European novelty score
    df_sorted = df.sort_values('european_novelty_score', ascending=False)
    
    print(f"European novelty score range: {min(novelty_scores):.1f} - {max(novelty_scores):.1f}")
    
    # Key columns for analysis
    key_columns = [
        'chrom', 'pos', 'ref', 'alt', 'hugo', 'so', 'hgvs.p', 'hgvs.c',
        'cadd.phred', 'revel.score', 'polyphen2.hdiv_pred', 'clinvar.sig',
        'gnomad.af_nfe', 'gnomad.af_fin', 'gnomad.af_afr', 'gnomad.af_eas',
        'thousandgenomes_european.gbr_af', 'thousandgenomes_european.ibs_af',
        'european_classification', 'european_novelty_score', 'population_specificity_score'
    ]
    
    available_columns = [col for col in key_columns if col in df.columns]
    df_analysis = df_sorted[available_columns].copy()
    
    # 1. Top European novel variants
    print(f"\n{'='*60}")
    print("TOP 20 EUROPEAN NOVEL VARIANTS")
    print(f"{'='*60}")
    
    top_variants = df_analysis.head(20)
    for i, (idx, variant) in enumerate(top_variants.iterrows(), 1):
        gene = variant.get('hugo', 'Unknown')
        chrom = variant.get('chrom', '?')
        pos = variant.get('pos', '?')
        so = variant.get('so', 'Unknown')
        euro_score = variant.get('european_novelty_score', 0)
        pop_score = variant.get('population_specificity_score', 0)
        classification = variant.get('european_classification', 'Unknown')
        euro_freq = variant.get('gnomad.af_nfe', 'N/A')
        
        print(f"{i:2d}. {gene} ({chrom}:{pos}) - {so}")
        print(f"    European Score: {euro_score:.1f}, Pop. Specificity: {pop_score:.1f}")
        print(f"    Classification: {classification}, European MAF: {euro_freq}")
    
    # 2. European classification analysis
    print(f"\n{'='*60}")
    print("EUROPEAN FREQUENCY CLASSIFICATION")
    print(f"{'='*60}")
    
    classification_counts = df_analysis['european_classification'].value_counts()
    print("European variant classifications:")
    for classification, count in classification_counts.items():
        percentage = (count / len(df_analysis)) * 100
        print(f"  {classification}: {count:,} ({percentage:.1f}%)")
    
    # 3. Population-specific variants
    print(f"\n{'='*60}")
    print("POPULATION-SPECIFIC VARIANTS")
    print(f"{'='*60}")
    
    pop_specific = df_analysis[df_analysis['european_classification'].str.contains('population_specific', na=False)]
    print(f"Population-specific variants: {len(pop_specific):,}")
    
    if len(pop_specific) > 0:
        print("\nTop 10 population-specific variants:")
        for i, (idx, variant) in enumerate(pop_specific.head(10).iterrows(), 1):
            gene = variant.get('hugo', 'Unknown')
            euro_freq = variant.get('gnomad.af_nfe', 'N/A')
            afr_freq = variant.get('gnomad.af_afr', 'N/A')
            eas_freq = variant.get('gnomad.af_eas', 'N/A')
            print(f"{i:2d}. {gene}: EUR={euro_freq}, AFR={afr_freq}, EAS={eas_freq}")
    
    # Save analysis results
    df_analysis.to_csv(folders['data'] / 'european_novel_variants.csv', index=False)
    print(f"\n✓ Saved European analysis to: {folders['data'] / 'european_novel_variants.csv'}")
    
    return df_analysis

def create_european_summary_tables(df_analysis, folders):
    """Create European population-specific summary tables"""
    print(f"\n{'='*60}")
    print("CREATING EUROPEAN VARIANT SUMMARY TABLES")
    print(f"{'='*60}")
    
    # 1. Top 100 European novel variants
    top_100 = df_analysis.head(100)
    top_100.to_csv(folders['data'] / 'top_100_european_novel_variants.csv', index=False)
    print(f"✓ Saved top 100 European novel variants")
    
    # 2. Ultra-high European novelty score variants
    ultra_high = df_analysis[df_analysis['european_novelty_score'] > 80]
    if len(ultra_high) > 0:
        ultra_high.to_csv(folders['data'] / 'ultra_high_european_novelty.csv', index=False)
        print(f"✓ Saved {len(ultra_high):,} ultra-high European novelty variants")
    
    # 3. Population-specific variants
    pop_specific = df_analysis[df_analysis['european_classification'].str.contains('population_specific', na=False)]
    if len(pop_specific) > 0:
        pop_specific.to_csv(folders['data'] / 'population_specific_variants.csv', index=False)
        print(f"✓ Saved {len(pop_specific):,} population-specific variants")
    
    # 4. Ultra-rare European variants
    ultra_rare = df_analysis[df_analysis['european_classification'] == 'ultra_rare_european']
    if len(ultra_rare) > 0:
        ultra_rare.to_csv(folders['data'] / 'ultra_rare_european_variants.csv', index=False)
        print(f"✓ Saved {len(ultra_rare):,} ultra-rare European variants")
    
    # 5. European frequency comparison table
    freq_comparison = []
    freq_cols = get_european_frequency_columns()
    
    for idx, row in df_analysis.head(50).iterrows():
        entry = {
            'gene': row.get('hugo', 'Unknown'),
            'chrom_pos': f"{row.get('chrom', '?')}:{row.get('pos', '?')}",
            'european_score': row.get('european_novelty_score', 0),
            'gnomad_nfe': row.get('gnomad.af_nfe', np.nan),
            'gnomad_fin': row.get('gnomad.af_fin', np.nan),
            'gbr_1000g': row.get('thousandgenomes_european.gbr_af', np.nan),
            'gnomad_afr': row.get('gnomad.af_afr', np.nan),
            'gnomad_eas': row.get('gnomad.af_eas', np.nan),
            'classification': row.get('european_classification', 'Unknown')
        }
        freq_comparison.append(entry)
    
    freq_df = pd.DataFrame(freq_comparison)
    freq_df.to_csv(folders['population_analysis'] / 'european_frequency_comparison.csv', index=False)
    print(f"✓ Saved European frequency comparison table")
    
    return top_100, ultra_high, pop_specific

def create_european_visualizations(df_analysis, folders):
    """Create European population-specific visualizations"""
    print(f"\n{'='*60}")
    print("CREATING EUROPEAN VISUALIZATIONS")
    print(f"{'='*60}")
    
    plt.style.use('seaborn-v0_8')
    
    # 1. European Novelty Score Distribution
    plt.figure(figsize=(14, 8))
    plt.hist(df_analysis['european_novelty_score'], bins=50, alpha=0.7, color='darkgreen', edgecolor='black')
    plt.axvline(x=df_analysis['european_novelty_score'].mean(), color='red', linestyle='--', 
                label=f'Mean: {df_analysis["european_novelty_score"].mean():.1f}')
    plt.axvline(x=df_analysis['european_novelty_score'].quantile(0.95), color='orange', linestyle='--',
                label=f'95th percentile: {df_analysis["european_novelty_score"].quantile(0.95):.1f}')
    plt.title('European Population-Specific Novelty Score Distribution', fontsize=16, fontweight='bold')
    plt.xlabel('European Novelty Score', fontsize=12)
    plt.ylabel('Number of Variants', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(folders['plots'] / 'european_novelty_score_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ European novelty score distribution saved")
    
    # 2. European Classification Pie Chart
    plt.figure(figsize=(12, 8))
    classification_counts = df_analysis['european_classification'].value_counts()
    colors = plt.cm.Set3(np.linspace(0, 1, len(classification_counts)))
    
    wedges, texts, autotexts = plt.pie(classification_counts.values, labels=classification_counts.index, 
                                      autopct='%1.1f%%', colors=colors, startangle=90)
    plt.title('European Frequency Classification of Novel Variants', fontsize=16, fontweight='bold')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(folders['plots'] / 'european_classification_pie.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ European classification pie chart saved")
    
    # 3. Population Frequency Comparison
    if all(col in df_analysis.columns for col in ['gnomad.af_nfe', 'gnomad.af_afr']):
        plt.figure(figsize=(12, 8))
        
        # Filter for variants with both frequencies
        comparison_data = df_analysis[
            df_analysis['gnomad.af_nfe'].notna() & 
            df_analysis['gnomad.af_afr'].notna()
        ].copy()
        
        if len(comparison_data) > 0:
            euro_freq = pd.to_numeric(comparison_data['gnomad.af_nfe'], errors='coerce')
            afr_freq = pd.to_numeric(comparison_data['gnomad.af_afr'], errors='coerce')
            
            # Remove zeros for log scale
            mask = (euro_freq > 0) & (afr_freq > 0)
            euro_freq = euro_freq[mask]
            afr_freq = afr_freq[mask]
            
            if len(euro_freq) > 0:
                scatter = plt.scatter(euro_freq, afr_freq, 
                                    c=comparison_data.loc[mask, 'european_novelty_score'], 
                                    cmap='viridis', alpha=0.6, s=30)
                plt.colorbar(scatter, label='European Novelty Score')
                plt.plot([1e-6, 1], [1e-6, 1], 'r--', alpha=0.5, label='Equal frequency')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('European Frequency (gnomAD NFE)', fontsize=12)
                plt.ylabel('African Frequency (gnomAD AFR)', fontsize=12)
                plt.title('European vs African Population Frequencies', fontsize=16, fontweight='bold')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(folders['plots'] / 'european_vs_african_frequencies.png', dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ European vs African frequency comparison saved")
    
    # 4. Population Specificity Scores
    if 'population_specificity_score' in df_analysis.columns:
        plt.figure(figsize=(12, 8))
        pop_scores = df_analysis['population_specificity_score']
        pop_scores_nonzero = pop_scores[pop_scores > 0]
        
        if len(pop_scores_nonzero) > 0:
            plt.hist(pop_scores_nonzero, bins=30, alpha=0.7, color='purple', edgecolor='black')
            plt.title('Population Specificity Score Distribution', fontsize=16, fontweight='bold')
            plt.xlabel('Population Specificity Score', fontsize=12)
            plt.ylabel('Number of Variants', fontsize=12)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(folders['plots'] / 'population_specificity_scores.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Population specificity scores saved")

def create_european_dashboard(df_analysis, folders):
    """Create comprehensive European analysis dashboard"""
    print("\nCreating European analysis dashboard...")
    
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle('European Population-Specific Novel Variants Dashboard', fontsize=20, fontweight='bold', y=0.95)
    
    # 1. European Novelty Score Distribution
    ax1 = plt.subplot(3, 3, 1)
    ax1.hist(df_analysis['european_novelty_score'], bins=30, alpha=0.7, color='darkgreen')
    ax1.set_title('European Novelty Scores', fontweight='bold')
    ax1.set_xlabel('Score')
    ax1.set_ylabel('Count')
    ax1.grid(True, alpha=0.3)
    
    # 2. European Classifications
    ax2 = plt.subplot(3, 3, 2)
    class_counts = df_analysis['european_classification'].value_counts().head(6)
    class_counts.plot(kind='bar', ax=ax2, color='steelblue')
    ax2.set_title('European Classifications', fontweight='bold')
    ax2.set_xlabel('Classification')
    ax2.set_ylabel('Count')
    ax2.tick_params(axis='x', rotation=45)
    
    # 3. Top Genes
    if 'hugo' in df_analysis.columns:
        ax3 = plt.subplot(3, 3, 3)
        gene_counts = df_analysis['hugo'].dropna().value_counts().head(8)
        gene_counts.plot(kind='barh', ax=ax3, color='lightcoral')
        ax3.set_title('Top Genes', fontweight='bold')
        ax3.set_xlabel('Novel Variants')
    
    # 4. European vs Global Frequencies
    if all(col in df_analysis.columns for col in ['gnomad.af_nfe', 'gnomad.af']):
        ax4 = plt.subplot(3, 3, 4)
        euro_data = pd.to_numeric(df_analysis['gnomad.af_nfe'], errors='coerce').dropna()
        if len(euro_data) > 0:
            euro_nonzero = euro_data[euro_data > 0]
            if len(euro_nonzero) > 0:
                ax4.hist(euro_nonzero, bins=20, alpha=0.7, color='blue')
                ax4.set_xscale('log')
                ax4.set_title('European Frequencies', fontweight='bold')
                ax4.set_xlabel('gnomAD NFE AF')
                ax4.set_ylabel('Count')
    
    # 5. Population Specificity
    if 'population_specificity_score' in df_analysis.columns:
        ax5 = plt.subplot(3, 3, 5)
        pop_scores = df_analysis['population_specificity_score']
        pop_nonzero = pop_scores[pop_scores > 0]
        if len(pop_nonzero) > 0:
            ax5.hist(pop_nonzero, bins=15, alpha=0.7, color='purple')
            ax5.set_title('Population Specificity', fontweight='bold')
            ax5.set_xlabel('Score')
            ax5.set_ylabel('Count')
    
    # 6. CADD Distribution
    if 'cadd.phred' in df_analysis.columns:
        ax6 = plt.subplot(3, 3, 6)
        cadd_data = pd.to_numeric(df_analysis['cadd.phred'], errors='coerce').dropna()
        ax6.hist(cadd_data, bins=20, alpha=0.7, color='red')
        ax6.axvline(x=25, color='black', linestyle='--')
        ax6.set_title('CADD Scores', fontweight='bold')
        ax6.set_xlabel('CADD Score')
        ax6.set_ylabel('Count')
    
    # 7. Summary Statistics
    ax7 = plt.subplot(3, 3, 7)
    summary_text = f"EUROPEAN NOVEL VARIANTS\n\n"
    summary_text += f"Total Variants: {len(df_analysis):,}\n"
    summary_text += f"Mean Euro Score: {df_analysis['european_novelty_score'].mean():.1f}\n"
    
    if 'hugo' in df_analysis.columns:
        unique_genes = df_analysis['hugo'].dropna().nunique()
        summary_text += f"Genes: {unique_genes:,}\n"
    
    pop_specific_count = len(df_analysis[df_analysis['european_classification'].str.contains('population_specific', na=False)])
    summary_text += f"Pop-specific: {pop_specific_count:,}\n"
    
    ultra_rare_count = len(df_analysis[df_analysis['european_classification'] == 'ultra_rare_european'])
    summary_text += f"Ultra-rare EUR: {ultra_rare_count:,}\n"
    
    ax7.text(0.1, 0.9, summary_text, transform=ax7.transAxes, 
             fontsize=11, verticalalignment='top', fontfamily='monospace')
    ax7.axis('off')
    
    # 8. Top European Variants Table
    ax8 = plt.subplot(3, 3, 8)
    top_5 = df_analysis.head(5)
    table_data = []
    for _, variant in top_5.iterrows():
        gene = str(variant.get('hugo', 'Unknown'))[:8]
        score = f"{variant.get('european_novelty_score', 0):.1f}"
        classification = str(variant.get('european_classification', 'Unknown'))[:12]
        table_data.append([gene, score, classification])
    
    table = ax8.table(cellText=table_data,
                     colLabels=['Gene', 'Score', 'Classification'],
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    ax8.set_title('Top 5 European Variants', fontweight='bold')
    ax8.axis('off')
    
    plt.tight_layout()
    plt.savefig(folders['plots'] / 'european_comprehensive_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ European comprehensive dashboard saved")

def generate_european_report(df_analysis, folders):
    """Generate detailed European population-specific report"""
    print("\nGenerating European population-specific report...")
    
    report_path = folders['reports'] / 'european_novel_variants_report.md'
    
    with open(report_path, 'w') as f:
        f.write("# European Population-Specific Novel Variants Report\n")
        f.write("## OpenCRAVAT WES Analysis - European Population Focus\n\n")
        
        f.write("### Executive Summary\n")
        f.write(f"- **Total European novel variants**: {len(df_analysis):,}\n")
        f.write(f"- **Mean European novelty score**: {df_analysis['european_novelty_score'].mean():.1f}\n")
        f.write(f"- **Population focus**: European genetic background\n")
        f.write(f"- **Frequency sources**: gnomAD NFE/FIN, 1000 Genomes European populations\n\n")
        
        # Population-specific variants
        pop_specific = df_analysis[df_analysis['european_classification'].str.contains('population_specific', na=False)]
        if len(pop_specific) > 0:
            f.write(f"- **Population-specific variants**: {len(pop_specific):,} (rare in Europeans, common elsewhere)\n\n")
        
        f.write("### Top 10 European Novel Variants\n\n")
        f.write("| Rank | Gene | Position | European Score | Classification | EUR Freq | Other Pop |\n")
        f.write("|------|------|----------|----------------|----------------|----------|----------|\n")
        
        top_10 = df_analysis.head(10)
        for i, (_, variant) in enumerate(top_10.iterrows(), 1):
            gene = variant.get('hugo', 'Unknown')
            chrom = variant.get('chrom', '?')
            pos = variant.get('pos', '?')
            score = variant.get('european_novelty_score', 0)
            classification = variant.get('european_classification', 'Unknown')
            eur_freq = variant.get('gnomad.af_nfe', 'N/A')
            afr_freq = variant.get('gnomad.af_afr', 'N/A')
            
            f.write(f"| {i} | **{gene}** | {chrom}:{pos} | {score:.1f} | {classification} | {eur_freq} | {afr_freq} |\n")
        
        f.write("\n### European Population Strategy\n")
        f.write("This analysis uses European population-specific approach:\n\n")
        f.write("1. **Primary European frequencies**: gnomAD Non-Finnish European, Finnish\n")
        f.write("2. **Secondary European**: 1000 Genomes European populations (GBR, IBS, TSI, CEU)\n")
        f.write("3. **Cross-population validation**: Comparison with African, East Asian populations\n")
        f.write("4. **Population-specific discovery**: Variants rare in Europeans but common elsewhere\n")
        f.write("5. **European novelty scoring**: Weighted for European genetic background\n\n")
        
        # Classification breakdown
        classification_counts = df_analysis['european_classification'].value_counts()
        f.write("### European Frequency Classifications\n\n")
        for classification, count in classification_counts.items():
            percentage = (count / len(df_analysis)) * 100
            f.write(f"- **{classification}**: {count:,} ({percentage:.1f}%)\n")
        
        if len(pop_specific) > 0:
            f.write("\n### Population-Specific Variants (Top 5)\n\n")
            f.write("Variants rare in Europeans but common in other populations:\n\n")
            for i, (_, variant) in enumerate(pop_specific.head(5).iterrows(), 1):
                gene = variant.get('hugo', 'Unknown')
                eur_freq = variant.get('gnomad.af_nfe', 'N/A')
                afr_freq = variant.get('gnomad.af_afr', 'N/A')
                eas_freq = variant.get('gnomad.af_eas', 'N/A')
                f.write(f"{i}. **{gene}**: EUR={eur_freq}, AFR={afr_freq}, EAS={eas_freq}\n")
        
        f.write("\n### Files Generated\n")
        f.write("- `european_novel_variants.csv` - Complete European analysis\n")
        f.write("- `top_100_european_novel_variants.csv` - Top European candidates\n")
        f.write("- `population_specific_variants.csv` - Population-specific variants\n")
        f.write("- `ultra_rare_european_variants.csv` - Ultra-rare in Europeans\n")
        f.write("- `european_frequency_comparison.csv` - Cross-population frequencies\n")
        f.write("- European population-specific visualizations\n")
    
    print(f"✓ European report saved to: {report_path}")

def main():
    """Main European analysis function"""
    filename = "your_opencravat_data.csv"  # Replace with your OpenCRAVAT annotated file
    
    print("European Population-Specific Novel Variants Discovery")
    print("=" * 80)
    print("Focus: European population genetic background")
    print("Strategy: European frequencies + cross-population validation")
    print("Goal: Population-specific novel pathogenic variant discovery")
    
    start_time = time.time()
    
    # Create output structure
    folders = create_output_structure()
    print(f"Results will be saved to: {folders['data'].parent}")
    
    # Load and filter European novel variants
    df_novel = load_and_filter_european_novel_variants(filename)
    
    if df_novel.empty:
        print("No European novel variants found")
        return
    
    # Analyze European variants
    df_analysis = analyze_european_novel_variants(df_novel, folders)
    
    # Create summary tables
    top_100, ultra_high, pop_specific = create_european_summary_tables(df_analysis, folders)
    
    # Create visualizations
    create_european_visualizations(df_analysis, folders)
    
    # Create dashboard
    create_european_dashboard(df_analysis, folders)
    
    # Generate report
    generate_european_report(df_analysis, folders)
    
    # Final summary
    total_time = time.time() - start_time
    print(f"\n{'='*80}")
    print("EUROPEAN NOVEL VARIANTS ANALYSIS COMPLETED!")
    print(f"{'='*80}")
    print(f"Analysis time: {total_time/60:.1f} minutes")
    print(f"Total European novel variants: {len(df_analysis):,}")
    print(f"Top European novelty score: {df_analysis['european_novelty_score'].max():.1f}")
    
    # Population-specific summary
    pop_specific_count = len(df_analysis[df_analysis['european_classification'].str.contains('population_specific', na=False)])
    ultra_rare_count = len(df_analysis[df_analysis['european_classification'] == 'ultra_rare_european'])
    
    print(f"Population-specific variants: {pop_specific_count:,}")
    print(f"Ultra-rare European variants: {ultra_rare_count:,}")
    print(f"Results folder: {folders['data'].parent}")

if __name__ == "__main__":
    main() 