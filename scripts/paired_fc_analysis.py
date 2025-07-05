#!/usr/bin/env python3
"""
paired_fc_analysis.py - Test FC predictions using paired nuclear/mitochondrial data
This is the smoking gun test for First-Classness theory
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
from pathlib import Path

def load_orbit_map(orbit_file):
    """Load orbit mapping from CSV"""
    orbit_df = pd.read_csv(orbit_file)
    return dict(zip(orbit_df['codon'], orbit_df['orbit_id']))

def convert_dna_to_rna(codon):
    """Convert DNA codon to RNA"""
    return codon.replace('T', 'U')

def calculate_fc_compliance(row, orbit_map, codon_cols):
    """Calculate FC compliance for a single organism"""
    
    # Group codons by orbit
    orbit_codons = {}
    for codon in codon_cols:
        rna_codon = convert_dna_to_rna(codon)
        if rna_codon in orbit_map:
            orbit = orbit_map[rna_codon]
            if orbit not in orbit_codons:
                orbit_codons[orbit] = []
            orbit_codons[orbit].append(codon)
    
    # Calculate RSCU for each orbit
    intra_dispersions = []
    orbit_means = []
    
    for orbit, codons in orbit_codons.items():
        if len(codons) > 1:  # Only multi-codon families
            family_total = sum(row[c] for c in codons)
            if family_total > 0:
                # Calculate RSCU values
                rscus = []
                for codon in codons:
                    expected = family_total / len(codons)
                    rscu = row[codon] / expected if expected > 0 else 0
                    rscus.append(rscu)
                
                # Within-orbit dispersion
                intra_dispersions.append(np.std(rscus))
                orbit_means.append(np.mean(rscus))
    
    # Calculate FC metrics
    sigma_intra = np.mean(intra_dispersions) if intra_dispersions else np.nan
    sigma_inter = np.std(orbit_means) if len(orbit_means) > 1 else np.nan
    fc_ratio = sigma_intra / sigma_inter if sigma_inter > 0 else np.nan
    
    return {
        'sigma_intra': sigma_intra,
        'sigma_inter': sigma_inter, 
        'fc_ratio': fc_ratio,
        'n_orbits': len(orbit_means)
    }

def paired_fc_analysis(nuclear_file, mito_file, nuclear_orbit_map, mito_orbit_map):
    """Perform paired FC analysis on nuclear vs mitochondrial data"""
    
    print("=== PAIRED FC ANALYSIS ===")
    
    # Load datasets
    nuclear_df = pd.read_csv(nuclear_file, sep='\t')
    mito_df = pd.read_csv(mito_file, sep='\t')
    
    print(f"Nuclear ciliates: {len(nuclear_df)} organisms")
    print(f"Mitochondrial ciliates: {len(mito_df)} organisms")
    
    # Get codon columns
    codon_cols = [col for col in nuclear_df.columns if len(col) == 3 and col.isalpha()]
    
    # Calculate FC compliance for all organisms
    print("\nCalculating FC compliance...")
    
    nuclear_results = []
    for idx, row in nuclear_df.iterrows():
        fc_metrics = calculate_fc_compliance(row, nuclear_orbit_map, codon_cols)
        nuclear_results.append({
            'Taxid': row['Taxid'],
            'Species': row.get('Species_Name', 'Unknown'),
            'Type': 'Nuclear',
            **fc_metrics
        })
    
    mito_results = []
    for idx, row in mito_df.iterrows():
        fc_metrics = calculate_fc_compliance(row, mito_orbit_map, codon_cols)
        mito_results.append({
            'Taxid': row['Taxid'],
            'Species': row.get('Species_Name', 'Unknown'),
            'Type': 'Mitochondrial',
            **fc_metrics
        })
    
    # Combine results
    all_results = pd.DataFrame(nuclear_results + mito_results)
    
    # Find paired organisms (same species, both nuclear and mitochondrial)
    nuclear_species = set(all_results[all_results['Type'] == 'Nuclear']['Species'])
    mito_species = set(all_results[all_results['Type'] == 'Mitochondrial']['Species'])
    paired_species = nuclear_species & mito_species
    
    print(f"\nPaired species analysis:")
    print(f"  Species with nuclear data: {len(nuclear_species)}")
    print(f"  Species with mitochondrial data: {len(mito_species)}")
    print(f"  Species with BOTH: {len(paired_species)}")
    
    # Statistical analysis
    nuclear_ratios = all_results[all_results['Type'] == 'Nuclear']['fc_ratio'].dropna()
    mito_ratios = all_results[all_results['Type'] == 'Mitochondrial']['fc_ratio'].dropna()
    
    print(f"\n=== OVERALL FC COMPLIANCE ===")
    print(f"Nuclear ciliates:")
    print(f"  Mean FC ratio: {nuclear_ratios.mean():.3f} ± {nuclear_ratios.std():.3f}")
    print(f"  Median FC ratio: {nuclear_ratios.median():.3f}")
    print(f"  Range: {nuclear_ratios.min():.3f} - {nuclear_ratios.max():.3f}")
    
    print(f"\nMitochondrial ciliates:")
    print(f"  Mean FC ratio: {mito_ratios.mean():.3f} ± {mito_ratios.std():.3f}")
    print(f"  Median FC ratio: {mito_ratios.median():.3f}")
    print(f"  Range: {mito_ratios.min():.3f} - {mito_ratios.max():.3f}")
    
    # Statistical test
    if len(nuclear_ratios) > 0 and len(mito_ratios) > 0:
        statistic, pvalue = stats.mannwhitneyu(nuclear_ratios, mito_ratios, alternative='less')
        print(f"\n=== FC PREDICTION TEST ===")
        print(f"Mann-Whitney U test (nuclear < mitochondrial):")
        print(f"  Statistic: {statistic}")
        print(f"  P-value: {pvalue:.2e}")
        
        if pvalue < 0.05:
            print("✓ FC PREDICTION CONFIRMED: Nuclear ciliates show significantly better FC compliance")
        else:
            print("✗ FC prediction not supported")
    
    # Paired organism analysis
    if paired_species:
        print(f"\n=== PAIRED ORGANISM ANALYSIS ===")
        paired_comparisons = []
        
        for species in paired_species:
            nuclear_data = all_results[(all_results['Species'] == species) & (all_results['Type'] == 'Nuclear')]
            mito_data = all_results[(all_results['Species'] == species) & (all_results['Type'] == 'Mitochondrial')]
            
            if len(nuclear_data) > 0 and len(mito_data) > 0:
                nuclear_fc = nuclear_data['fc_ratio'].iloc[0]
                mito_fc = mito_data['fc_ratio'].iloc[0]
                
                if not (np.isnan(nuclear_fc) or np.isnan(mito_fc)):
                    paired_comparisons.append({
                        'Species': species,
                        'Nuclear_FC': nuclear_fc,
                        'Mitochondrial_FC': mito_fc,
                        'Difference': nuclear_fc - mito_fc,
                        'FC_Advantage': mito_fc / nuclear_fc if nuclear_fc > 0 else np.nan
                    })
        
        if paired_comparisons:
            paired_df = pd.DataFrame(paired_comparisons)
            print(f"Paired comparisons: {len(paired_df)} species")
            
            # Test if nuclear is consistently better
            improvements = paired_df['Difference'] < 0  # Nuclear better (lower ratio)
            n_improved = improvements.sum()
            
            print(f"Species where nuclear < mitochondrial: {n_improved}/{len(paired_df)} ({100*n_improved/len(paired_df):.1f}%)")
            
            # Wilcoxon signed-rank test for paired data
            if len(paired_df) > 5:
                statistic, pvalue = stats.wilcoxon(paired_df['Nuclear_FC'], paired_df['Mitochondrial_FC'], alternative='less')
                print(f"Wilcoxon signed-rank test: p = {pvalue:.2e}")
                
                if pvalue < 0.05:
                    print("✓ SMOKING GUN: Nuclear consistently better than mitochondrial within species")
                else:
                    print("○ Trend suggests nuclear advantage but not statistically significant")
            
            # Show top examples
            paired_df_sorted = paired_df.sort_values('FC_Advantage', ascending=False)
            print(f"\nTop examples of FC advantage:")
            for _, row in paired_df_sorted.head(5).iterrows():
                print(f"  {row['Species']}: {row['FC_Advantage']:.2f}x better (nuclear={row['Nuclear_FC']:.3f}, mito={row['Mitochondrial_FC']:.3f})")
    
    # Visualization
    create_fc_comparison_plots(all_results, paired_species)
    
    return all_results, paired_species

def create_fc_comparison_plots(results_df, paired_species):
    """Create comprehensive FC comparison plots"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: FC ratio distribution
    nuclear_data = results_df[results_df['Type'] == 'Nuclear']['fc_ratio'].dropna()
    mito_data = results_df[results_df['Type'] == 'Mitochondrial']['fc_ratio'].dropna()
    
    ax1.hist(nuclear_data, bins=30, alpha=0.7, label='Nuclear (Ciliate code)', color='green')
    ax1.hist(mito_data, bins=30, alpha=0.7, label='Mitochondrial (Standard code)', color='red')
    ax1.axvline(nuclear_data.median(), color='green', linestyle='--', alpha=0.8)
    ax1.axvline(mito_data.median(), color='red', linestyle='--', alpha=0.8)
    ax1.set_xlabel('FC Compliance Ratio (σ_intra/σ_inter)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('FC Compliance Distribution')
    ax1.legend()
    ax1.set_yscale('log')
    
    # Plot 2: Box plot comparison
    plot_data = [nuclear_data, mito_data]
    ax2.boxplot(plot_data, labels=['Nuclear\n(Ciliate)', 'Mitochondrial\n(Standard)'])
    ax2.set_ylabel('FC Compliance Ratio')
    ax2.set_title('FC Compliance Comparison')
    ax2.set_yscale('log')
    
    # Plot 3: Paired comparisons (if available)
    if paired_species:
        paired_data = []
        for species in paired_species:
            nuclear_fc = results_df[(results_df['Species'] == species) & (results_df['Type'] == 'Nuclear')]['fc_ratio']
            mito_fc = results_df[(results_df['Species'] == species) & (results_df['Type'] == 'Mitochondrial')]['fc_ratio']
            
            if len(nuclear_fc) > 0 and len(mito_fc) > 0:
                nuclear_val = nuclear_fc.iloc[0]
                mito_val = mito_fc.iloc[0]
                if not (np.isnan(nuclear_val) or np.isnan(mito_val)):
                    paired_data.append((nuclear_val, mito_val))
        
        if paired_data:
            nuclear_paired, mito_paired = zip(*paired_data)
            ax3.scatter(nuclear_paired, mito_paired, alpha=0.7, s=50)
            
            # Add diagonal line (y=x)
            max_val = max(max(nuclear_paired), max(mito_paired))
            ax3.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Equal FC compliance')
            
            ax3.set_xlabel('Nuclear FC Ratio')
            ax3.set_ylabel('Mitochondrial FC Ratio')
            ax3.set_title(f'Paired Comparison ({len(paired_data)} species)')
            ax3.legend()
            ax3.set_xscale('log')
            ax3.set_yscale('log')
    
    # Plot 4: FC ratio by organism type
    sns.violinplot(data=results_df, x='Type', y='fc_ratio', ax=ax4)
    ax4.set_ylabel('FC Compliance Ratio')
    ax4.set_title('FC Compliance by Genetic Code Type')
    ax4.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('ciliate_fc_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    if len(sys.argv) != 5:
        print("Usage: python paired_fc_analysis.py <nuclear_file> <mito_file> <nuclear_orbit_map> <mito_orbit_map>")
        sys.exit(1)
    
    nuclear_file, mito_file, nuclear_orbit_file, mito_orbit_file = sys.argv[1:5]
    
    # Load orbit maps
    nuclear_orbit_map = load_orbit_map(nuclear_orbit_file)
    mito_orbit_map = load_orbit_map(mito_orbit_file)
    
    # Perform analysis
    results, paired_species = paired_fc_analysis(nuclear_file, mito_file, nuclear_orbit_map, mito_orbit_map)
    
    print(f"\n=== SUMMARY ===")
    print(f"This analysis tests the core FC prediction:")
    print(f"  Non-standard genetic codes that eliminate problematic orbits")
    print(f"  should show better FC compliance than standard codes")
    print(f"  even when comparing the same organisms.")

if __name__ == "__main__":
    main()