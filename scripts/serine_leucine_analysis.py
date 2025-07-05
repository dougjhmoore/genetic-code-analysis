#!/usr/bin/env python3
"""
Serine vs Leucine: Fine structure analysis of 6-fold codon families
Why does Serine follow FC while Leucine violates it?
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('../CUTG/AGG/refseq_codon_species.tsv', sep='\t')

# Define the two families with their sub-groups
FAMILIES = {
    'Serine': {
        'UC_group': ['TCT', 'TCC', 'TCA', 'TCG'],  # UC* subfamily
        'AG_group': ['AGT', 'AGC'],               # AG* subfamily
        'all': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
    },
    'Leucine': {
        'TT_group': ['TTA', 'TTG'],                           # UU* subfamily  
        'CT_group': ['CTT', 'CTC', 'CTA', 'CTG'],            # CU* subfamily
        'all': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
    }
}

def analyze_family_structure(family_name, groups, sample_size=500):
    """Analyze the internal structure of a 6-fold family"""
    print(f"\n{'='*60}")
    print(f"{family_name.upper()} FAMILY STRUCTURE ANALYSIS")
    print(f"{'='*60}")
    
    sample_df = df.sample(sample_size)
    
    # Calculate relative usage within each sub-group and across groups
    within_subgroup_cvs = []
    between_subgroup_cvs = []
    
    for _, row in sample_df.iterrows():
        # Get counts for all codons in family
        all_counts = [row[codon] for codon in groups['all']]
        total_family = sum(all_counts)
        
        if total_family == 0:
            continue
            
        # Calculate frequencies within each sub-group
        subgroup_means = []
        subgroup_cvs = []
        
        for subgroup_name, codons in groups.items():
            if subgroup_name == 'all':
                continue
                
            counts = [row[codon] for codon in codons]
            subgroup_total = sum(counts)
            
            if subgroup_total == 0:
                continue
                
            # Frequencies within this sub-group
            freqs = np.array(counts) / subgroup_total
            mean_freq = np.mean(freqs)
            cv = np.std(freqs, ddof=1) / mean_freq if mean_freq > 0 else 0
            
            subgroup_cvs.append(cv)
            subgroup_means.append(subgroup_total / total_family)  # Fraction of family usage
        
        # Average CV within sub-groups
        if subgroup_cvs:
            within_subgroup_cvs.append(np.mean(subgroup_cvs))
            
        # CV between sub-group usage levels
        if len(subgroup_means) > 1:
            between_cv = np.std(subgroup_means) / np.mean(subgroup_means) if np.mean(subgroup_means) > 0 else 0
            between_subgroup_cvs.append(between_cv)
    
    # Print results
    print(f"Sample size: {len(within_subgroup_cvs)} species")
    print(f"\nWithin sub-groups CV: {np.median(within_subgroup_cvs):.3f} ± {np.std(within_subgroup_cvs):.3f}")
    print(f"Between sub-groups CV: {np.median(between_subgroup_cvs):.3f} ± {np.std(between_subgroup_cvs):.3f}")
    
    # Show typical usage patterns
    print(f"\nTypical codon usage patterns:")
    
    # Get median usage for each codon
    medians = {}
    for codon in groups['all']:
        values = []
        for _, row in sample_df.iterrows():
            total = sum(row[c] for c in groups['all'])
            if total > 0:
                values.append(row[codon] / total)
        if values:
            medians[codon] = np.median(values)
    
    # Group by subfamily
    for subgroup_name, codons in groups.items():
        if subgroup_name == 'all':
            continue
        print(f"\n  {subgroup_name}:")
        for codon in codons:
            if codon in medians:
                print(f"    {codon}: {medians[codon]:.3f}")
    
    return within_subgroup_cvs, between_subgroup_cvs

# Analyze both families
ser_within, ser_between = analyze_family_structure('Serine', FAMILIES['Serine'])
leu_within, leu_between = analyze_family_structure('Leucine', FAMILIES['Leucine'])

# Compare the results
print(f"\n{'='*60}")
print("COMPARATIVE SUMMARY")
print(f"{'='*60}")

print(f"\nSerine (CV = 0.650, GOOD FC compliance):")
print(f"  Within sub-groups:  {np.median(ser_within):.3f}")
print(f"  Between sub-groups: {np.median(ser_between):.3f}")

print(f"\nLeucine (CV = 0.991, BAD FC compliance):")
print(f"  Within sub-groups:  {np.median(leu_within):.3f}")
print(f"  Between sub-groups: {np.median(leu_between):.3f}")

print(f"\nHYPOTHESIS:")
print(f"If Serine's sub-groups are balanced → low between-group CV")
print(f"If Leucine's sub-groups are imbalanced → high between-group CV")

ratio_ser = np.median(ser_between) / np.median(ser_within) if np.median(ser_within) > 0 else float('inf')
ratio_leu = np.median(leu_between) / np.median(leu_within) if np.median(leu_within) > 0 else float('inf')

print(f"\nBetween/Within ratios:")
print(f"  Serine:  {ratio_ser:.2f}")
print(f"  Leucine: {ratio_leu:.2f}")

if ratio_leu > ratio_ser:
    print(f"\n✓ CONFIRMED: Leucine shows more between-group imbalance!")
else:
    print(f"\n✗ UNEXPECTED: Need different explanation...")