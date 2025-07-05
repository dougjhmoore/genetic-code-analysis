#!/usr/bin/env python3
"""
Orbit-by-orbit First-Classness analysis
Test whether 6-fold orbits (especially Serine) violate FC more than others
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('../CUTG/AGG/refseq_codon_species.tsv', sep='\t')
codon_cols = [c for c in df.columns if len(c) == 3]

# Define orbit families with their expected FC compliance
ORBIT_FAMILIES = {
    # Size 1 (perfect FC - no variance possible)
    'Met1': ['AUG'],
    'Trp1': ['UGG'],
    
    # Size 2 (strong FC expected)
    'Phe2': ['UUU', 'UUC'],
    'Tyr2': ['UAU', 'UAC'], 
    'Cys2': ['UGU', 'UGC'],
    'His2': ['CAU', 'CAC'],
    'Gln2': ['CAA', 'CAG'],
    'Asn2': ['AAU', 'AAC'],
    'Lys2': ['AAA', 'AAG'],
    'Asp2': ['GAU', 'GAC'],
    'Glu2': ['GAA', 'GAG'],
    
    # Size 3 (moderate FC expected)
    'Stop3': ['UAA', 'UAG', 'UGA'],
    'Ile3': ['AUU', 'AUC', 'AUA'],
    
    # Size 4 (good FC expected)
    'Val4': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Pro4': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Thr4': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Ala4': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Gly4': ['GGU', 'GGC', 'GGA', 'GGG'],
    
    # Size 6 (WEAK FC expected - split families!)
    'Ser6': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],  # UC* vs AG*
    'Leu6': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],  # UU* vs CU*
    'Arg6': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # CG* vs AG*
}

# Convert to DNA codons (matching your data)
DNA_FAMILIES = {}
for name, codons in ORBIT_FAMILIES.items():
    DNA_FAMILIES[name] = [c.replace('U', 'T') for c in codons]

def calculate_within_orbit_variance(data_row, family_codons):
    """Calculate coefficient of variation within a codon family"""
    if len(family_codons) <= 1:
        return np.nan
    
    counts = [data_row[codon] for codon in family_codons if codon in data_row.index]
    if len(counts) <= 1 or sum(counts) == 0:
        return np.nan
    
    # Convert to frequencies
    freqs = np.array(counts) / sum(counts)
    
    # Calculate coefficient of variation (std/mean)
    mean_freq = np.mean(freqs)
    if mean_freq == 0:
        return np.nan
    
    cv = np.std(freqs, ddof=1) / mean_freq
    return cv

# Analyze each orbit family
results = {}
sample_size = min(1000, len(df))  # Use subset for speed
sample_df = df.sample(sample_size)

print(f"Analyzing {sample_size} species across {len(DNA_FAMILIES)} orbit families...")

for family_name, family_codons in DNA_FAMILIES.items():
    print(f"Processing {family_name}...")
    
    cv_values = []
    for _, row in sample_df.iterrows():
        cv = calculate_within_orbit_variance(row, family_codons)
        if np.isfinite(cv):
            cv_values.append(cv)
    
    if cv_values:
        results[family_name] = {
            'median_cv': np.median(cv_values),
            'mean_cv': np.mean(cv_values),
            'std_cv': np.std(cv_values),
            'n_species': len(cv_values),
            'values': cv_values
        }

# Sort by orbit size and median CV
orbit_sizes = {'1': [], '2': [], '3': [], '4': [], '6': []}
for name, stats in results.items():
    size = name[-1]  # Last character is the size
    if size in orbit_sizes:
        orbit_sizes[size].append((name, stats['median_cv']))

print("\n" + "="*60)
print("ORBIT-BY-ORBIT FIRST-CLASSNESS COMPLIANCE")
print("="*60)
print("Lower CV = Better FC compliance (more uniform usage)")
print()

for size in ['1', '2', '3', '4', '6']:
    if orbit_sizes[size]:
        print(f"\n{size}-FOLD ORBITS:")
        sorted_orbits = sorted(orbit_sizes[size], key=lambda x: x[1])
        for name, median_cv in sorted_orbits:
            stats = results[name]
            print(f"  {name:8s}: CV = {median_cv:.3f} ± {stats['std_cv']:.3f} (n={stats['n_species']})")

# Special focus on 6-fold families
print("\n" + "="*60) 
print("6-FOLD FAMILY DETAILED ANALYSIS")
print("="*60)
print("Prediction: Ser6 should have HIGHEST CV (worst FC compliance)")
print()

sixfold_results = [(name, stats) for name, stats in results.items() if name.endswith('6')]
sixfold_sorted = sorted(sixfold_results, key=lambda x: x[1]['median_cv'], reverse=True)

for i, (name, stats) in enumerate(sixfold_sorted):
    rank = "WORST" if i == 0 else "BEST" if i == len(sixfold_sorted)-1 else f"#{i+1}"
    print(f"  {rank:5s} {name}: CV = {stats['median_cv']:.3f}")

# Test prediction: Is Serine the worst?
ser_cv = results.get('Ser6', {}).get('median_cv', None)
if ser_cv:
    other_6fold = [stats['median_cv'] for name, stats in results.items() 
                   if name.endswith('6') and name != 'Ser6']
    if other_6fold:
        is_worst = ser_cv > max(other_6fold)
        print(f"\nPREDICTION TEST:")
        print(f"Serine CV = {ser_cv:.3f}")
        print(f"Other 6-fold max = {max(other_6fold):.3f}")
        print(f"Serine is worst FC violator: {is_worst}")

print("\n" + "="*60)