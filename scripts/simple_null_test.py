#!/usr/bin/env python3
"""
Simple null model test for First-Classness
"""
import numpy as np
import pandas as pd

# Your real result
OBSERVED = 4.810

# Load your actual codon data
df = pd.read_csv('../CUTG/AGG/refseq_codon_species.tsv', sep='\t')
codon_cols = [c for c in df.columns if len(c) == 3 and c.replace('T','U') in 
              [a+b+c for a in "UCAG" for b in "UCAG" for c in "UCAG"]]

# Convert to frequencies and take sample of rows
sample_size = min(1000, len(df))  # Use 1000 species for speed
sample_data = df[codon_cols].sample(sample_size).values

# Real orbit assignments (matching your orbit_map.csv)
real_orbits = [2,2,6,6,6,6,6,6,2,2,3,3,2,2,3,1,  # U row
               6,6,6,6,4,4,4,4,2,2,2,2,6,6,6,6,  # C row  
               3,3,3,1,4,4,4,4,2,2,2,2,6,6,6,6,  # A row
               4,4,4,4,4,4,4,4,2,2,2,2,4,4,4,4]  # G row

print(f"Testing {sample_size} species against {OBSERVED}")

# Run null model
null_ratios = []
rng = np.random.default_rng(2025)

for trial in range(1000):
    if trial % 100 == 0:
        print(f"  Trial {trial}...")
    
    # Shuffle orbit assignments
    fake_orbits = real_orbits.copy()
    rng.shuffle(fake_orbits)
    
    trial_ratios = []
    for row in sample_data:
        if row.sum() == 0:
            continue
            
        # Calculate RSCU with fake orbits
        freq = row / row.sum()
        rscu = freq * fake_orbits
        
        # Group by fake orbit
        groups = {}
        for i, orbit in enumerate(fake_orbits):
            groups.setdefault(orbit, []).append(rscu[i])
        
        # Calculate σ_intra / σ_inter
        intra_vals = []
        inter_vals = []
        for orbit, values in groups.items():
            if len(values) > 1:
                mean_val = np.mean(values)
                intra_vals.extend([v - mean_val for v in values])
                inter_vals.extend([mean_val] * len(values))
        
        if len(intra_vals) > 1 and len(inter_vals) > 1:
            ratio = np.std(intra_vals, ddof=1) / np.std(inter_vals, ddof=1)
            if np.isfinite(ratio):
                trial_ratios.append(ratio)
    
    if trial_ratios:
        null_ratios.append(np.median(trial_ratios))

# Calculate statistics
null_ratios = np.array(null_ratios)
null_mean = np.mean(null_ratios)
null_std = np.std(null_ratios)
z_score = (OBSERVED - null_mean) / null_std
p_value = (null_ratios >= OBSERVED).mean()

print(f"\nResults:")
print(f"Null mean: {null_mean:.3f}")
print(f"Null std:  {null_std:.3f}")
print(f"Z-score:   {z_score:.2f}")
print(f"P-value:   {p_value:.4f}")
print(f"Observed:  {OBSERVED}")