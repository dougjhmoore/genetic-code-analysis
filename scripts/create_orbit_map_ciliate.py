#!/usr/bin/env python3
"""
Create orbit_map_ciliate.csv for ciliate genetic code analysis
Ciliate code changes: UAA, UAG (stops) → Glutamine
"""

# Complete ciliate genetic code orbit assignments
orbit_assignments = [
    # UNN codons
    ('UUU', 2), ('UUC', 2),  # Phe
    ('UUA', 6), ('UUG', 6),  # Leu
    ('UCU', 6), ('UCC', 6), ('UCA', 6), ('UCG', 6),  # Ser
    ('UAU', 2), ('UAC', 2),  # Tyr
    ('UAA', 4), ('UAG', 4),  # CHANGED: Stop → Gln (now 4-fold)
    ('UGU', 2), ('UGC', 2),  # Cys
    ('UGA', 1),              # CHANGED: Sole remaining stop (1-fold)
    ('UGG', 1),              # Trp
    
    # CNN codons
    ('CUU', 6), ('CUC', 6), ('CUA', 6), ('CUG', 6),  # Leu
    ('CCU', 4), ('CCC', 4), ('CCA', 4), ('CCG', 4),  # Pro
    ('CAU', 2), ('CAC', 2),  # His
    ('CAA', 4), ('CAG', 4),  # CHANGED: Gln (now 4-fold with UAA,UAG)
    ('CGU', 6), ('CGC', 6), ('CGA', 6), ('CGG', 6),  # Arg
    
    # ANN codons
    ('AUU', 3), ('AUC', 3), ('AUA', 3),  # Ile
    ('AUG', 1),              # Met
    ('ACU', 4), ('ACC', 4), ('ACA', 4), ('ACG', 4),  # Thr
    ('AAU', 2), ('AAC', 2),  # Asn
    ('AAA', 2), ('AAG', 2),  # Lys
    ('AGU', 6), ('AGC', 6),  # Ser
    ('AGA', 6), ('AGG', 6),  # Arg
    
    # GNN codons
    ('GUU', 4), ('GUC', 4), ('GUA', 4), ('GUG', 4),  # Val
    ('GCU', 4), ('GCC', 4), ('GCA', 4), ('GCG', 4),  # Ala
    ('GAU', 2), ('GAC', 2),  # Asp
    ('GAA', 2), ('GAG', 2),  # Glu
    ('GGU', 4), ('GGC', 4), ('GGA', 4), ('GGG', 4),  # Gly
]

# Write to CSV file
with open('orbit_map_ciliate.csv', 'w') as f:
    for codon, orbit in orbit_assignments:
        f.write(f"{codon},{orbit}\n")

print("Created orbit_map_ciliate.csv")
print("Key changes from standard genetic code:")
print("- UAA, UAG: Stop → Gln (creates 4-fold Gln family)")
print("- UGA: Sole remaining stop (1-fold)")
print("- Total orbit structure: 1-3-4-6 preserved with modifications")