#!/usr/bin/env python3
"""
extract_ciliates.py - Extract ciliate data for FC testing
Creates nuclear and mitochondrial datasets from GenBank codon usage data
"""

import pandas as pd
import sys
from pathlib import Path

def extract_ciliate_data(codon_file, index_file, output_dir):
    """Extract ciliate nuclear and mitochondrial data"""
    
    print("=== CILIATE DATA EXTRACTION ===")
    
    # Load index to identify ciliate organisms
    print(f"Loading index: {index_file}")
    index_df = pd.read_csv(index_file, sep='\t', names=['Row', 'Database', 'NA', 'Taxid', 'Species', 'Organelle'])
    
    # Identify ciliates
    ciliate_genera = ['Paramecium', 'Tetrahymena', 'Oxytricha', 'Stentor']
    ciliate_mask = index_df['Species'].str.contains('|'.join(ciliate_genera), case=False, na=False)
    ciliate_index = index_df[ciliate_mask]
    
    print(f"Found {len(ciliate_index)} ciliate entries:")
    print(ciliate_index['Organelle'].value_counts())
    
    # Separate nuclear and mitochondrial
    nuclear_ciliates = ciliate_index[ciliate_index['Organelle'] == 'genomic']
    mito_ciliates = ciliate_index[ciliate_index['Organelle'] == 'mitochondrion']
    
    print(f"\nNuclear ciliates: {len(nuclear_ciliates)}")
    print(f"Mitochondrial ciliates: {len(mito_ciliates)}")
    
    # Load codon usage data
    print(f"\nLoading codon usage data: {codon_file}")
    codon_df = pd.read_csv(codon_file, sep='\t')
    
    # Extract ciliate data by Taxid
    nuclear_taxids = set(nuclear_ciliates['Taxid'].astype(str))
    mito_taxids = set(mito_ciliates['Taxid'].astype(str))
    
    # Match by Taxid and Organelle
    nuclear_data = codon_df[
        (codon_df['Taxid'].astype(str).isin(nuclear_taxids)) & 
        (codon_df['Organelle'] == 'genomic')
    ].copy()
    
    mito_data = codon_df[
        (codon_df['Taxid'].astype(str).isin(mito_taxids)) & 
        (codon_df['Organelle'] == 'mitochondrion')
    ].copy()
    
    print(f"\nExtracted nuclear ciliate data: {len(nuclear_data)} organisms")
    print(f"Extracted mitochondrial ciliate data: {len(mito_data)} organisms")
    
    # Add species names for better identification
    taxid_to_species = dict(zip(ciliate_index['Taxid'].astype(str), ciliate_index['Species']))
    nuclear_data['Species_Name'] = nuclear_data['Taxid'].astype(str).map(taxid_to_species)
    mito_data['Species_Name'] = mito_data['Taxid'].astype(str).map(taxid_to_species)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save datasets
    nuclear_file = output_path / 'ciliate_nuclear_codon_usage.tsv'
    mito_file = output_path / 'ciliate_mitochondrial_codon_usage.tsv'
    
    nuclear_data.to_csv(nuclear_file, sep='\t', index=False)
    mito_data.to_csv(mito_file, sep='\t', index=False)
    
    print(f"\n✓ Nuclear ciliate data saved: {nuclear_file}")
    print(f"✓ Mitochondrial ciliate data saved: {mito_file}")
    
    # Create species overlap analysis
    nuclear_species = set(nuclear_data['Species_Name'].dropna())
    mito_species = set(mito_data['Species_Name'].dropna())
    overlap_species = nuclear_species & mito_species
    
    print(f"\n=== SPECIES OVERLAP ANALYSIS ===")
    print(f"Nuclear-only species: {len(nuclear_species - mito_species)}")
    print(f"Mitochondrial-only species: {len(mito_species - nuclear_species)}")
    print(f"Both nuclear AND mitochondrial: {len(overlap_species)}")
    
    if overlap_species:
        print(f"\nSpecies with both datasets (perfect for FC testing):")
        for species in sorted(overlap_species):
            print(f"  - {species}")
    
    # Save overlap list for paired analysis
    if overlap_species:
        overlap_file = output_path / 'ciliate_paired_species.txt'
        with open(overlap_file, 'w') as f:
            for species in sorted(overlap_species):
                f.write(f"{species}\n")
        print(f"\n✓ Paired species list saved: {overlap_file}")
    
    return nuclear_file, mito_file, overlap_species

def create_ciliate_orbit_map(output_dir):
    """Create orbit mapping for ciliate genetic code"""
    
    print("\n=== CREATING CILIATE ORBIT MAP ===")
    
    # Ciliate genetic code changes:
    # UAA: Stop → Gln 
    # UAG: Stop → Gln
    # UGA: remains Stop (sole terminator)
    
    orbit_map = {
        # 1-fold orbits
        'AUG': ('Met', 1),
        'UGG': ('Trp', 1), 
        'UGA': ('Stop', 1),  # Sole stop in ciliates
        
        # 2-fold orbits (standard)
        'UUU': ('Phe', 2), 'UUC': ('Phe', 2),
        'UAU': ('Tyr', 2), 'UAC': ('Tyr', 2),
        'UGU': ('Cys', 2), 'UGC': ('Cys', 2),
        'CAU': ('His', 2), 'CAC': ('His', 2),
        'AAU': ('Asn', 2), 'AAC': ('Asn', 2),
        'AAA': ('Lys', 2), 'AAG': ('Lys', 2),
        'GAU': ('Asp', 2), 'GAC': ('Asp', 2),
        'GAA': ('Glu', 2), 'GAG': ('Glu', 2),
        
        # 3-fold orbit
        'AUU': ('Ile', 3), 'AUC': ('Ile', 3), 'AUA': ('Ile', 3),
        
        # 4-fold orbits (standard)
        'GUU': ('Val', 4), 'GUC': ('Val', 4), 'GUA': ('Val', 4), 'GUG': ('Val', 4),
        'CCU': ('Pro', 4), 'CCC': ('Pro', 4), 'CCA': ('Pro', 4), 'CCG': ('Pro', 4),
        'ACU': ('Thr', 4), 'ACC': ('Thr', 4), 'ACA': ('Thr', 4), 'ACG': ('Thr', 4),
        'GCU': ('Ala', 4), 'GCC': ('Ala', 4), 'GCA': ('Ala', 4), 'GCG': ('Ala', 4),
        'GGU': ('Gly', 4), 'GGC': ('Gly', 4), 'GGA': ('Gly', 4), 'GGG': ('Gly', 4),
        
        # 4-fold Gln orbit (ciliate-specific: includes UAA, UAG)
        'CAA': ('Gln', 4), 'CAG': ('Gln', 4), 'UAA': ('Gln', 4), 'UAG': ('Gln', 4),
        
        # 6-fold orbits (standard)
        'UUA': ('Leu', 6), 'UUG': ('Leu', 6), 'CUU': ('Leu', 6), 'CUC': ('Leu', 6), 'CUA': ('Leu', 6), 'CUG': ('Leu', 6),
        'UCU': ('Ser', 6), 'UCC': ('Ser', 6), 'UCA': ('Ser', 6), 'UCG': ('Ser', 6), 'AGU': ('Ser', 6), 'AGC': ('Ser', 6),
        'CGU': ('Arg', 6), 'CGC': ('Arg', 6), 'CGA': ('Arg', 6), 'CGG': ('Arg', 6), 'AGA': ('Arg', 6), 'AGG': ('Arg', 6)
    }
    
    # Convert to DataFrame for fc_checker.py
    orbit_df = pd.DataFrame([
        {'codon': codon, 'amino_acid': aa, 'orbit_size': size, 'orbit_id': f"{aa}{size}"}
        for codon, (aa, size) in orbit_map.items()
    ])
    
    output_path = Path(output_dir)
    orbit_file = output_path / 'orbit_map_ciliate.csv'
    orbit_df.to_csv(orbit_file, index=False)
    
    print(f"✓ Ciliate orbit map created: {orbit_file}")
    
    # Print orbit summary
    orbit_summary = orbit_df.groupby('orbit_size').size()
    print(f"\nCiliate genetic code orbit structure:")
    for size, count in orbit_summary.items():
        print(f"  {size}-fold orbits: {count}")
    
    print(f"\nKey difference from standard code:")
    print(f"  - UAA, UAG: Stop → Gln (converts Stop3 to Gln4)")
    print(f"  - UGA: remains sole Stop (Stop3 → Stop1)")
    print(f"  - Expected FC improvement due to Stop3 elimination")
    
    return orbit_file

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_ciliates.py <codon_file> <index_file> <output_dir>")
        print("Example: python extract_ciliates.py ../CUTG/AGG/genbank_codon_species.tsv ../CUTG/INDEX/genbank_species_index.tsv ciliate_data")
        sys.exit(1)
    
    codon_file, index_file, output_dir = sys.argv[1:4]
    
    # Extract ciliate data
    nuclear_file, mito_file, overlap_species = extract_ciliate_data(codon_file, index_file, output_dir)
    
    # Create ciliate-specific orbit map
    orbit_file = create_ciliate_orbit_map(output_dir)
    
    # Generate test commands
    print(f"\n=== NEXT STEPS: FC TESTING ===")
    print(f"1. Test nuclear ciliates (non-standard code):")
    print(f"   python fc_checker.py {nuclear_file} --map {orbit_file} --progress 50")
    
    print(f"\n2. Test mitochondrial ciliates (standard code):")
    print(f"   python fc_checker.py {mito_file} --map orbit_map_standard.csv --progress 50")
    
    print(f"\n3. Expected FC predictions:")
    print(f"   - Nuclear ciliates: σ_intra/σ_inter < 2.0 (better than standard)")
    print(f"   - Mitochondrial ciliates: σ_intra/σ_inter ≈ 5.0 (similar to your previous result)")
    print(f"   - Smoking gun: Nuclear significantly better than mitochondrial")
    
    if overlap_species:
        print(f"\n4. Paired analysis opportunity:")
        print(f"   - {len(overlap_species)} species have both nuclear and mitochondrial data")
        print(f"   - Perfect for within-organism FC comparison")

if __name__ == "__main__":
    main()