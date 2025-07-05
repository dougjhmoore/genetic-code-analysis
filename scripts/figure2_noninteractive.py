#!/usr/bin/env python3
"""
figure2_noninteractive.py - Generate Figure 2 panels (headless version)
====================================================================

Non-interactive version that saves all plots without displaying them.
Perfect for WSL/headless environments.

Usage:
    python figure2_noninteractive.py [data_file]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Publication-quality settings
plt.rcParams.update({
    'figure.figsize': (15, 10),
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'serif'
})

# Color scheme from handover document
ORBIT_COLORS = {
    1: '#800080',  # Purple for 1-fold
    2: '#008000',  # Green for 2-fold  
    3: '#FFA500',  # Orange for 3-fold
    4: '#0000FF',  # Blue for 4-fold
    6: '#FF0000',  # Red for 6-fold (Leu/Arg)
    'Ser6': '#FF8C00'  # Orange for Serine 6-fold (distinct from other 6-fold)
}

# RNA codon families for detailed analysis
CODON_FAMILIES = {
    'Met': ['AUG'],
    'Trp': ['UGG'],
    'Phe': ['UUU', 'UUC'],
    'Tyr': ['UAU', 'UAC'],
    'Cys': ['UGU', 'UGC'],
    'His': ['CAU', 'CAC'],
    'Gln': ['CAA', 'CAG'],
    'Asn': ['AAU', 'AAC'],
    'Lys': ['AAA', 'AAG'],
    'Asp': ['GAU', 'GAC'],
    'Glu': ['GAA', 'GAG'],
    'Stop': ['UAA', 'UAG', 'UGA'],
    'Ile': ['AUU', 'AUC', 'AUA'],
    'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Gly': ['GGU', 'GGC', 'GGA', 'GGG'],
    'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
}

def load_codon_data(filepath):
    """Load codon usage data with proper preprocessing."""
    try:
        df = pd.read_csv(filepath, sep='\t', low_memory=False)
        print(f"Loaded {len(df)} species from {filepath}")
        
        # Find codon columns and convert DNA->RNA if needed
        codon_cols = [col for col in df.columns if len(col) == 3]
        
        # Convert to RNA notation
        rename_dict = {}
        for col in codon_cols:
            if 'T' in col:
                rename_dict[col] = col.replace('T', 'U')
        
        if rename_dict:
            df = df.rename(columns=rename_dict)
            codon_cols = [rename_dict.get(col, col) for col in codon_cols]
        
        print(f"Found {len(codon_cols)} codon columns")
        return df, codon_cols
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None

def calculate_family_cvs(df, codon_cols, max_species=5000):
    """Calculate coefficient of variation for each amino acid family."""
    
    # Sample for performance if needed
    if len(df) > max_species:
        df_sample = df.sample(n=max_species, random_state=42)
        print(f"Sampling {max_species} species for analysis...")
    else:
        df_sample = df
    
    family_cv_data = {}
    
    for family_name, codons in CODON_FAMILIES.items():
        available_codons = [c for c in codons if c in codon_cols]
        
        if len(available_codons) > 1:  # Need multiple codons for CV calculation
            cvs = []
            
            for idx, row in df_sample.iterrows():
                # Get codon counts for this family
                counts = []
                for codon in available_codons:
                    if codon in row.index and not pd.isna(row[codon]):
                        counts.append(max(0, row[codon]))
                
                if len(counts) == len(available_codons) and sum(counts) > 0:
                    # Calculate coefficient of variation
                    mean_count = np.mean(counts)
                    std_count = np.std(counts, ddof=1)
                    if mean_count > 0:
                        cv = std_count / mean_count
                        cvs.append(cv)
            
            if cvs:
                family_cv_data[family_name] = cvs
                orbit_size = len(available_codons)
                print(f"{family_name} ({orbit_size}-fold): {len(cvs)} species, median CV = {np.median(cvs):.3f}")
    
    return family_cv_data

def create_figure2a_orbit_violins(family_cv_data, output_path):
    """Create Panel A: Individual orbit violin plots with color gradient."""
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Order families by expected FC compliance (from handover: 2 < 4 ≈ Ser6 < 3 < Leu6/Arg6)
    ordered_families = []
    
    # 2-fold families (best compliance) - Green
    for family in ['Phe', 'Tyr', 'Cys', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu']:
        if family in family_cv_data:
            ordered_families.append((family, 2))
    
    # 4-fold families (good compliance) - Blue  
    for family in ['Val', 'Pro', 'Thr', 'Ala', 'Gly']:
        if family in family_cv_data:
            ordered_families.append((family, 4))
    
    # Serine 6-fold (excellent compliance) - Orange
    if 'Ser' in family_cv_data:
        ordered_families.append(('Ser', 'Ser6'))
    
    # 3-fold families (poor compliance) - Yellow
    for family in ['Stop', 'Ile']:
        if family in family_cv_data:
            ordered_families.append((family, 3))
    
    # Other 6-fold families (worst compliance) - Red
    for family in ['Leu', 'Arg']:
        if family in family_cv_data:
            ordered_families.append((family, 6))
    
    positions = []
    labels = []
    
    for i, (family, orbit_key) in enumerate(ordered_families):
        if family in family_cv_data:
            data = family_cv_data[family]
            
            # Create violin plot
            parts = ax.violinplot([data], positions=[i+1], widths=0.8, 
                                showmeans=True, showmedians=True)
            
            # Color by orbit type
            color = ORBIT_COLORS.get(orbit_key, '#808080')
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
                pc.set_edgecolor('black')
                pc.set_linewidth(0.5)
            
            # Style the statistical markers
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans'):
                if partname in parts:
                    parts[partname].set_edgecolor('black')
                    parts[partname].set_linewidth(1)
            
            positions.append(i+1)
            orbit_size = len(CODON_FAMILIES[family])
            labels.append(f'{family}\n({orbit_size})')
            
            # Add median value as text
            median_cv = np.median(data)
            ax.text(i+1, median_cv + 0.05, f'{median_cv:.2f}', 
                   ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Coefficient of Variation (CV)')
    ax.set_title('First-Classness Compliance by Orbit Family\n(Systematic Gradient: 2-fold < 4-fold ≈ Ser6 < 3-fold < Leu6/Arg6)')
    
    # Add legend
    legend_elements = []
    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=ORBIT_COLORS[2], alpha=0.7, label='2-fold orbits'))
    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=ORBIT_COLORS[4], alpha=0.7, label='4-fold orbits'))
    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=ORBIT_COLORS['Ser6'], alpha=0.7, label='Serine 6-fold'))
    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=ORBIT_COLORS[3], alpha=0.7, label='3-fold orbits'))
    legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=ORBIT_COLORS[6], alpha=0.7, label='Leu/Arg 6-fold'))
    
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output_path / 'Figure2A_orbit_violins.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'Figure2A_orbit_violins.svg', bbox_inches='tight')
    plt.close()  # Close figure to free memory
    print(f"✓ Saved Figure 2A: {output_path / 'Figure2A_orbit_violins.png'}")

def create_figure2b_sixfold_structure(family_cv_data, df, codon_cols, output_path):
    """Create Panel B: Serine vs Leucine fine structure comparison."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Serine subfamilies: UC* (4 codons) vs AG* (2 codons) 
    ser_uc = ['UCU', 'UCC', 'UCA', 'UCG']
    ser_ag = ['AGU', 'AGC']
    
    # Leucine subfamilies: CU* (4 codons) vs UU* (2 codons)
    leu_cu = ['CUU', 'CUC', 'CUA', 'CUG'] 
    leu_uu = ['UUA', 'UUG']
    
    # Sample subset for detailed analysis
    df_subset = df.sample(n=min(1000, len(df)), random_state=42)
    
    # Analyze subfamily balance for Serine
    if 'Ser' in family_cv_data:
        ser_uc_available = [c for c in ser_uc if c in codon_cols]
        ser_ag_available = [c for c in ser_ag if c in codon_cols]
        
        uc_fractions = []
        ag_fractions = []
        
        for _, row in df_subset.iterrows():
            uc_total = sum(row[c] for c in ser_uc_available if c in row.index and not pd.isna(row[c]))
            ag_total = sum(row[c] for c in ser_ag_available if c in row.index and not pd.isna(row[c]))
            total = uc_total + ag_total
            
            if total > 0:
                uc_fractions.append(uc_total / total)
                ag_fractions.append(ag_total / total)
        
        # Box plots for Serine subfamilies
        if uc_fractions and ag_fractions:
            bp1 = ax1.boxplot([uc_fractions, ag_fractions], 
                       labels=['UC* (4 codons)', 'AG* (2 codons)'],
                       patch_artist=True)
            
            # Color the boxes
            bp1['boxes'][0].set_facecolor('lightblue')
            bp1['boxes'][1].set_facecolor('lightcoral')
            
            ax1.set_title('Serine Fine Structure\n(Balanced 4+2 Composition)')
            ax1.set_ylabel('Subfamily Fraction')
            ax1.set_ylim(0, 1)
            
            # Add balance information
            balance_ratio = np.mean(ag_fractions) / np.mean(uc_fractions) if np.mean(uc_fractions) > 0 else 0
            ax1.text(0.5, 0.95, f'Balance Ratio: {balance_ratio:.2f}\nCV = {np.median(family_cv_data["Ser"]):.3f}', 
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # Analyze subfamily balance for Leucine  
    if 'Leu' in family_cv_data:
        leu_cu_available = [c for c in leu_cu if c in codon_cols]
        leu_uu_available = [c for c in leu_uu if c in codon_cols]
        
        cu_fractions = []
        uu_fractions = []
        
        for _, row in df_subset.iterrows():
            cu_total = sum(row[c] for c in leu_cu_available if c in row.index and not pd.isna(row[c]))
            uu_total = sum(row[c] for c in leu_uu_available if c in row.index and not pd.isna(row[c]))
            total = cu_total + uu_total
            
            if total > 0:
                cu_fractions.append(cu_total / total)
                uu_fractions.append(uu_total / total)
        
        # Box plots for Leucine subfamilies
        if cu_fractions and uu_fractions:
            bp2 = ax2.boxplot([cu_fractions, uu_fractions],
                       labels=['CU* (4 codons)', 'UU* (2 codons)'],
                       patch_artist=True)
            
            # Color the boxes
            bp2['boxes'][0].set_facecolor('lightgreen')
            bp2['boxes'][1].set_facecolor('lightyellow')
            
            ax2.set_title('Leucine Fine Structure\n(Imbalanced Geometric Structure)')
            ax2.set_ylabel('Subfamily Fraction')
            ax2.set_ylim(0, 1)
            
            # Add imbalance information
            imbalance_ratio = np.mean(uu_fractions) / np.mean(cu_fractions) if np.mean(cu_fractions) > 0 else 0
            ax2.text(0.5, 0.95, f'Imbalance Ratio: {imbalance_ratio:.2f}\nCV = {np.median(family_cv_data["Leu"]):.3f}', 
                    transform=ax2.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path / 'Figure2B_sixfold_structure.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'Figure2B_sixfold_structure.svg', bbox_inches='tight')
    plt.close()  # Close figure to free memory
    print(f"✓ Saved Figure 2B: {output_path / 'Figure2B_sixfold_structure.png'}")

def create_figure2c_variant_codes(output_path):
    """Create Panel C: Genetic code variant validation bar chart."""
    
    # Data from handover document - proven results
    variants = ['Standard\nNuclear', 'Mitochondrial\nVariant', 'Ciliate\nNuclear']
    fc_ratios = [4.809, 4.949, 5.734]
    colors = ['#4472C4', '#70AD47', '#FFC000']  # Professional color scheme
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    bars = ax.bar(variants, fc_ratios, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
    
    # Add value labels on bars
    for bar, ratio in zip(bars, fc_ratios):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{ratio:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    # Add improvement percentages
    mito_improvement = ((4.949 - 4.809) / 4.809) * 100
    ciliate_improvement = ((5.734 - 4.809) / 4.809) * 100
    
    ax.text(1, 4.949 + 0.2, f'+{mito_improvement:.1f}%', ha='center', va='bottom', 
            fontsize=10, color='green', fontweight='bold')
    ax.text(2, 5.734 + 0.2, f'+{ciliate_improvement:.1f}%', ha='center', va='bottom', 
            fontsize=10, color='green', fontweight='bold')
    
    ax.set_ylabel('σ_intra / σ_inter')
    ax.set_title('Genetic Code Variant Validation\n(FC Compliance Improvements Through Orbit Mergers)')
    ax.set_ylim(0, 6.5)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add explanatory text
    ax.text(0.02, 0.98, 
            'Higher ratios = better FC compliance\nCiliate code shows 19.2% improvement', 
            transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8),
            fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path / 'Figure2C_variant_codes.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_path / 'Figure2C_variant_codes.svg', bbox_inches='tight')
    plt.close()  # Close figure to free memory
    print(f"✓ Saved Figure 2C: {output_path / 'Figure2C_variant_codes.png'}")

def main():
    # Default to the correct path based on your directory structure
    if len(sys.argv) < 2:
        data_file = "../CUTG/AGG/refseq_codon_species.tsv"
        print(f"Using default data file: {data_file}")
    else:
        data_file = sys.argv[1]
    
    output_dir = Path('fc_plots')
    output_dir.mkdir(exist_ok=True)
    
    print("🎯 GENERATING FIGURE 2 - ORBIT ANALYSIS AND FINE STRUCTURE (NON-INTERACTIVE)")
    print("=" * 70)
    
    # Load data
    print("\n📊 Loading codon usage data...")
    df, codon_cols = load_codon_data(data_file)
    
    if df is None:
        print("❌ Failed to load data")
        return
    
    # Calculate family CVs
    print("\n🧮 Calculating coefficient of variation for each family...")
    family_cv_data = calculate_family_cvs(df, codon_cols)
    
    if not family_cv_data:
        print("❌ No CV data calculated")
        return
    
    print(f"\n✅ Successfully calculated CVs for {len(family_cv_data)} families")
    
    # Generate all three panels
    print("\n🎨 Generating Figure 2 panels...")
    
    print("\n📊 Panel A: Orbit violin plots...")
    create_figure2a_orbit_violins(family_cv_data, output_dir)
    
    print("\n🔬 Panel B: Six-fold fine structure...")
    create_figure2b_sixfold_structure(family_cv_data, df, codon_cols, output_dir)
    
    print("\n📈 Panel C: Variant code validation...")
    create_figure2c_variant_codes(output_dir)
    
    print("\n🎉 FIGURE 2 GENERATION COMPLETE!")
    print("=" * 70)
    print(f"📁 All files saved to: {output_dir.absolute()}")
    print("\n📋 Generated files:")
    print("  • Figure2A_orbit_violins.png/svg")
    print("  • Figure2B_sixfold_structure.png/svg") 
    print("  • Figure2C_variant_codes.png/svg")
    print("\n✅ Ready for publication!")

if __name__ == '__main__':
    main()