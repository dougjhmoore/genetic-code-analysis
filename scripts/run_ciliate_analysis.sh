#!/bin/bash
# run_ciliate_analysis.sh - Execute the complete ciliate FC analysis

echo "🧬 FIRST-CLASSNESS CILIATE ANALYSIS"
echo "==================================="
echo "This analysis tests the crucial FC prediction:"
echo "Genetic codes that eliminate problematic orbits should show improved FC compliance"
echo ""

# Check we're in the right directory
if [ ! -f "fc_checker.py" ]; then
    echo "❌ Error: fc_checker.py not found"
    echo "Please run from the scripts directory containing fc_checker.py"
    exit 1
fi

# Step 1: Extract ciliate data if needed
echo "Step 1: Checking ciliate data..."
if [ ! -d "ciliate_data" ] || [ ! -f "ciliate_data/ciliate_nuclear_codon_usage.tsv" ]; then
    echo "📊 Extracting ciliate data from GenBank..."
    python extract_ciliates.py \
        ../CUTG/AGG/genbank_codon_species.tsv \
        ../CUTG/INDEX/genbank_species_index.tsv \
        ciliate_data
    
    if [ $? -eq 0 ]; then
        echo "✅ Ciliate data extraction successful"
    else
        echo "❌ Ciliate data extraction failed"
        exit 1
    fi
else
    echo "✅ Ciliate data already exists"
fi

# Step 2: Create standard orbit map for comparisons
echo ""
echo "Step 2: Creating standard genetic code orbit map..."
cat > orbit_map_standard.csv << 'EOF'
codon,orbit
UUU,2
UUC,2
UUA,6
UUG,6
UCU,6
UCC,6
UCA,6
UCG,6
UAU,2
UAC,2
UAA,3
UAG,3
UGU,2
UGC,2
UGA,3
UGG,1
CUU,6
CUC,6
CUA,6
CUG,6
CCU,4
CCC,4
CCA,4
CCG,4
CAU,2
CAC,2
CAA,2
CAG,2
CGU,6
CGC,6
CGA,6
CGG,6
AUU,3
AUC,3
AUA,3
AUG,1
ACU,4
ACC,4
ACA,4
ACG,4
AAU,2
AAC,2
AAA,2
AAG,2
AGU,6
AGC,6
AGA,6
AGG,6
GUU,4
GUC,4
GUA,4
GUG,4
GCU,4
GCC,4
GCA,4
GCG,4
GAU,2
GAC,2
GAA,2
GAG,2
GGU,4
GGC,4
GGA,4
GGG,4
EOF
echo "✅ Standard orbit map created"

# Step 3: Verify all required files exist
echo ""
echo "Step 3: Verifying required files..."
required_files=(
    "ciliate_data/ciliate_nuclear_codon_usage.tsv"
    "ciliate_data/ciliate_mitochondrial_codon_usage.tsv" 
    "ciliate_data/orbit_map_ciliate.csv"
    "orbit_map_standard.csv"
    "fc_checker.py"
)

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file"
    else
        echo "❌ Missing: $file"
        exit 1
    fi
done

# Step 4: Run the complete FC analysis
echo ""
echo "🔬 RUNNING COMPLETE FC ANALYSIS"
echo "==============================="

# Create results directory
mkdir -p ciliate_results

echo ""
echo "Analysis 1: Ciliate Nuclear (UAA/UAG → Gln variant)"
echo "FC Prediction: σ_intra/σ_inter < 3.0 (major improvement expected)"
echo "Running analysis on 114 ciliate nuclear genomes..."
python fc_checker.py ciliate_data/ciliate_nuclear_codon_usage.tsv \
    --map ciliate_data/orbit_map_ciliate.csv \
    --progress 50 \
    --null 1000 \
    --plots \
    --out ciliate_results/nuclear \
    --quiet

echo ""
echo "Analysis 2: Ciliate Mitochondrial (standard genetic code)"  
echo "Expected: σ_intra/σ_inter ≈ 5.0 (similar to other mitochondria)"
echo "Running analysis on 254 ciliate mitochondrial genomes..."
python fc_checker.py ciliate_data/ciliate_mitochondrial_codon_usage.tsv \
    --map orbit_map_standard.csv \
    --progress 100 \
    --null 1000 \
    --plots \
    --out ciliate_results/mitochondrial \
    --quiet

echo ""
echo "Analysis 3: Paired Comparison (Within-Organism Test)"
echo "The smoking gun: Nuclear vs mitochondrial within same species"
echo "Running paired analysis on 59 species with both datasets..."

if [ -f "paired_fc_analysis.py" ]; then
    python paired_fc_analysis.py \
        ciliate_data/ciliate_nuclear_codon_usage.tsv \
        ciliate_data/ciliate_mitochondrial_codon_usage.tsv \
        ciliate_data/orbit_map_ciliate.csv \
        orbit_map_standard.csv
else
    echo "⚠️  paired_fc_analysis.py not found - skipping paired analysis"
fi

echo ""
echo "📊 RESULTS SUMMARY"
echo "=================="

# Extract and display key results
echo "Key FC Compliance Results:"

if [ -f "ciliate_results/nuclear/violin_sigma_ratio.png" ]; then
    echo "✅ Ciliate nuclear analysis completed"
    echo "   Results saved to: ciliate_results/nuclear/"
fi

if [ -f "ciliate_results/mitochondrial/violin_sigma_ratio.png" ]; then
    echo "✅ Ciliate mitochondrial analysis completed"  
    echo "   Results saved to: ciliate_results/mitochondrial/"
fi

echo ""
echo "📈 COMPARISON WITH PREVIOUS RESULTS"
echo "==================================="
echo "Previous Results (from your RefSeq analysis):"
echo "  Nuclear Standard:      σ_intra/σ_inter = 4.809  (70,425 organisms)"
echo "  Mitochondrial Variant: σ_intra/σ_inter = 4.949  (441 organisms)"
echo ""
echo "New Ciliate Results:"
echo "  [Results will be displayed above from the fc_checker output]"
echo ""
echo "🎯 FC THEORY VALIDATION"
echo "======================="
echo "If Ciliate Nuclear < Standard Nuclear:"
echo "  ✅ FC THEORY CONFIRMED - genetic code variants improve compliance"
echo "If Ciliate Nuclear ≈ Standard Nuclear:"
echo "  ⚠️  FC theory needs refinement"
echo ""
echo "Next step: Integrate these results into your paper's"
echo "'Genetic Code Variants Validate FC Predictions' section"

echo ""
echo "🏁 ANALYSIS COMPLETE!"
echo "All results saved to ciliate_results/ directory"