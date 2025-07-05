#!/bin/bash
# debug_ciliate_analysis.sh - Debug and re-run the ciliate analysis

echo "🔍 DEBUGGING CILIATE ANALYSIS"
echo "==============================="

# Check current directory and files
echo "Current directory: $(pwd)"
echo "Available files:"
ls -la

echo ""
echo "Checking ciliate data files:"
if [ -d "ciliate_data" ]; then
    echo "📁 ciliate_data directory exists"
    ls -la ciliate_data/
else
    echo "❌ ciliate_data directory not found"
    exit 1
fi

echo ""
echo "🧪 MANUAL FC ANALYSIS - STEP BY STEP"
echo "===================================="

# Test 1: Ciliate Nuclear (The crucial test!)
echo "Test 1: Ciliate Nuclear FC Analysis"
echo "File: ciliate_data/ciliate_nuclear_codon_usage.tsv"
echo "Orbit map: ciliate_data/orbit_map_ciliate.csv"
echo "Expected: σ_intra/σ_inter < 3.0 (FC improvement predicted)"
echo ""

if [ -f "ciliate_data/ciliate_nuclear_codon_usage.tsv" ] && [ -f "ciliate_data/orbit_map_ciliate.csv" ]; then
    echo "Running ciliate nuclear analysis..."
    python fc_checker.py ciliate_data/ciliate_nuclear_codon_usage.tsv \
        --map ciliate_data/orbit_map_ciliate.csv \
        --progress 25 \
        --null 1000
    echo ""
else
    echo "❌ Required files missing for ciliate nuclear analysis"
fi

# Test 2: Ciliate Mitochondrial (Control test)
echo "Test 2: Ciliate Mitochondrial FC Analysis"
echo "File: ciliate_data/ciliate_mitochondrial_codon_usage.tsv"
echo "Orbit map: orbit_map_standard.csv"
echo "Expected: σ_intra/σ_inter ≈ 5.0 (similar to other mitochondria)"
echo ""

# Create standard orbit map if missing
if [ ! -f "orbit_map_standard.csv" ]; then
    echo "Creating standard orbit map..."
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
fi

if [ -f "ciliate_data/ciliate_mitochondrial_codon_usage.tsv" ] && [ -f "orbit_map_standard.csv" ]; then
    echo "Running ciliate mitochondrial analysis..."
    python fc_checker.py ciliate_data/ciliate_mitochondrial_codon_usage.tsv \
        --map orbit_map_standard.csv \
        --progress 50 \
        --null 1000
    echo ""
else
    echo "❌ Required files missing for ciliate mitochondrial analysis"
fi

echo "📊 COMPARISON SUMMARY"
echo "===================="
echo "Your Previous Results:"
echo "  Nuclear Standard (RefSeq):      σ_intra/σ_inter = 4.809"
echo "  Mitochondrial Variant (RefSeq): σ_intra/σ_inter = 4.949"
echo ""
echo "New Ciliate Results (shown above):"
echo "  Ciliate Nuclear:      [see output above]"
echo "  Ciliate Mitochondrial: [see output above]"
echo ""

echo "🎯 INTERPRETATION GUIDE"
echo "======================="
echo "Ciliate Nuclear Result Interpretation:"
echo "  • If < 2.0: ✅ EXCELLENT FC improvement (strong validation)"
echo "  • If 2.0-3.5: ✅ GOOD FC improvement (moderate validation)"
echo "  • If 3.5-4.5: ⚠️ MARGINAL improvement (weak evidence)"
echo "  • If > 4.5: ❌ NO improvement (theory needs refinement)"
echo ""
echo "Statistical Significance:"
echo "  • Look for 'Null-model Z' and 'empirical P' values"
echo "  • Z > 2 and P < 0.05 indicates significant FC compliance"
echo ""

echo "🔬 WHAT THIS MEANS FOR YOUR PAPER"
echo "=================================="
echo "This analysis tests the KEY FC prediction:"
echo "  'Genetic codes that eliminate problematic orbits show improved compliance'"
echo ""
echo "If ciliate nuclear shows improvement:"
echo "  ✅ Proves FC is architectural principle, not evolutionary accident"
echo "  ✅ Validates mathematical derivation of genetic code structure"
echo "  ✅ Provides quantitative framework for synthetic biology applications"
echo ""
echo "Ready to integrate into your 'Genetic Code Variants' section!"