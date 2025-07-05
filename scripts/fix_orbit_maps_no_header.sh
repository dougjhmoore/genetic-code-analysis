#!/bin/bash
# fix_orbit_maps_no_header.sh - Create orbit maps WITHOUT headers for fc_checker.py

echo "🔧 CREATING HEADERLESS ORBIT MAPS"
echo "=================================="

# Create standard orbit map WITHOUT header (64 lines exactly)
echo "Creating standard orbit map (no header)..."
cat > orbit_map_standard_final.csv << 'EOF'
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

# Create ciliate orbit map WITHOUT header (UAA/UAG → Gln = orbit 4)
echo "Creating ciliate orbit map (no header)..."
cat > orbit_map_ciliate_final.csv << 'EOF'
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
UAA,4
UAG,4
UGU,2
UGC,2
UGA,1
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
CAA,4
CAG,4
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

# Verify both files have exactly 64 lines
standard_lines=$(wc -l < orbit_map_standard_final.csv)
ciliate_lines=$(wc -l < orbit_map_ciliate_final.csv)

echo "Standard orbit map: $standard_lines lines"
echo "Ciliate orbit map: $ciliate_lines lines"

if [ "$standard_lines" -eq 64 ] && [ "$ciliate_lines" -eq 64 ]; then
    echo "✅ Both orbit maps have exactly 64 lines (no headers)"
else
    echo "❌ Line count error"
    exit 1
fi

echo ""
echo "🎯 RUNNING CILIATE FC ANALYSIS (FINAL ATTEMPT)"
echo "==============================================="

# Test 1: Ciliate Nuclear - THE CRUCIAL TEST
echo "🧬 CILIATE NUCLEAR FC ANALYSIS"
echo "==============================="
echo "Testing FC prediction: UAA/UAG → Gln should improve compliance"
echo "Expected: σ_intra/σ_inter < 3.0 (significant improvement over 4.809)"
echo ""

python fc_checker.py ciliate_data/ciliate_nuclear_codon_usage.tsv \
    --map orbit_map_ciliate_final.csv \
    --progress 25 \
    --null 1000

echo ""
echo "🧬 CILIATE MITOCHONDRIAL FC ANALYSIS"  
echo "===================================="
echo "Control test with standard genetic code"
echo "Expected: σ_intra/σ_inter ≈ 5.0 (similar to RefSeq mitochondrial)"
echo ""

python fc_checker.py ciliate_data/ciliate_mitochondrial_codon_usage.tsv \
    --map orbit_map_standard_final.csv \
    --progress 50 \
    --null 1000

echo ""
echo "🎉 CILIATE ANALYSIS COMPLETE!"
echo "============================="
echo ""
echo "📊 YOUR COMPLETE FC RESULTS:"
echo "1. Nuclear Standard:      σ_intra/σ_inter = 4.809"
echo "2. Mitochondrial Variant: σ_intra/σ_inter = 4.949" 
echo "3. Ciliate Nuclear:       σ_intra/σ_inter = [see above]"
echo "4. Ciliate Mitochondrial: σ_intra/σ_inter = [see above]"
echo ""
echo "🎯 FC THEORY VALIDATION:"
echo "If Ciliate Nuclear < 3.0: ✅ FC THEORY CONFIRMED!"
echo "If Ciliate Nuclear ≈ 4.8: ⚠️ Theory needs refinement"
echo ""
echo "This completes your empirical validation!"
echo "Ready to integrate into your paper's 'Genetic Code Variants' section."