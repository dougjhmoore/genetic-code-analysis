#!/bin/bash
# fix_orbit_maps.sh - Fix orbit map formatting and run analysis

echo "🔧 FIXING ORBIT MAP FORMATTING"
echo "==============================="

# Create properly formatted standard orbit map (all 64 RNA codons)
echo "Creating standard genetic code orbit map..."
cat > orbit_map_standard_fixed.csv << 'EOF'
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

# Create properly formatted ciliate orbit map (UAA/UAG → Gln)
echo "Creating ciliate genetic code orbit map..."
cat > orbit_map_ciliate_fixed.csv << 'EOF'
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

echo "✅ Fixed orbit maps created"

# Verify the maps have exactly 64 codons
echo "Verifying orbit maps..."
standard_count=$(tail -n +2 orbit_map_standard_fixed.csv | wc -l)
ciliate_count=$(tail -n +2 orbit_map_ciliate_fixed.csv | wc -l)

echo "Standard orbit map: $standard_count codons"
echo "Ciliate orbit map: $ciliate_count codons"

if [ "$standard_count" -eq 64 ] && [ "$ciliate_count" -eq 64 ]; then
    echo "✅ Both orbit maps have exactly 64 codons"
else
    echo "❌ Orbit map formatting error"
    exit 1
fi

echo ""
echo "🧬 RUNNING CILIATE FC ANALYSIS"
echo "=============================="

# Test 1: Ciliate Nuclear (The crucial test!)
echo "🎯 TEST 1: CILIATE NUCLEAR FC ANALYSIS"
echo "======================================"
echo "This tests the KEY FC prediction:"
echo "Genetic codes that eliminate Stop3 orbit should show improved compliance"
echo ""
echo "Genetic code changes in ciliates:"
echo "  UAA: Stop → Glutamine"
echo "  UAG: Stop → Glutamine" 
echo "  UGA: Remains sole stop"
echo ""
echo "Orbit structure changes:"
echo "  Standard: Stop3 {UAA,UAG,UGA} + Gln2 {CAA,CAG}"
echo "  Ciliate:  Stop1 {UGA} + Gln4 {CAA,CAG,UAA,UAG}"
echo ""
echo "FC Prediction: σ_intra/σ_inter should be < 3.0 (major improvement)"
echo ""

if [ -f "ciliate_data/ciliate_nuclear_codon_usage.tsv" ]; then
    echo "Running analysis on 114 ciliate nuclear genomes..."
    python fc_checker.py ciliate_data/ciliate_nuclear_codon_usage.tsv \
        --map orbit_map_ciliate_fixed.csv \
        --progress 25
    
    echo ""
    echo "🔍 CILIATE NUCLEAR RESULT ANALYSIS:"
    echo "If σ_intra/σ_inter < 3.0: ✅ FC THEORY VALIDATED!"
    echo "If σ_intra/σ_inter ≈ 4.8: ⚠️ No improvement (unexpected)"
else
    echo "❌ Ciliate nuclear data file not found"
fi

echo ""
echo "🎯 TEST 2: CILIATE MITOCHONDRIAL FC ANALYSIS"
echo "==========================================="
echo "Control test using standard genetic code"
echo "Expected: σ_intra/σ_inter ≈ 5.0 (similar to your RefSeq mitochondrial result)"
echo ""

if [ -f "ciliate_data/ciliate_mitochondrial_codon_usage.tsv" ]; then
    echo "Running analysis on 254 ciliate mitochondrial genomes..."
    python fc_checker.py ciliate_data/ciliate_mitochondrial_codon_usage.tsv \
        --map orbit_map_standard_fixed.csv \
        --progress 50
else
    echo "❌ Ciliate mitochondrial data file not found"
fi

echo ""
echo "📊 COMPLETE RESULTS COMPARISON"
echo "=============================="
echo "Your Complete FC Analysis Results:"
echo ""
echo "1. Nuclear Standard (RefSeq):      σ_intra/σ_inter = 4.809  (baseline)"
echo "2. Mitochondrial Variant (RefSeq):  σ_intra/σ_inter = 4.949  (slight increase)"
echo "3. Ciliate Nuclear (GenBank):       σ_intra/σ_inter = [see above]"
echo "4. Ciliate Mitochondrial (GenBank): σ_intra/σ_inter = [see above]"
echo ""

echo "🎯 FC THEORY VALIDATION SUMMARY"
echo "==============================="
echo "The ciliate nuclear result is CRUCIAL for your paper because:"
echo ""
echo "✅ If improved (< 3.0): PROVES FC is architectural principle"
echo "   → Genetic code variants that eliminate problematic orbits"
echo "   → show predictable improvements in FC compliance"
echo "   → This validates your mathematical derivation!"
echo ""
echo "⚠️ If no improvement (≈ 4.8): Theory needs refinement"
echo "   → May indicate additional constraints beyond orbit structure"
echo "   → Still valuable scientific insight for future research"
echo ""
echo "📝 READY FOR PAPER INTEGRATION"
echo "==============================" 
echo "Add the ciliate results to your 'Genetic Code Variants' section:"
echo "This provides the quantitative validation of FC architectural predictions"