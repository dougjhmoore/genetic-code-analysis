#!/bin/bash
# extract_ciliate_nuclear_result.sh - Get the crucial ciliate nuclear FC result

echo "🔍 EXTRACTING CILIATE NUCLEAR RESULT"
echo "====================================="

echo "Running ONLY the ciliate nuclear analysis to capture the result clearly..."
echo ""
echo "🎯 CILIATE NUCLEAR FC ANALYSIS"
echo "=============================="
echo "This is THE crucial test for FC theory:"
echo "Genetic code variants that eliminate Stop3 orbit should show improved compliance"
echo ""
echo "Ciliate changes: UAA/UAG Stop → Gln (creates Gln4 orbit, eliminates Stop3)"
echo "FC Prediction: σ_intra/σ_inter should be significantly < 4.809"
echo ""

# Run just the ciliate nuclear analysis
python fc_checker.py ciliate_data/ciliate_nuclear_codon_usage.tsv \
    --map orbit_map_ciliate_final.csv \
    --progress 25 \
    --quiet

echo ""
echo "📊 COMPLETE RESULTS SUMMARY"
echo "==========================="
echo "1. Nuclear Standard (RefSeq):      σ_intra/σ_inter = 4.809  (70,425 organisms)"
echo "2. Mitochondrial Variant (RefSeq):  σ_intra/σ_inter = 4.949  (441 organisms)"  
echo "3. Ciliate Mitochondrial (GenBank): σ_intra/σ_inter = 5.044  (254 organisms)"
echo "4. Ciliate Nuclear (GenBank):       σ_intra/σ_inter = [SEE ABOVE]"

echo ""
echo "🎯 FC THEORY VALIDATION ANALYSIS"
echo "================================"
echo ""
echo "✅ CONFIRMED PREDICTIONS:"
echo "  • Ciliate mitochondrial (5.044) ≈ RefSeq mitochondrial (4.949) ✓"
echo "  • Mitochondrial constraint effects consistent across datasets ✓"
echo ""
echo "🔬 KEY TEST RESULT:"
echo "  • Ciliate Nuclear = [VALUE FROM ABOVE]"
echo ""
echo "INTERPRETATION:"
echo "  • If < 2.0:  🎉 EXCELLENT FC validation (major breakthrough!)"
echo "  • If 2.0-3.5: ✅ GOOD FC validation (solid evidence)"
echo "  • If 3.5-4.5: ⚠️ MARGINAL improvement (weak evidence)" 
echo "  • If > 4.5:  ❌ NO improvement (theory needs refinement)"
echo ""
echo "This result determines whether FC is a fundamental organizing principle!"