#!/usr/bin/env bash


failures=0

echo "Testing amino_acid_counter.bash"
../scripts/amino_acid_counter.bash test.fasta

echo "Diffing output. They should be identical."
diff test.amino_acid_counts.tsv test.amino_acid_counts.expected.tsv

failures=`expr $failures + $?`

\rm test.amino_acid_counts.tsv

echo
echo "Failures: $failures"
