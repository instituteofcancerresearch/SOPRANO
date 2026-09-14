#!/bin/bash -e
#
# Build the >=30 amino acid transcript length files that SOPRANO's OFF mode
# uses, from the unfiltered ones beside them.
#
# run_localSSBselection_vLOCAL_MOD4OFF.sh substitutes these for the full
# length files when restricting the target BED:
#
#   cut -f1 $BED.tmp | sort -u |
#       fgrep -w -f - $SUPA/ensemble_transcript_protein_min30.length
#
# The two files express one filter in their own units. The protein file is in
# amino acids, so the threshold is 30. The transcript file is in bases, and a
# 30 amino acid protein needs 90 coding bases plus a stop codon, so the
# equivalent threshold is 91. Both drop the same transcripts.
#
# Derived by reproducing Beatriz Monterde's own files exactly: 104763 rows in,
# 102781 out, identical on both sides.
#
# Regenerate after changing Ensembl release -- site counts move with it, see
# the note in src/SOPRANO/core/dnds.py.
#
# Usage: scripts/make_min30_lengths.sh [aux_dir]

AUX="${1:-$(cd "$(dirname "$0")/.." && pwd)/data/aux_soprano}"

MIN_AA=30
MIN_NT=$(( MIN_AA * 3 + 1 ))

for spec in "ensemble_transcript_protein.length:$MIN_AA" \
            "ensemble_transcript.length:$MIN_NT"; do
  src="$AUX/${spec%%:*}"
  thr="${spec##*:}"
  dst="${src%.length}_min30.length"

  [ -f "$src" ] || { echo "Missing $src" >&2; exit 1; }

  awk -F'\t' -v t="$thr" '$2 >= t' "$src" > "$dst"

  printf '%-42s %7s -> %7s rows  (column 2 >= %s)\n' \
    "$(basename "$dst")" "$(wc -l < "$src" | tr -d ' ')" \
    "$(wc -l < "$dst" | tr -d ' ')" "$thr"
done
