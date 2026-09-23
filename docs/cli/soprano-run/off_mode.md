# OFF mode

`--off_mode` measures **OFF-target** selection across a cohort of patients.

## Why it exists

Patients have different HLA types, so their immunopeptidomes differ, and so do
the regions that count as OFF-target for each of them. Comparing OFF-target
selection across a cohort needs one region they can all be measured against.

The approach is to take the **complement of the intersection** of the cohort's
immunopeptidomes — everything not presented by every patient — and use that as
the target region. SOPRANO then reports its usual ON-target numbers over that
region, and because of what the region is, **those ON columns are the cohort's
OFF-target estimates**.

Nothing in the output is renamed. The interpretation comes from the target
region you pass, not from the column headers.

## What the switch changes

| | |
|---|---|
| Target BED | restricted to transcripts carrying at least one mutation in the input |
| Transcript lengths | the `_min30` files are used, in place of the unfiltered ones |
| Empty intermediates | abort rather than continue, if the filtered BED or either length file is empty |
| `Pvalue` | reported as `NA` when there is no intron correction |

The `Pvalue` behaviour is deliberately asymmetric. Without intron correction
the ON-versus-OFF comparison is not meaningful in this mode, so no p-value is
reported. With intron correction both p-values are reported as usual.

The `_min30` substitution applies only where you have not named your own
`--transcript` or `--protein_transcript`. A path you pass explicitly is always
used as given. Regenerate the files with `scripts/make_min30_lengths.sh`.

## Usage

Two runs give a cohort both sides. The ON run is an ordinary SOPRANO run
against the merged immunopeptidome; the OFF run adds the switch and targets
the complement of the intersection.

```shell
# ON-target selection
soprano-run \
  -i cohort_ON.anno \
  -b immunopeptidome_merged.bed \
  -o results -n cohort_ON \
  --use_ssb192 --keep_drivers

# OFF-target selection
soprano-run \
  -i cohort_OFF.anno \
  -b immunopeptidome_intersection_complement.bed \
  -o results -n cohort_OFF \
  --use_ssb192 --keep_drivers --off_mode
```

Build the two immunopeptidomes by merging the cohort's individual ones, and by
intersecting them and taking the complement of that intersection.

## Provenance

Ported from `run_localSSBselection_vLOCAL_MOD4OFF.sh` by Beatriz Monterde, on
branch `fix_issue_3` of [luisgls/SOPRANO](https://github.com/luisgls/SOPRANO),
together with its two modified R scripts. That script is
`run_localSSBselection_v4.sh` with 38 substantive lines changed out of 471, so
it is implemented here as a switch on the existing pipeline rather than as a
separate code path.
