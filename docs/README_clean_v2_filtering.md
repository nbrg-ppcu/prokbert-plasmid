# NCBI Plasmid Filtering

This document summarizes the data sources, sequence processing, chromosome-containment QC, and filtering rules used to build the NCBI plasmid database.


## Data Sources And Processing

Representative plasmid-containing assemblies were downloaded as NCBI Datasets genome ZIP packages. 
NCBI genome sequence-report metadata was used to identify plasmid sequence records. Sequence DB rows labelled as plasmid became candidate plasmids; non-plasmid rows from complete-genome assemblies were used to build a chromosome reference for containment QC.

| Item | Count |
|---|---:|
| Downloaded representative plasmid-containing assemblies | 45,185 |
| Original high-confidence plasmid records | 161,655 |
| Clean-v2 retained plasmids | 126,056 |
| Clean-v2 excluded plasmids | 35,599 |
| Clean-v2 FASTA total bp | 9,860,846,760 |
| Clean-v2 host-mapping rows retained | 220,095 |

## Chromosome Reference

The chromosome reference was built from assemblies annotated as `Complete Genome` and `latest` in the NCBI assembly summary. Only non-plasmid sequence DB rows were exported. The minimum chromosome-reference contig length was `200,000 bp`.

| Item | Count |
|---|---:|
| Complete/latest assemblies in metadata | 392,889 |
| Downloaded manifest assemblies | 45,185 |
| Downloaded complete assemblies resolved | 37,602 |
| Downloaded complete assemblies selected | 37,602 |
| Sequence DB rows scanned | 729,036 |
| Non-plasmid chromosome contigs exported | 38,999 |
| Chromosome reference bp | 165,380,469,963 |
| Skipped plasmid-label rows | 100,989 |
| Skipped short non-plasmid rows | 866 |

The chromosome minimap2 index was built from the exported chromosome FASTA using `-x asm5`.

## Terms

| Term | Definition |
|---|---|
| `pid` | Alignment identity: matching bases divided by aligned block length. |
| `qcov` | Fraction of the plasmid query covered after unioning query intervals. |
| `covered_bp` | Non-overlapping plasmid query bp covered by one or more chromosome alignments. |
| `single-contig qcov` | Plasmid query coverage by one chromosome reference contig. |
| `total chromosome qcov` | Plasmid query coverage after aggregating broad chromosome hits. |
| `broad-HQ hit` | Chromosome contig-level hit with `qcov >= 0.80` and `pid >= 0.70`. |
| `HQ recurrent hit` | Chromosome contig-level hit with `qcov >= 0.95` and `pid >= 0.99`. |
| `NR` | skani-derived non-redundant representative database built after clean-v2 filtering. |

## Filtering Rules

Length rules:

| Rule | Action |
|---|---|
| `seq_length < 2,000` | Exclude |
| `seq_length > 4,000,000` | Exclude |
| `2,000 <= seq_length <= 4,000,000` | Retain unless a chromosome-like rule applies |

v1 chromosome-containment rules retained in clean v2:

| Rule source | Action |
|---|---|
| `large_plasmid_chromosomal_categories.tsv` | Exclude all rows in conservative mode |
| `single_assembly_99/97/95/90/80` category files | Exclude |
| `few_large_block_hybrid` category file | Exclude |
| `partial_hybrid_candidate` category file | Exclude |
| `strict_drop_candidates_gt50kb_99.tsv` | Exclude |
| `strict_drop_candidates_gt50kb_99_bacterial.tsv` | Exclude |
| `greyzone_10_50kb_99_manual_review.tsv` | Exclude |

The `large_plasmid_chromosomal_categories.tsv` file was created from the original plasmid-vs-complete-genome-chromosome PAF, restricted to large plasmids `>50 kb` and long high-quality chromosome alignments. A long high-quality alignment was defined as:

```text
PAF alignment block length >= 20,000 bp
pid >= 0.97
```

For each plasmid, query intervals were unioned so repeated or overlapping alignments did not inflate coverage. Coverage was summarized in two ways:

| Metric | Meaning |
|---|---|
| `best_assembly_query_cov_frac` | Best fraction of the plasmid covered by long-HQ hits within one chromosome assembly. |
| `total_query_cov_frac` | Fraction of the plasmid covered by long-HQ chromosome hits across all chromosome assemblies/contigs. |
| `total_merged_query_block_count` | Number of merged plasmid query blocks covered by long-HQ chromosome hits. |

Large-plasmid categorical definitions:

| Category | Definition | Interpretation |
|---|---|---|
| `single_assembly_99` | `best_assembly_query_cov_frac >= 0.99` | Almost the whole plasmid is present in one chromosome assembly. |
| `single_assembly_97` | `best_assembly_query_cov_frac >= 0.97` and `<0.99` | Nearly whole-plasmid chromosome containment in one assembly. |
| `single_assembly_95` | `best_assembly_query_cov_frac >= 0.95` and `<0.97` | Strong near-complete chromosome containment in one assembly. |
| `single_assembly_90` | `best_assembly_query_cov_frac >= 0.90` and `<0.95` | Most of the plasmid is found in one chromosome assembly. |
| `single_assembly_80` | `best_assembly_query_cov_frac >= 0.80` and `<0.90` | Large single-assembly chromosome-like signal. |
| `few_large_block_hybrid` | `total_query_cov_frac >= 0.80` and `total_merged_query_block_count <= 5`, while no single assembly reaches `0.80` | Most of the plasmid is chromosome-like, but spread over a few large blocks/assemblies. |
| `fragmented_chromosomal_signal` | `total_query_cov_frac >= 0.80` and more than `5` merged blocks, while no single assembly reaches `0.80` | Broad but fragmented chromosome signal; exported for review but not part of the clean-v2 retained v1 category-file list. |
| `partial_hybrid_candidate` | `total_query_cov_frac >= 0.50` and `<0.80` | At least half of the plasmid has chromosome-like long-HQ support. |
| `below_50pct_chromosomal_long_hq` | `total_query_cov_frac < 0.50` | Long-HQ chromosome hits exist, but cover less than half of the plasmid. |

Clean v2 used conservative mode for `large_plasmid_chromosomal_categories.tsv`, so every plasmid present in that table was excluded, including the `below_50pct_chromosomal_long_hq` category. The separate category files are retained for auditability and show which category each excluded plasmid belonged to.

The strict/manual chromosome tables were defined from the best plasmid-to-chromosome containment hit table:

| Table | Definition | Interpretation |
|---|---|---|
| `strict_drop_candidates_gt50kb_99.tsv` | Plasmid length `>50 kb`, `qcov >= 0.99`, `pid >= 0.99`, and `query_covered_fraction >= 0.99` against a chromosome reference target. | Large plasmids almost fully contained in a chromosome at very high identity. |
| `strict_drop_candidates_gt50kb_99_bacterial.tsv` | Same as above, additionally requiring the chromosome target to be bacterial. | High-confidence large bacterial chromosome-containment subset. |
| `greyzone_10_50kb_99_manual_review.tsv` | Plasmid length `10-50 kb`, `qcov >= 0.99`, `pid >= 0.99`, and `query_covered_fraction >= 0.99`. | Medium-size elements that are almost perfectly chromosome-contained; biologically ambiguous but excluded in clean v2. |

New clean-v2 chromosome-like rules:

| Rule | Action |
|---|---|
| Long plasmid `>50 kb`, `total_chromosome_qcov >= 0.97`, and at least `20` broad-HQ chromosome contigs | Exclude |
| Long plasmid `>50 kb`, `best_single_contig_qcov >= 0.98`, and at least `20` broad-HQ chromosome contigs | Exclude |
| Long plasmid `>50 kb` present in the loose95 chromosome-containment table | Exclude |
| Medium plasmid `10-50 kb` with at least `5` HQ recurrent chromosome contig hits | Exclude |

## Filtering Results

Overall clean-v2 filtering:

| Metric | Count |
|---|---:|
| Input plasmids | 161,655 |
| Retained plasmids | 126,056 |
| Excluded plasmids, union | 35,599 |
| v1 chromosome exclusion union | 10,463 |
| v2 extra chromosome exclusion union | 1,896 |
| Chromosome exclusion union, v1 or v2 | 11,007 |
| Short `<2 kb` exclusion union | 24,582 |
| Too long `>4 Mb` exclusion union | 21 |
| Excluded with multiple reasons | 3,704 |

Per-rule exclusion counts are not additive because one plasmid can satisfy multiple rules. The `Only this rule` column counts plasmids for which that rule was the sole exclusion reason.

| Rule family | Rule | Plasmids matched by rule | Only this rule |
|---|---|---:|---:|
| Length | `too_short_lt_2000` | 24,582 | 24,582 |
| Length | `too_long_gt_4000000` | 21 | 10 |
| v1 chromosome | `chromosome_category_table_all_rows_conservative` | 8,951 | 6,146 |
| v1 chromosome | `greyzone_10_50kb_99_chromosome_containment` | 1,512 | 960 |
| v1 chromosome | `chromosome_category_partial_hybrid_candidate` | 1,496 | 0 |
| v1 chromosome | `chromosome_category_single_assembly_99` | 759 | 0 |
| v1 chromosome | `strict_gt50kb_99_chromosome_containment` | 509 | 0 |
| v1 chromosome | `strict_gt50kb_99_bacterial_chromosome_containment` | 501 | 0 |
| v1 chromosome | `chromosome_category_single_assembly_80` | 147 | 0 |
| v1 chromosome | `chromosome_category_few_large_block_hybrid` | 113 | 0 |
| v1 chromosome | `chromosome_category_single_assembly_90` | 91 | 0 |
| v1 chromosome | `chromosome_category_single_assembly_97` | 57 | 0 |
| v1 chromosome | `chromosome_category_single_assembly_95` | 30 | 0 |
| v2 extra chromosome | `clean_v2_long_broad_total_qcov97_broadhq20` | 909 | 10 |
| v2 extra chromosome | `clean_v2_long_broad_best_contig_qcov98_broadhq20` | 892 | 0 |
| v2 extra chromosome | `medium_10kb_50kb_recurrent_chromosome_qcov95_pid99_n5` | 678 | 126 |
| v2 extra chromosome | `clean_v2_long_loose95` | 441 | 61 |

Short-plasmid handling:

| Length bin | Input | Retained | Excluded |
|---|---:|---:|---:|
| `<1 kb` | 14,817 | 0 | 14,817 |
| `1-2 kb` | 9,765 | 0 | 9,765 |
| `2-10 kb` | 42,639 | 42,639 | 0 |

v2 extra chromosome-like rule counts:

| Rule | Count |
|---|---:|
| `clean_v2_long_broad_best_contig_qcov98_broadhq20` | 892 |
| `clean_v2_long_broad_total_qcov97_broadhq20` | 909 |
| `clean_v2_long_loose95` | 441 |
| `medium_10kb_50kb_recurrent_chromosome_qcov95_pid99_n5` | 678 |
| v2 extra chromosome union | 1,896 |

Length-bin filtering summary:

| Length bin | Input | Retained | Excluded | v1 chromosome | v2 chromosome | short | too long |
|---|---:|---:|---:|---:|---:|---:|---:|
| `<1 kb` | 14,817 | 0 | 14,817 | 0 | 0 | 14,817 | 0 |
| `1-2 kb` | 9,765 | 0 | 9,765 | 0 | 0 | 9,765 | 0 |
| `2-10 kb` | 42,639 | 42,639 | 0 | 0 | 0 | 0 | 0 |
| `10-50 kb` | 35,922 | 34,284 | 1,638 | 1,512 | 678 | 0 | 0 |
| `50-300 kb` | 51,894 | 43,870 | 8,024 | 7,675 | 1,033 | 0 | 0 |
| `300-500 kb` | 3,559 | 2,592 | 967 | 948 | 70 | 0 | 0 |
| `500 kb-1 Mb` | 1,998 | 1,834 | 164 | 129 | 72 | 0 | 0 |
| `1-4 Mb` | 1,040 | 837 | 203 | 194 | 34 | 0 | 0 |
| `>4 Mb` | 21 | 0 | 21 | 5 | 9 | 0 | 21 |

## NR Database Creation

The clean-v2 NR database is built from the retained clean-v2 all-original plasmids. The NR step does not create new biological filters; it collapses highly similar plasmid sequences into representative clusters while preserving all member-to-assembly evidence.

NR input preparation:

| Input/output | Description | Count |
|---|---|---:|
| `fasta/plasmids.clean_v2.fasta` | Clean-v2 all-original plasmid FASTA; one record per retained `PL...` plasmid ID. | 126,056 |
| `seqrecords/*.fasta` | One FASTA file per plasmid for `skani triangle`. | 126,056 |
| `seqrecords/fasta_files.txt` | List of per-plasmid FASTA files used by skani. | 126,056 |
| `seqrecords/lengths.tsv` | Plasmid ID/file and sequence length table used by clustering. | 126,056 |
| `tables/plasmid_nr_input.feather` | Clean-v2 plasmid metadata used to materialize representative records. | 126,056 |

The skani profile tag is:

```text
skani_minaf80_ani99_afsym97p5_affull99_afpart95
```

Skani is run on the individual FASTA files with:

```bash
skani triangle \
  -l fasta_files.txt \
  -E \
  -t <threads> \
  -c 30 \
  -m 100 \
  --min-af 80 \
  --faster-small \
  -o skani/int_edges__skani_minaf80.tsv
```

The clean-v2 Slurm wrapper used `93` skani threads by default for the triangle step. The `-E` flag writes edge-level pairwise results needed for later filtering. The `--min-af 80` skani prefilter keeps candidate edges with at least moderate alignment fraction before the stricter local edge filtering below.

Skani edge filtering:

| Parameter | Value | Meaning |
|---|---:|---|
| `ANI` | `99.0` | Minimum average nucleotide identity for retained skani edges. |
| `AF_SYM` | `97.5` | Symmetric near-duplicate rule: both directions must have at least this alignment fraction. |
| `AF_FULL` | `99.0` | Containment rule: the larger alignment fraction must be at least this value. |
| `AF_PART` | `95.0` | Containment rule: the smaller alignment fraction must still be at least this value. |

An edge is retained when:

```text
ANI >= 99.0
AND (
  min(AF_ref, AF_query) >= 97.5
  OR (
    max(AF_ref, AF_query) >= 99.0
    AND min(AF_ref, AF_query) >= 95.0
  )
)
```

Retained edges are labeled as:

| Edge reason | Definition |
|---|---|
| `SYMMETRIC` | `min(AF_ref, AF_query) >= 97.5`; sequences are near full-length duplicates of each other. |
| `CONTAINED` | containment rule passes but symmetric rule does not; one sequence is nearly contained in the other. |

If multiple skani rows exist for the same ordered plasmid pair, the best row is retained by preferring `SYMMETRIC` over `CONTAINED`, then higher `minAF`, then higher `ANI`.

Clustering algorithm:

1. Load plasmid lengths from `seqrecords/lengths.tsv`.
2. Orient each retained edge from longer plasmid to shorter plasmid; equal-length ties are broken by plasmid ID.
3. Iterate plasmids from longest to shortest.
4. If a plasmid is not already assigned, make it a representative.
5. Assign its retained shorter neighbors that are still unassigned to that representative.
6. Do not reassign already assigned plasmids.

This makes the representative the longest available plasmid in each local similarity group. The procedure is deterministic for a fixed edge table and length table.
