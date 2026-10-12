#!/usr/bin/env python3
"""
snpDetect.py - Reference-Free SNP Detection from Raw Sequencing Reads.

Developed for TATFAR outbreak analysis (Clovis Galiez / Ensimag Course).
Algorithmic complexity: Strictly linear time O(N) in total sequence length.
Memory complexity: Proportional to genome size O(G), filtered via coverage threshold.

Method:
1. Stream reads from both FASTA files and count canonical k-mers.
2. Apply statistical coverage cutoff to filter out sequencing errors (Poisson model).
3. Identify differential solid k-mers (present in A, absent in B, and vice-versa).
4. Resolve sequence bubbles via k-mer pairing and overlap extension.
5. Report mutations with flanking genomic context, allele calls, and read depths.
"""

import os
import sys
import argparse
import time
from collections import Counter
from tqdm import tqdm


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement."""
    rc = reverse_complement(kmer)
    return kmer if kmer < rc else rc


def parse_fasta(filepath: str, max_reads: int = None):
    """
    Generator yielding sequences from a FASTA file in linear time O(N).
    Handles single-line and multi-line FASTA records efficiently.
    """
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        seq_buffer = []
        n_reads = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if seq_buffer:
                    yield "".join(seq_buffer).upper()
                    n_reads += 1
                    if max_reads and n_reads >= max_reads:
                        return
                    seq_buffer = []
            else:
                seq_buffer.append(line)
        if seq_buffer:
            yield "".join(seq_buffer).upper()


def count_canonical_kmers(filepath: str, k: int = 31, max_reads: int = None,
                          desc: str = "Counting k-mers", show_progress: bool = True):
    """
    Count canonical k-mers from a FASTA file in linear time O(N).
    Displays a real-time progress bar using tqdm.
    """
    counts = Counter()
    total_reads = 0
    total_bases = 0

    if max_reads:
        pbar = tqdm(total=max_reads, desc=desc, unit="reads",
                    disable=not show_progress, dynamic_ncols=True)
    else:
        file_size = os.path.getsize(filepath) if os.path.exists(filepath) else 0
        pbar = tqdm(total=file_size, desc=desc, unit="B", unit_scale=True,
                    unit_divisor=1024, disable=not show_progress, dynamic_ncols=True)

    bytes_read = 0
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        seq_buffer = []
        for line in f:
            line_bytes = len(line.encode("utf-8", errors="ignore"))
            bytes_read += line_bytes
            if not max_reads and show_progress:
                pbar.update(line_bytes)

            stripped = line.strip()
            if not stripped:
                continue
            if stripped[0] == ">":
                if seq_buffer:
                    seq = "".join(seq_buffer).upper()
                    total_reads += 1
                    total_bases += len(seq)
                    L = len(seq)
                    if L >= k:
                        for i in range(L - k + 1):
                            kmer = seq[i : i + k]
                            if "N" not in kmer:
                                counts[canonical_kmer(kmer)] += 1
                    if max_reads:
                        if show_progress and total_reads % 500 == 0:
                            pbar.update(500)
                        if total_reads >= max_reads:
                            seq_buffer = []
                            break
                    seq_buffer = []
            else:
                seq_buffer.append(stripped)

        if seq_buffer and (not max_reads or total_reads < max_reads):
            seq = "".join(seq_buffer).upper()
            total_reads += 1
            total_bases += len(seq)
            L = len(seq)
            if L >= k:
                for i in range(L - k + 1):
                    kmer = seq[i : i + k]
                    if "N" not in kmer:
                        counts[canonical_kmer(kmer)] += 1

    pbar.n = total_reads if max_reads else bytes_read
    pbar.refresh()
    pbar.close()
    if show_progress:
        print(f"    Loaded {total_reads:,} reads ({total_bases:,} bp, {len(counts):,} unique canonical k-mers)")
    return counts, total_reads, total_bases


def estimate_coverage_threshold(counts: Counter, user_min_count: int = None) -> int:
    """
    Determine the cutoff threshold t to distinguish genomic k-mers from error k-mers.
    If user_min_count is not provided, estimate it from the k-mer frequency spectrum.
    In Illumina data with >=10x coverage, the spectrum has a trough between
    the low-count error peak (counts 1-3) and the genomic peak.
    """
    if user_min_count is not None and user_min_count > 0:
        return user_min_count

    # Build histogram of multiplicities
    mult_freq = Counter(counts.values())
    # Find local minimum between count 2 and 20
    trough_count = 5
    min_freq = float("inf")
    for c in range(2, min(25, max(mult_freq.keys(), default=2))):
        freq = mult_freq.get(c, 0)
        if freq < min_freq:
            min_freq = freq
            trough_count = c
        elif freq > min_freq * 1.5:
            # Reached ascending edge towards genomic peak
            break

    return max(4, trough_count)


def find_differential_kmers(counts_a: Counter, counts_b: Counter,
                            min_cov_a: int, min_cov_b: int,
                            max_error: int = 1):
    """
    Find solid k-mers specific to A (absent or near-zero in B)
    and solid k-mers specific to B (absent or near-zero in A).
    """
    diff_a = {
        km: c for km, c in counts_a.items()
        if c >= min_cov_a and counts_b.get(km, 0) <= max_error
    }
    diff_b = {
        km: c for km, c in counts_b.items()
        if c >= min_cov_b and counts_a.get(km, 0) <= max_error
    }
    return diff_a, diff_b


def assemble_kmer_chains(kmers_dict: dict, k: int = 31):
    """
    Assemble overlapping differential k-mers into continuous contigs.
    Since differential k-mers form a bubble of length ~ k bases flanking the mutation,
    chaining them reconstructs the local sequence surrounding the mutation.
    """
    # Build prefix/suffix adjacency index
    # We allow forward and reverse complement connections
    all_kmers = set()
    for km in kmers_dict:
        all_kmers.add(km)
        all_kmers.add(reverse_complement(km))

    prefix_map = {}
    suffix_map = {}
    for km in all_kmers:
        pref = km[:-1]
        suff = km[1:]
        prefix_map.setdefault(pref, []).append(km)
        suffix_map.setdefault(suff, []).append(km)

    visited = set()
    contigs = []

    for km in sorted(all_kmers, key=lambda x: kmers_dict.get(canonical_kmer(x), 0), reverse=True):
        if km in visited:
            continue

        # Grow right
        chain = [km]
        visited.add(km)
        visited.add(reverse_complement(km))

        curr = km
        while True:
            target_pref = curr[1:]
            next_kmers = [cand for cand in prefix_map.get(target_pref, [])
                          if cand not in visited]
            if len(next_kmers) == 1:
                nxt = next_kmers[0]
                chain.append(nxt)
                visited.add(nxt)
                visited.add(reverse_complement(nxt))
                curr = nxt
            else:
                break

        # Grow left
        curr = chain[0]
        while True:
            target_suff = curr[:-1]
            prev_kmers = [cand for cand in suffix_map.get(target_suff, [])
                          if cand not in visited]
            if len(prev_kmers) == 1:
                prv = prev_kmers[0]
                chain.insert(0, prv)
                visited.add(prv)
                visited.add(reverse_complement(prv))
                curr = prv
            else:
                break

        # Reconstruct sequence
        assembled = chain[0] + "".join(step[-1] for step in chain[1:])
        avg_depth = sum(kmers_dict.get(canonical_kmer(x), 0) for x in chain) / len(chain)
        contigs.append({
            "seq": assembled,
            "n_kmers": len(chain),
            "avg_depth": avg_depth,
            "kmers": chain
        })

    return contigs


def align_and_call_variants(contigs_a: list, contigs_b: list, k: int = 31):
    """
    Align assembled differential bubble contigs between Strain A and Strain B
    by matching identical flanking anchors on left and right.
    """
    variants = []
    anchor_len = 8

    # Index contigs_b by 8-bp prefix and 8-bp suffix
    b_index = {}
    for cb in contigs_b:
        seq_b = cb["seq"]
        if len(seq_b) < anchor_len:
            continue
        orientations = [(seq_b, False), (reverse_complement(seq_b), True)]
        for cand_b, is_rc in orientations:
            pref = cand_b[:anchor_len]
            b_index.setdefault(pref, []).append((cb, cand_b, is_rc))

    for ca in contigs_a:
        seq_a = ca["seq"]
        if len(seq_a) < anchor_len:
            continue
        pref_a = seq_a[:anchor_len]
        candidates = b_index.get(pref_a, [])

        for cb, cand_b, is_rc in candidates:
            # Find maximum common prefix (left flank)
            prefix_len = 0
            max_prefix = min(len(seq_a), len(cand_b))
            while prefix_len < max_prefix and seq_a[prefix_len] == cand_b[prefix_len]:
                prefix_len += 1

            # Find maximum common suffix (right flank)
            suffix_len = 0
            max_suffix = min(len(seq_a) - prefix_len, len(cand_b) - prefix_len)
            while (suffix_len < max_suffix and
                   seq_a[len(seq_a) - 1 - suffix_len] == cand_b[len(cand_b) - 1 - suffix_len]):
                suffix_len += 1

            allele_a = seq_a[prefix_len : len(seq_a) - suffix_len]
            allele_b = cand_b[prefix_len : len(cand_b) - suffix_len]

            # We require at least some common flanking context (e.g. >= 8 bp)
            if prefix_len >= 8 and suffix_len >= 8 and (allele_a or allele_b):
                left_flank = seq_a[:prefix_len]
                right_flank = seq_a[len(seq_a) - suffix_len:]

                var_entry = {
                    "left_flank": left_flank,
                    "allele_a": allele_a,
                    "allele_b": allele_b,
                    "right_flank": right_flank,
                    "full_seq_a": seq_a,
                    "full_seq_b": cand_b,
                    "depth_a": ca["avg_depth"],
                    "depth_b": cb["avg_depth"],
                    "is_rc": is_rc,
                    "type": "SNP" if len(allele_a) == len(allele_b) == 1 else
                            "MNP" if len(allele_a) == len(allele_b) else
                            "INDEL"
                }

                # Avoid duplicate reporting of same variant in reverse complement
                norm_key = tuple(sorted([
                    (var_entry["allele_a"], var_entry["allele_b"], left_flank[-10:], right_flank[:10]),
                    (reverse_complement(var_entry["allele_a"]),
                     reverse_complement(var_entry["allele_b"]),
                     reverse_complement(right_flank)[:10],
                     reverse_complement(left_flank)[-10:])
                ]))

                if not any(v.get("_norm_key") == norm_key for v in variants):
                    var_entry["_norm_key"] = norm_key
                    variants.append(var_entry)

    return variants


def detect_snps(path_a: str, path_b: str, k: int = 31,
                max_reads: int = None, min_coverage: int = None,
                max_error_count: int = 1, verbose: bool = True):
    """
    Main detection pipeline running in linear time O(N).
    """
    t0 = time.time()
    counts_a, reads_a, bases_a = count_canonical_kmers(
        path_a, k=k, max_reads=max_reads, desc="Strain A (WT)", show_progress=verbose
    )
    counts_b, reads_b, bases_b = count_canonical_kmers(
        path_b, k=k, max_reads=max_reads, desc="Strain B (Var)", show_progress=verbose
    )

    thresh_a = estimate_coverage_threshold(counts_a, min_coverage)
    thresh_b = estimate_coverage_threshold(counts_b, min_coverage)

    if verbose:
        print(f"[*] Coverage cutoffs: Strain A >= {thresh_a}, Strain B >= {thresh_b}")

    diff_a, diff_b = find_differential_kmers(counts_a, counts_b, thresh_a, thresh_b, max_error_count)

    if verbose:
        print(f"[*] Candidate differential k-mers: {len(diff_a)} in Strain A, {len(diff_b)} in Strain B")

    contigs_a = assemble_kmer_chains(diff_a, k=k)
    contigs_b = assemble_kmer_chains(diff_b, k=k)

    if verbose:
        print(f"[*] Assembled {len(contigs_a)} local contig(s) in Strain A, {len(contigs_b)} in Strain B")

    variants = align_and_call_variants(contigs_a, contigs_b, k=k)

    if verbose:
        print(f"[*] Found {len(variants)} variant locus/loci in {time.time() - t0:.2f}s total execution time.\n")

    return variants


def format_report(variants: list, name_a: str = "Strain A", name_b: str = "Strain B") -> str:
    """Format called variants into a human-readable report."""
    if not variants:
        return "No confident variants detected between the two datasets."

    lines = []
    lines.append("=" * 80)
    lines.append(f"VARIATION CALLING REPORT: {name_a} (Wild-Type) vs {name_b} (Variant)")
    lines.append("=" * 80)

    for idx, var in enumerate(variants, start=1):
        lines.append(f"\n[Variant {idx}] Type: {var['type']}")
        lines.append(f"  {name_a} Allele: {var['allele_a']} (mean read depth: {var['depth_a']:.1f}x)")
        lines.append(f"  {name_b} Allele: {var['allele_b']} (mean read depth: {var['depth_b']:.1f}x)")

        l_ctx = var["left_flank"][-25:] if len(var["left_flank"]) > 25 else var["left_flank"]
        r_ctx = var["right_flank"][:25] if len(var["right_flank"]) > 25 else var["right_flank"]

        lines.append(f"\n  Genomic Context (Forward orientation):")
        lines.append(f"    {name_a}: ... {l_ctx} [{var['allele_a']}] {r_ctx} ...")
        lines.append(f"    {name_b}: ... {l_ctx} [{var['allele_b']}] {r_ctx} ...")

        # Also display reverse complement
        rc_allele_a = reverse_complement(var["allele_a"])
        rc_allele_b = reverse_complement(var["allele_b"])
        rc_l = reverse_complement(r_ctx)
        rc_r = reverse_complement(l_ctx)

        lines.append(f"  Complementary Strand:")
        lines.append(f"    {name_a}: ... {rc_l} [{rc_allele_a}] {rc_r} ...")
        lines.append(f"    {name_b}: ... {rc_l} [{rc_allele_b}] {rc_r} ...")
        lines.append("-" * 80)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Reference-free de novo SNP and variant detector from Illumina reads."
    )
    parser.add_argument("readsA", help="Path to FASTA reads of Strain A (e.g. Wild-type)")
    parser.add_argument("readsB", help="Path to FASTA reads of Strain B (e.g. Variant)")
    parser.add_argument("-k", "--kmer", type=int, default=31, help="k-mer length (default: 31)")
    parser.add_argument("-m", "--max-reads", type=int, default=None,
                        help="Maximum reads to sample from each file (e.g. 250000 for ultra-fast mode)")
    parser.add_argument("-c", "--min-coverage", type=int, default=None,
                        help="Coverage cutoff threshold (default: auto-detected)")
    parser.add_argument("-e", "--max-error", type=int, default=1,
                        help="Maximum count allowed in opposite strain to tolerate sequencing noise (default: 1)")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Optional file path to write results")

    args = parser.parse_args()

    variants = detect_snps(
        args.readsA,
        args.readsB,
        k=args.kmer,
        max_reads=args.max_reads,
        min_coverage=args.min_coverage,
        max_error_count=args.max_error,
        verbose=True
    )

    report = format_report(variants, name_a="Wild-Type", name_b="Variant")
    print(report)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report + "\n")
        print(f"\nReport saved to: {args.output}")


if __name__ == "__main__":
    main()
