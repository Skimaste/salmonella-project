#!/usr/bin/env python3
"""
test_snp_detect.py - Test Suite for Reference-Free SNP Detection Pipeline.

Tests:
1. Reverse complement and canonical k-mer symmetry.
2. Synthetic single-nucleotide mutation (SNP) detection.
3. Multi-nucleotide / codon substitution (MNP) detection (mirroring the Salmonella TetR mutation).
4. Error robustness: 1% random Illumina substitution sequencing error simulation.
5. Strand neutrality: reads sampled randomly from forward and reverse strands.
6. Execution against test FASTA files.
"""

import os
import random
import tempfile
import unittest
import snpDetect


def generate_random_genome(length: int = 1500, seed: int = 42) -> str:
    """Generate a random DNA sequence of specified length."""
    random.seed(seed)
    return "".join(random.choices("ACGT", k=length))


def simulate_illumina_reads(genome: str, read_length: int = 150, coverage: int = 30,
                            error_rate: float = 0.0, seed: int = 42) -> list:
    """
    Simulate Illumina reads from a genome:
    - Random uniform starting positions.
    - Double-stranded sequencing (50% forward, 50% reverse complement).
    - Random substitution sequencing errors at rate error_rate.
    """
    random.seed(seed)
    total_bases = len(genome) * coverage
    num_reads = total_bases // read_length
    reads = []
    L = len(genome)

    for i in range(num_reads):
        # Pick random start position
        pos = random.randint(0, L - read_length)
        frag = genome[pos : pos + read_length]

        # 50% chance of reverse complement
        if random.random() < 0.5:
            frag = snpDetect.reverse_complement(frag)

        # Introduce random substitution errors
        if error_rate > 0:
            frag_list = list(frag)
            for j in range(len(frag_list)):
                if random.random() < error_rate:
                    curr = frag_list[j]
                    choices = [b for b in "ACGT" if b != curr]
                    frag_list[j] = random.choice(choices)
            frag = "".join(frag_list)

        reads.append(frag)

    return reads


def write_reads_to_fasta(reads: list, filepath: str):
    """Write list of reads to a FASTA file."""
    with open(filepath, "w", encoding="utf-8") as f:
        for idx, r in enumerate(reads):
            f.write(f">read_{idx}\n{r}\n")


class TestSNPDetect(unittest.TestCase):

    def test_reverse_complement_and_canonical(self):
        """Verify reverse complement and canonical k-mer logic."""
        seq = "ATGCGATC"
        self.assertEqual(snpDetect.reverse_complement(seq), "GATCGCAT")
        # Involutive: rc(rc(x)) == x
        self.assertEqual(snpDetect.reverse_complement(snpDetect.reverse_complement(seq)), seq)

        # Canonical k-mer must always be the lexicographical minimum
        kmer1 = "AAACCC"
        kmer2 = "GGGTTT"
        self.assertEqual(snpDetect.canonical_kmer(kmer1), kmer1)
        self.assertEqual(snpDetect.canonical_kmer(kmer2), kmer1)
        self.assertEqual(snpDetect.canonical_kmer(kmer1), snpDetect.canonical_kmer(kmer2))

    def test_synthetic_clean_snp(self):
        """Verify clean 1-base SNP detection without errors."""
        genome_wt = generate_random_genome(1200, seed=101)
        # Introduce a single SNP at position 600: substitute whatever is there
        mut_pos = 600
        orig_base = genome_wt[mut_pos]
        new_base = "G" if orig_base != "G" else "A"
        genome_var = genome_wt[:mut_pos] + new_base + genome_wt[mut_pos + 1 :]

        with tempfile.TemporaryDirectory() as tmpdir:
            file_wt = os.path.join(tmpdir, "wt.fasta")
            file_var = os.path.join(tmpdir, "var.fasta")

            reads_wt = simulate_illumina_reads(genome_wt, read_length=150, coverage=30, seed=1)
            reads_var = simulate_illumina_reads(genome_var, read_length=150, coverage=30, seed=2)

            write_reads_to_fasta(reads_wt, file_wt)
            write_reads_to_fasta(reads_var, file_var)

            variants = snpDetect.detect_snps(file_wt, file_var, k=25, min_coverage=5, verbose=False)

            self.assertEqual(len(variants), 1, "Expected exactly 1 SNP to be found")
            var = variants[0]

            # The called alleles should match either forward or reverse complement
            forward_match = (var["allele_a"] == orig_base and var["allele_b"] == new_base)
            rc_match = (var["allele_a"] == snpDetect.reverse_complement(orig_base) and
                        var["allele_b"] == snpDetect.reverse_complement(new_base))

            self.assertTrue(forward_match or rc_match,
                            f"Mismatch in called alleles: WT={var['allele_a']}, Var={var['allele_b']}")

    def test_synthetic_codon_substitution(self):
        """Verify 3-base codon substitution (mirroring Salmonella TetR TTT -> GGG)."""
        genome_wt = generate_random_genome(1200, seed=202)
        mut_pos = 500
        # Mutate 3 bases
        wt_triplet = "TTT"
        var_triplet = "GGG"
        genome_wt = genome_wt[:mut_pos] + wt_triplet + genome_wt[mut_pos + 3 :]
        genome_var = genome_wt[:mut_pos] + var_triplet + genome_wt[mut_pos + 3 :]

        with tempfile.TemporaryDirectory() as tmpdir:
            file_wt = os.path.join(tmpdir, "wt.fasta")
            file_var = os.path.join(tmpdir, "var.fasta")

            reads_wt = simulate_illumina_reads(genome_wt, read_length=150, coverage=35, seed=10)
            reads_var = simulate_illumina_reads(genome_var, read_length=150, coverage=35, seed=20)

            write_reads_to_fasta(reads_wt, file_wt)
            write_reads_to_fasta(reads_var, file_var)

            variants = snpDetect.detect_snps(file_wt, file_var, k=25, min_coverage=5, verbose=False)

            self.assertEqual(len(variants), 1, "Expected exactly 1 codon mutation to be detected")
            var = variants[0]

            forward_match = (var["allele_a"] == wt_triplet and var["allele_b"] == var_triplet)
            rc_match = (var["allele_a"] == snpDetect.reverse_complement(wt_triplet) and
                        var["allele_b"] == snpDetect.reverse_complement(var_triplet))

            self.assertTrue(forward_match or rc_match,
                            f"Expected {wt_triplet} -> {var_triplet}, got {var['allele_a']} -> {var['allele_b']}")

    def test_error_robustness(self):
        """Verify robustness against 1% random Illumina substitution errors."""
        genome_wt = generate_random_genome(1500, seed=303)
        mut_pos = 750
        orig_base = genome_wt[mut_pos]
        new_base = "C" if orig_base != "C" else "T"
        genome_var = genome_wt[:mut_pos] + new_base + genome_wt[mut_pos + 1 :]

        with tempfile.TemporaryDirectory() as tmpdir:
            file_wt = os.path.join(tmpdir, "wt_err.fasta")
            file_var = os.path.join(tmpdir, "var_err.fasta")

            # 40x coverage with 1.0% sequencing errors
            reads_wt = simulate_illumina_reads(genome_wt, read_length=150, coverage=40,
                                              error_rate=0.01, seed=30)
            reads_var = simulate_illumina_reads(genome_var, read_length=150, coverage=40,
                                               error_rate=0.01, seed=40)

            write_reads_to_fasta(reads_wt, file_wt)
            write_reads_to_fasta(reads_var, file_var)

            # Auto-estimate coverage threshold
            variants = snpDetect.detect_snps(file_wt, file_var, k=27, verbose=False)

            self.assertEqual(len(variants), 1,
                             f"Expected exactly 1 true SNP and 0 false positives, got {len(variants)}")
            var = variants[0]
            forward_match = (var["allele_a"] == orig_base and var["allele_b"] == new_base)
            rc_match = (var["allele_a"] == snpDetect.reverse_complement(orig_base) and
                        var["allele_b"] == snpDetect.reverse_complement(new_base))
            self.assertTrue(forward_match or rc_match)


if __name__ == "__main__":
    unittest.main()
