#!/usr/bin/env python3

import os
import re
import argparse
from collections import defaultdict


def safe_filename(name):
    name = name.replace("/", "_").replace("\\", "_")
    name = re.sub(r"\s+", "_", name)
    return name


def wrap_fasta(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def parse_maf(maf_file, remove_gaps=False):
    """
    Returns:
        genome_to_contigs[genome_name] = list of (contig_name, sequence)

    In Parsnp MAF, fields[1] may look like:
        spades_ESBL5780_S5_L001#NODE_1_length_317861_cov_41_164578

    This script groups by only:
        spades_ESBL5780_S5_L001
    """
    genome_to_contigs = defaultdict(list)
    current_cluster = None
    cluster_line_counts = defaultdict(int)

    with open(maf_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            if line.startswith("a "):
                fields = line.split()
                current_cluster = None

                for field in fields:
                    if field.startswith("cluster="):
                        current_cluster = field.split("=", 1)[1]
                        break

                if current_cluster is None:
                    current_cluster = "unknown"

                continue

            if line.startswith("s"):
                fields = line.split()

                if len(fields) < 7:
                    print(f"Warning: skipping malformed s-line:\n{line}")
                    continue

                full_name = fields[1]

                # Keep only genome name before "#"
                genome_name = full_name.split("#", 1)[0]

                # Optional: keep contig info for the FASTA header
                contig_name_original = full_name.split("#", 1)[1] if "#" in full_name else "unknown_contig"

                seq = fields[-1].upper()

                if remove_gaps:
                    seq = seq.replace("-", "")

                cluster_line_counts[(genome_name, current_cluster)] += 1
                count = cluster_line_counts[(genome_name, current_cluster)]

                if count == 1:
                    fasta_contig_name = f"cluster_{current_cluster}|{contig_name_original}"
                else:
                    fasta_contig_name = f"cluster_{current_cluster}_part{count}|{contig_name_original}"

                genome_to_contigs[genome_name].append((fasta_contig_name, seq))

    return genome_to_contigs


def write_fastas(genome_to_contigs, outdir):
    os.makedirs(outdir, exist_ok=True)

    for genome_name, contigs in genome_to_contigs.items():
        fasta_name = safe_filename(genome_name) + ".fna"
        fasta_path = os.path.join(outdir, fasta_name)

        with open(fasta_path, "w") as out:
            for contig_name, seq in contigs:
                out.write(f">{contig_name}\n")
                out.write(wrap_fasta(seq) + "\n")

    print(f"Wrote {len(genome_to_contigs)} genome FASTA files to: {outdir}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract each genome's aligned MAF block sequences into one multi-contig FASTA per genome."
    )
    parser.add_argument(
        "maf",
        help="Input MAF file from Parsnp"
    )
    parser.add_argument(
        "-o", "--outdir",
        default="core_genome",
        help="Output folder name. Default: core_genome"
    )
    parser.add_argument(
        "--remove-gaps",
        action="store_true",
        help="Remove '-' gap characters from aligned sequences before writing FASTA."
    )

    args = parser.parse_args()

    genome_to_contigs = parse_maf(args.maf, remove_gaps=args.remove_gaps)
    write_fastas(genome_to_contigs, args.outdir)


if __name__ == "__main__":
    main()