import argparse
import logging

from phylo_read_classifier.orf.orfm_runner import predict_orfs


def main():

    parser = argparse.ArgumentParser(
        description="Phylogenetic read classifier"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # ORF prediction
    orf_parser = subparsers.add_parser(
        "orf",
        help="Predict ORFs using OrfM",
    )

    orf_parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input FASTQ or FASTQ.GZ file",
    )

    orf_parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output amino acid FASTA file",
    )

    orf_parser.add_argument(
        "-m",
        "--min-length",
        type=int,
        default=96,
        help="Minimum ORF length in nucleotides",
    )

    orf_parser.add_argument(
        "--genetic-code",
        type=int,
        default=11,
        help="NCBI translation table (default: 11)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    if args.command == "orf":

        predict_orfs(
            input_fastq=args.input,
            output_fasta=args.output,
            min_length=args.min_length,
            genetic_code=args.genetic_code,
        )


if __name__ == "__main__":
    main()