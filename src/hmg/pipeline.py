"""Run the HMG metagenomic ORF-to-placement workflow for one sample."""

import json
import logging
import os
from pathlib import Path

from hmg.orf.orfm_runner import predict_orfs
from hmg.search.mmseqs_runner import search_orfs
from hmg.placement.candidate_export import export_candidate_sequences
from hmg.placement.epa_runner import place_sequences
from hmg.placement.tree_annotation import annotate_tree
from hmg.placement.placement_classification import classify_placements

logger = logging.getLogger(__name__)

REQUIRED_REFERENCE_FILES = (
    "references.faa",
    "tree.nwk",
    "seed_alignment.fasta",
    "reference.hmm",
    "groups.tsv",
)


def discover_reference_packages(reference_dir: str | Path) -> dict[str, dict[str, Path]]:
    """Find and validate all protein-specific reference packages."""
    root = Path(reference_dir).resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)

    packages = {}
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        missing = [name for name in REQUIRED_REFERENCE_FILES if not (folder / name).is_file()]
        if missing:
            raise ValueError(f"Incomplete reference package {folder.name}: {missing}")
        packages[folder.name] = {name: folder / name for name in REQUIRED_REFERENCE_FILES}

    if not packages:
        raise ValueError(f"No reference packages found in {root}")
    return packages


def prepare_search_references(
    packages: dict[str, dict[str, Path]], output_dir: str | Path
) -> list[Path]:
    """Create protein-named symlinks because MMseqs2 module infers family from filename."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for protein, files in packages.items():
        source = files["references.faa"].resolve()
        target = output_dir / f"{protein}.faa"
        if target.is_symlink():
            if target.resolve() != source:
                target.unlink()
        elif target.exists():
            raise FileExistsError(f"Expected symlink, found existing file: {target}")
        if not target.exists():
            target.symlink_to(source)
        paths.append(target)
    return paths


def run_pipeline(
    input_fastq: str | Path,
    reference_dir: str | Path,
    output_dir: str | Path,
    *,
    threads: int = 8,
    extraction_workers: int = 4,
    min_orf_length: int = 96,
    genetic_code: int = 11,
    sensitivity: float = 5.7,
    max_evalue: float = 1e-5,
    max_hits: int = 10,
    model: str = "LG+G4",
    min_lwr: float = 0.75,
) -> dict:
    """Run OrfM, MMseqs2, candidate extraction, EPA-ng, and group classification."""
    input_fastq = Path(input_fastq).resolve()
    output_dir = Path(output_dir).resolve()
    if not input_fastq.is_file():
        raise FileNotFoundError(input_fastq)
    if threads < 1 or extraction_workers < 1:
        raise ValueError("threads and extraction_workers must be >= 1")

    packages = discover_reference_packages(reference_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Discovered %d reference packages", len(packages))

    reference_fastas = prepare_search_references(
        packages, output_dir / "reference_links"
    )

    # 1. ORF prediction: subprocess streams amino acid FASTA to disk.
    orf_fasta = predict_orfs(
        input_fastq=input_fastq,
        output_fasta=output_dir / "01_orfs" / "all_orfs.faa",
        min_length=min_orf_length,
        genetic_code=genetic_code,
    )

    # 2. MMseqs2 search: one combined search across all protein families.
    best_hits = search_orfs(
        orf_fasta=orf_fasta,
        reference_files=reference_fastas,
        threads=threads,
        sensitivity=sensitivity,
        max_evalue=max_evalue,
        max_hits=max_hits,
    )
    logger.info("Best hits for %d reads", len(best_hits))

    # Persist the assignments for auditability.
    search_dir = output_dir / "02_search"
    search_dir.mkdir(parents=True, exist_ok=True)
    best_hits_path = search_dir / "best_hits.json"
    with best_hits_path.open("w") as handle:
        json.dump(best_hits, handle, indent=2)

    # 3. Extract only selected ORF sequences; one FASTA per protein family.
    candidate_files = export_candidate_sequences(
        best_hits=best_hits,
        orf_fasta=orf_fasta,
        output_dir=output_dir / "03_candidates",
        workers=extraction_workers,
    )

    # 4-5. Each protein uses its own tree, alignment, HMM, and group definitions.
    results = {}
    for protein, query_fasta in sorted(candidate_files.items()):
        if protein not in packages:
            raise KeyError(f"No reference package for protein: {protein}")
        package = packages[protein]
        placement_dir = output_dir / "04_placement" / protein
        logger.info("Placing %s sequences", protein)

        jplace_file = place_sequences(
            tree=package["tree.nwk"],
            seed_alignment=package["seed_alignment.fasta"],
            hmm=package["reference.hmm"],
            query_fasta=query_fasta,
            model=model,
            output_dir=placement_dir,
            threads=threads,
        )
        edge_annotations = annotate_tree(
            jplace_file=jplace_file,
            groups_file=package["groups.tsv"],
        )
        classifications = classify_placements(
            jplace_file=jplace_file,
            edge_annotations=edge_annotations,
            min_lwr=min_lwr,
        )
        results[protein] = {
            "jplace_file": str(jplace_file),
            "classifications": classifications,
        }

        classification_dir = output_dir / "05_classification"
        classification_dir.mkdir(parents=True, exist_ok=True)
        with (classification_dir / f"{protein}.json").open("w") as handle:
            json.dump(classifications, handle, indent=2)

    return {
        "orf_fasta": str(orf_fasta),
        "best_hits_file": str(best_hits_path),
        "candidate_files": {key: str(path) for key, path in candidate_files.items()},
        "proteins": results,
    }

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "HMG: Functional and phylogenetic classification "
            "of metagenomic reads using OrfM, MMseqs2 and EPA-ng."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --------------------------------------------------
    # Required arguments
    # --------------------------------------------------

    parser.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        help="Input metagenome in FASTQ or FASTQ.GZ format",
    )

    parser.add_argument(
        "-r", "--references",
        required=True,
        type=Path,
        help="Directory containing protein reference packages",
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output directory",
    )

    # --------------------------------------------------
    # General parameters
    # --------------------------------------------------

    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads for MMseqs2 and EPA-ng",
    )

    parser.add_argument(
        "--extraction-workers",
        type=int,
        default=4,
        help="Number of parallel candidate extraction processes",
    )

    # --------------------------------------------------
    # OrfM parameters
    # --------------------------------------------------

    parser.add_argument(
        "--min-orf-length",
        type=int,
        default=96,
        help="Minimum ORF length in nucleotides",
    )

    parser.add_argument(
        "--genetic-code",
        type=int,
        default=11,
        help="NCBI genetic code for ORF translation",
    )

    # --------------------------------------------------
    # MMseqs2 parameters
    # --------------------------------------------------

    parser.add_argument(
        "--sensitivity",
        type=float,
        default=5.7,
        help="MMseqs2 search sensitivity",
    )

    parser.add_argument(
        "--max-evalue",
        type=float,
        default=1e-5,
        help="Maximum MMseqs2 E-value",
    )

    parser.add_argument(
        "--max-hits",
        type=int,
        default=10,
        help="Maximum number of MMseqs2 hits per ORF",
    )

    # --------------------------------------------------
    # EPA-ng parameters
    # --------------------------------------------------

    parser.add_argument(
        "--model",
        type=str,
        default="LG+G4",
        help="Amino acid substitution model for EPA-ng",
    )

    # --------------------------------------------------
    # Classification parameters
    # --------------------------------------------------

    parser.add_argument(
        "--min-lwr",
        type=float,
        default=0.75,
        help="Minimum summed likelihood weight ratio for classification",
    )

    args = parser.parse_args()

    # --------------------------------------------------
    # Validate parameters
    # --------------------------------------------------

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    if args.extraction_workers < 1:
        parser.error("--extraction-workers must be >= 1")

    if args.max_hits < 1:
        parser.error("--max-hits must be >= 1")

    if args.max_evalue <= 0:
        parser.error("--max-evalue must be > 0")

    if not 0 <= args.min_lwr <= 1:
        parser.error("--min-lwr must be between 0 and 1")

    if not args.input.is_file():
        parser.error(f"Input FASTQ not found: {args.input}")

    if not args.references.is_dir():
        parser.error(
            f"Reference directory not found: {args.references}"
        )

    # --------------------------------------------------
    # Logging
    # --------------------------------------------------

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    # --------------------------------------------------
    # Run pipeline
    # --------------------------------------------------

    try:
        results = run_pipeline(
            input_fastq=args.input,
            reference_dir=args.references,
            output_dir=args.output,
            threads=args.threads,
            extraction_workers=args.extraction_workers,
            min_orf_length=args.min_orf_length,
            genetic_code=args.genetic_code,
            sensitivity=args.sensitivity,
            max_evalue=args.max_evalue,
            max_hits=args.max_hits,
            model=args.model,
            min_lwr=args.min_lwr,
        )

    except Exception:
        logging.exception("HMG pipeline failed")
        raise SystemExit(1)

    logging.info("Pipeline completed successfully")
    logging.info("ORF FASTA: %s", results["orf_fasta"])
    logging.info("Best hits: %s", results["best_hits_file"])
    logging.info(
        "Protein families processed: %d",
        len(results["proteins"]),
    )


if __name__ == "__main__":
    main()