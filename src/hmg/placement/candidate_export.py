"""
Extract best-hit ORF sequences grouped by protein family.

Protein families are processed in parallel using multiprocessing.
The maximum number of concurrent processes is configurable.

Each worker:
- Reads the ORF FASTA sequentially.
- Extracts sequences belonging to one protein family.
- Writes one protein-specific FASTA file.
- Does not load the entire ORF FASTA into memory.
"""

import logging
import os
import re
import tempfile

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from hmg.search.mmseqs_runner import read_fasta


logger = logging.getLogger(__name__)


def _extract_protein(
    protein: str,
    selected_orfs: set[str],
    orf_fasta: Path,
    output_fasta: Path,
) -> tuple[str, Path, int]:
    """
    Extract candidate ORFs for one protein family.

    This function runs in a separate process.
    """

    found_orfs = set()
    temporary_path = None

    try:
        # Write to a temporary file first.
        # The final output is created only after successful extraction.
        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=output_fasta.parent,
            prefix=f".{output_fasta.stem}_",
            suffix=".tmp",
            delete=False,
        ) as output:

            temporary_path = Path(output.name)

            for orf_id, sequence in read_fasta(orf_fasta):

                if orf_id not in selected_orfs:
                    continue

                if orf_id in found_orfs:
                    raise ValueError(
                        f"Duplicate ORF identifier: {orf_id}"
                    )

                output.write(
                    f">{orf_id}\n{sequence}\n"
                )

                found_orfs.add(orf_id)

                # Stop when all required sequences have been extracted.
                if len(found_orfs) == len(selected_orfs):
                    break

        # Validate extraction
        missing_orfs = selected_orfs - found_orfs

        if missing_orfs:

            preview = sorted(missing_orfs)[:10]

            raise ValueError(
                f"{len(missing_orfs)} ORFs missing for "
                f"protein {protein}: {preview}"
            )

        # Atomically replace the final output file
        os.replace(temporary_path, output_fasta)

        return protein, output_fasta, len(found_orfs)

    except Exception:

        # Remove incomplete temporary output
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

        raise


def export_candidate_sequences(
    best_hits: dict,
    orf_fasta: str | Path,
    output_dir: str | Path,
    workers: int = 4,
) -> dict[str, Path]:
    """
    Extract best-hit ORFs grouped by protein family.

    Protein families are processed in parallel.

    Parameters
    ----------
    best_hits : dict
        Dictionary returned by search_orfs().

    orf_fasta : str | Path
        FASTA file containing all predicted ORFs.

    output_dir : str | Path
        Directory for protein-specific FASTA files.

    workers : int, default=4
        Maximum number of parallel extraction processes.

    Returns
    -------
    dict[str, Path]
        Mapping of protein families to their FASTA files.

        Example:
        {
            "SoxB": Path("candidates/SoxB.faa"),
            "DsrA": Path("candidates/DsrA.faa"),
        }
    """

    orf_fasta = Path(orf_fasta).resolve()
    output_dir = Path(output_dir).resolve()

    # --------------------------------------------------
    # 1. Validate inputs
    # --------------------------------------------------

    if not orf_fasta.is_file():
        raise FileNotFoundError(
            f"ORF FASTA not found: {orf_fasta}"
        )

    if workers < 1:
        raise ValueError(
            "workers must be >= 1"
        )

    if not best_hits:
        logger.warning("No best-hit sequences to export.")
        return {}

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------
    # 2. Group selected ORFs by protein family
    # --------------------------------------------------

    protein_orfs = defaultdict(set)
    assigned_orfs = set()

    for read_id, hit in best_hits.items():

        protein = hit["protein"]
        orf_id = hit["orf_id"]

        if not isinstance(protein, str) or not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9_.-]*", protein
        ):
            raise ValueError(
                f"Invalid protein name: {protein!r}"
            )

        if protein in {".", ".."}:
            raise ValueError(
                f"Invalid protein name: {protein!r}"
            )

        if orf_id in assigned_orfs:
            raise ValueError(
                f"ORF assigned more than once: {orf_id}"
            )

        assigned_orfs.add(orf_id)
        protein_orfs[protein].add(orf_id)

    logger.info(
        "Exporting %d candidate ORFs across %d protein families",
        len(assigned_orfs),
        len(protein_orfs),
    )

    # --------------------------------------------------
    # 3. Parallel processing
    # --------------------------------------------------

    output_files = {}

    max_workers = min(
        workers,
        len(protein_orfs),
    )

    with ProcessPoolExecutor(
        max_workers=max_workers
    ) as executor:

        futures = {}

        for protein, selected_orfs in sorted(
            protein_orfs.items()
        ):

            output_fasta = output_dir / f"{protein}.faa"

            future = executor.submit(
                _extract_protein,
                protein,
                selected_orfs,
                orf_fasta,
                output_fasta,
            )

            futures[future] = protein

        # --------------------------------------------------
        # 4. Collect results
        # --------------------------------------------------

        for future in as_completed(futures):

            protein, output_fasta, count = future.result()

            output_files[protein] = output_fasta

            logger.info(
                "%s: exported %d ORFs to %s",
                protein,
                count,
                output_fasta,
            )

    logger.info(
        "Completed extraction of %d protein families",
        len(output_files),
    )

    return dict(sorted(output_files.items()))