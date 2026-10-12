import gzip
import re
import shutil
import subprocess
import tempfile

from pathlib import Path


# OrfM appends:
# _startPosition_frameNumber_orfNumber
ORFM_ID_PATTERN = re.compile(
    r"^(.+)_\d+_[1-6]_\d+$"
)


def get_read_id(orf_id: str) -> str:
    """Recover the original read ID from an OrfM ORF ID."""

    match = ORFM_ID_PATTERN.fullmatch(orf_id)

    if match is None:
        raise ValueError(
            f"Invalid OrfM identifier: {orf_id}"
        )

    return match.group(1)


def read_fasta(path):
    """
    Iterate through a FASTA file without loading
    all sequences into memory.

    Supports uncompressed and gzip-compressed FASTA.
    """

    path = Path(path)

    opener = (
        gzip.open
        if path.suffix == ".gz"
        else open
    )

    with opener(path, "rt") as handle:

        header = None
        sequence_parts = []

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].split()[0]
                sequence_parts = []

            else:
                if header is None:
                    raise ValueError(
                        f"Sequence without FASTA header in {path}"
                    )

                sequence_parts.append(line)

        if header is not None:
            yield header, "".join(sequence_parts)


def prepare_references(
    reference_files,
    output_fasta,
):
    """
    Combine reference FASTA files.

    Each sequence receives a unique internal identifier.
    The original protein family and reference identifier
    are retained in a mapping dictionary.
    """

    reference_mapping = {}
    seen_proteins = set()
    counter = 0

    with open(output_fasta, "w") as output:

        for reference_file in reference_files:

            reference_file = Path(reference_file)

            if not reference_file.is_file():
                raise FileNotFoundError(reference_file)

            # SoxB.faa -> SoxB
            protein = reference_file.name

            if protein.endswith(".gz"):
                protein = protein[:-3]

            protein = Path(protein).stem

            if protein in seen_proteins:
                raise ValueError(
                    f"Duplicate protein family: {protein}"
                )

            seen_proteins.add(protein)

            for original_id, sequence in read_fasta(
                reference_file
            ):

                counter += 1

                internal_id = f"REF{counter:09d}"

                reference_mapping[internal_id] = {
                    "protein": protein,
                    "reference_id": original_id,
                }

                output.write(
                    f">{internal_id}\n{sequence}\n"
                )

    if counter == 0:
        raise ValueError("No reference sequences found.")

    return reference_mapping


def parse_best_hits(
    hits_file,
    reference_mapping,
):
    """
    Parse MMseqs2 results and retain the best-scoring
    alignment across all ORFs belonging to each read.

    The hits file is processed line by line.
    """

    best_hits = {}

    with open(hits_file) as handle:

        for line in handle:

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 4:
                raise ValueError(
                    f"Unexpected MMseqs2 result: {line}"
                )

            orf_id, target_id, bits, evalue = fields

            read_id = get_read_id(orf_id)

            bitscore = float(bits)
            evalue = float(evalue)

            reference = reference_mapping[target_id]

            candidate = {
                "protein": reference["protein"],
                "reference_id": reference["reference_id"],
                "orf_id": orf_id,
                "bitscore": bitscore,
                "evalue": evalue,
            }

            previous = best_hits.get(read_id)

            if previous is None:
                best_hits[read_id] = candidate
                continue

            # Highest bitscore wins.
            # E-value and identifiers resolve ties.
            candidate_rank = (
                -bitscore,
                evalue,
                candidate["protein"],
                candidate["reference_id"],
                orf_id,
            )

            previous_rank = (
                -previous["bitscore"],
                previous["evalue"],
                previous["protein"],
                previous["reference_id"],
                previous["orf_id"],
            )

            if candidate_rank < previous_rank:
                best_hits[read_id] = candidate

    return best_hits


def search_orfs(
    orf_fasta,
    reference_files,
    threads=8,
    sensitivity=5.7,
    max_evalue=1e-5,
    max_hits=10,
):
    """
    Search OrfM-predicted ORFs against multiple protein
    reference sets using MMseqs2.

    Returns
    -------
    dict
        Original read ID -> best-scoring hit information.
    """

    orf_fasta = Path(orf_fasta)

    if not orf_fasta.is_file():
        raise FileNotFoundError(orf_fasta)

    if not reference_files:
        raise ValueError("No reference files provided.")

    executable = shutil.which("mmseqs")

    if executable is None:
        raise RuntimeError(
            "MMseqs2 executable not found in PATH."
        )

    with tempfile.TemporaryDirectory(
        prefix="phylo_mmseqs_"
    ) as temporary_dir:

        temporary_dir = Path(temporary_dir)

        combined_fasta = (
            temporary_dir / "combined_references.faa"
        )

        hits_file = temporary_dir / "hits.tsv"

        tmp_dir = temporary_dir / "tmp"

        # 1. Combine reference sequences
        reference_mapping = prepare_references(
            reference_files,
            combined_fasta,
        )

        # 2. Run MMseqs2
        command = [
            executable,
            "easy-search",
            str(orf_fasta),
            str(combined_fasta),
            str(hits_file),
            str(tmp_dir),
            "--threads",
            str(threads),
            "-s",
            str(sensitivity),
            "-e",
            str(max_evalue),
            "--max-seqs",
            str(max_hits),
            "--format-output",
            "query,target,bits,evalue",
        ]

        subprocess.run(
            command,
            check=True,
        )

        # 3. Select best hit per original read
        best_hits = parse_best_hits(
            hits_file,
            reference_mapping,
        )

    return best_hits