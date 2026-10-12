"""EPA-ng phylogenetic placement using only the Python standard library.

Input: aligned reference FASTA, its matching HMM, query protein FASTA, Newick tree.
External programs: hmmalign (HMMER3), epa-ng.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

LOGGER = logging.getLogger(__name__)


def _fasta_records(path: Path):
    """Yield (identifier, aligned_sequence) without loading all records."""
    identifier, fragments = None, []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(fragments)
                identifier = line[1:].split()[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier in {path}")
                if not re.fullmatch(r"[^\s]+", identifier):
                    raise ValueError(f"Invalid FASTA identifier: {identifier!r}")
                fragments = []
            else:
                if identifier is None:
                    raise ValueError(f"FASTA sequence before header: {path}")
                fragments.append(line)
    if identifier is not None:
        yield identifier, "".join(fragments)


def _read_reference(path: Path):
    sequences = {}
    length = None
    for name, seq in _fasta_records(path):
        seq = seq.upper().replace(".", "-")
        if not seq or any(c.isspace() for c in seq):
            raise ValueError(f"Invalid reference sequence: {name}")
        if name in sequences:
            raise ValueError(f"Duplicate reference ID: {name}")
        if length is None:
            length = len(seq)
        elif len(seq) != length:
            raise ValueError("Reference FASTA is not a multiple sequence alignment")
        sequences[name] = seq
    if len(sequences) < 3:
        raise ValueError("At least three reference sequences required")
    return sequences, length


def _write_stockholm(references: dict[str, str], path: Path):
    """Convert FASTA MSA into the Stockholm format required by --mapali."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        for name, seq in references.items():
            if not re.fullmatch(r"[^\s]+", name):
                raise ValueError(f"Reference ID contains whitespace: {name!r}")
            handle.write(f"{name} {seq}\n")
        handle.write("//\n")


def _read_stockholm(path: Path):
    """Read Stockholm alignment, including interleaved records (no Biopython)."""
    fragments: dict[str, list[str]] = {}
    ended = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped == "//":
                ended = True
                break
            pieces = stripped.split()
            if len(pieces) != 2:
                raise ValueError(f"Unexpected Stockholm alignment line: {stripped!r}")
            name, block = pieces
            fragments.setdefault(name, []).append(block)
    if not ended:
        raise ValueError("Incomplete Stockholm alignment: missing //")
    records = {name: "".join(parts) for name, parts in fragments.items()}
    lengths = {len(sequence) for sequence in records.values()}
    if not records or len(lengths) != 1:
        raise ValueError("Inconsistent Stockholm alignment lengths")
    return records


def _alignment_queries(
    stockholm_path: Path,
    references: dict[str, str],
    reference_length: int,
    output_path: Path,
    expected_query_ids: set[str],
):
    records = _read_stockholm(stockholm_path)
    if not set(references).issubset(records):
        raise ValueError("hmmalign output does not include all reference sequences")
    actual_query_ids = set(records) - set(references)
    if actual_query_ids != expected_query_ids:
        raise ValueError(
            "hmmalign query IDs differ from input: "
            f"missing={list(expected_query_ids - actual_query_ids)[:5]}, "
            f"extra={list(actual_query_ids - expected_query_ids)[:5]}"
        )

    # A column with gaps in every reference is newly inserted for queries.
    # Keep only columns which existed in the original reference alignment.
    # This requires that original reference columns are not all-gap columns.
    aligned_refs = [records[name] for name in references]
    keep = [
        i
        for i in range(len(aligned_refs[0]))
        if any(seq[i] not in "-." for seq in aligned_refs)
    ]
    if len(keep) != reference_length:
        raise ValueError(
            "Cannot preserve seed alignment coordinates: original all-gap "
            "columns or an incompatible HMM/seed alignment."
        )
    for name, original in references.items():
        reconstructed = "".join(records[name][i] for i in keep).upper().replace(".", "-")
        if reconstructed != original:
            raise ValueError(
                f"hmmalign changed reference alignment coordinates for {name}. "
                "Ensure that the HMM was built from exactly this seed alignment."
            )

    with output_path.open("w", encoding="utf-8") as out:
        for name in sorted(actual_query_ids):
            aligned = "".join(records[name][i] for i in keep).upper().replace(".", "-")
            if not aligned.replace("-", ""):
                raise ValueError(f"No aligned residues for query {name}")
            out.write(f">{name}\n{aligned}\n")
    return output_path


def place_sequences(
    tree: str | Path,
    seed_alignment: str | Path,
    hmm: str | Path,
    query_fasta: str | Path,
    model: str,
    output_dir: str | Path,
    threads: int = 8,
) -> Path:
    """Align query proteins against seed HMM and run EPA-ng; return JPLACE path.

    The reference seed alignment is aligned protein FASTA. The profile HMM
    MUST have been built from that identical alignment (hmmalign --mapali).
    All original seed alignment columns must include >=1 reference residue.
    """
    paths = {
        "tree": Path(tree).resolve(),
        "seed_alignment": Path(seed_alignment).resolve(),
        "hmm": Path(hmm).resolve(),
        "query_fasta": Path(query_fasta).resolve(),
    }
    for label, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if not model:
        raise ValueError("EPA-ng model must be supplied")
    for executable in ("hmmalign", "epa-ng"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"Executable not found in PATH: {executable}")

    references, ref_length = _read_reference(paths["seed_alignment"])
    queries = set()
    for name, seq in _fasta_records(paths["query_fasta"]):
        if name in queries or name in references:
            raise ValueError(f"Duplicate/colliding query identifier: {name}")
        if not seq:
            raise ValueError(f"Empty query sequence: {name}")
        queries.add(name)
    if not queries:
        raise ValueError("Query FASTA contains no sequences")

    target = Path(output_dir).resolve()
    target.mkdir(parents=True, exist_ok=True)
    aligned_queries = target / "aligned_queries.fasta"
    jplace_file = target / "epa_result.jplace"
    if jplace_file.exists():
        raise FileExistsError(f"Existing EPA-ng result: {jplace_file}")

    with tempfile.TemporaryDirectory(prefix="hmg_hmmalign_") as temp:
        tmp = Path(temp)
        seed_sto = tmp / "seed.sto"
        combined_sto = tmp / "combined.sto"
        _write_stockholm(references, seed_sto)
        command = [
            "hmmalign", "--amino", "--mapali", str(seed_sto),
            str(paths["hmm"]), str(paths["query_fasta"]),
        ]
        LOGGER.info("Running: %s", " ".join(command))
        with combined_sto.open("w", encoding="utf-8") as out:
            subprocess.run(command, stdout=out, check=True)
        _alignment_queries(combined_sto, references, ref_length, aligned_queries, queries)

    model_arg = str(Path(model).resolve()) if Path(model).is_file() else model
    command = [
        "epa-ng", "--tree", str(paths["tree"]),
        "--ref-msa", str(paths["seed_alignment"]),
        "--query", str(aligned_queries), "--model", model_arg,
        "--threads", str(threads), "--outdir", str(target),
    ]
    LOGGER.info("Running: %s", " ".join(command))
    subprocess.run(command, check=True)
    if not jplace_file.is_file() or not jplace_file.stat().st_size:
        raise RuntimeError(f"EPA-ng did not produce a valid JPLACE file: {jplace_file}")
    return jplace_file