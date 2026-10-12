import logging
import shutil
import subprocess

from pathlib import Path


logger = logging.getLogger(__name__)


def predict_orfs(
    input_fastq: str | Path,
    output_fasta: str | Path,
    min_length: int = 96,
    genetic_code: int = 11,
) -> Path:
    """
    Predict and translate ORFs from a FASTQ file using OrfM.

    Supports compressed (.fastq.gz) and uncompressed (.fastq)
    input files.

    OrfM output is streamed directly to a FASTA file.
    Sequence data is never loaded into Python memory.

    Parameters
    ----------
    input_fastq : str | Path
        Input FASTQ file.
    output_fasta : str | Path
        Output FASTA file containing translated ORFs.
    min_length : int
        Minimum ORF length in nucleotides (multiple of 3).
    genetic_code : int
        NCBI translation table.

    Returns
    -------
    Path
        Path to the generated ORF FASTA file.
    """

    input_fastq = Path(input_fastq)
    output_fasta = Path(output_fasta)

    # Validate input
    if not input_fastq.is_file():
        raise FileNotFoundError(
            f"Input FASTQ not found: {input_fastq}"
        )

    if not input_fastq.name.endswith((".fastq", ".fastq.gz")):
        raise ValueError(
            "Input must be .fastq or .fastq.gz"
        )

    if min_length <= 0 or min_length % 3 != 0:
        raise ValueError(
            "Minimum ORF length must be positive and divisible by 3."
        )

    # Check OrfM availability
    executable = shutil.which("orfm")

    if executable is None:
        raise RuntimeError(
            "OrfM executable not found in PATH."
        )

    # Create output directory
    output_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Construct command
    command = [
        executable,
        "-m",
        str(min_length),
        "-c",
        str(genetic_code),
        str(input_fastq),
    ]

    logger.info("Running OrfM on %s", input_fastq)

    # Stream OrfM stdout directly to disk
    with output_fasta.open("wb") as output_handle:

        subprocess.run(
            command,
            stdout=output_handle,
            check=True,
        )

    logger.info(
        "ORF prediction completed: %s",
        output_fasta,
    )

    return output_fasta