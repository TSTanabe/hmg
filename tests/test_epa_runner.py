from pathlib import Path
import pytest
from src.hmg.placement.epa_runner import _read_stockholm, _alignment_queries, _read_reference, _write_stockholm


def test_extract_alignment_preserves_seed(tmp_path):
    seed = tmp_path / 'seed.fasta'
    seed.write_text('>A\nAC-G\n>B\nACTG\n>C\nACT-\n')
    refs, n = _read_reference(seed)
    stockholm = tmp_path / 'hmmalign.sto'
    stockholm.write_text('# STOCKHOLM 1.0\nA AC--G\nB ACT-G\nC ACT--\nQ ACW-G\n//\n')
    query = tmp_path / 'queries.fasta'
    _alignment_queries(stockholm, refs, n, query, {'Q'})
    assert query.read_text() == '>Q\nACWG\n'


def test_nonmatching_seed_rejected(tmp_path):
    seed = tmp_path / 'seed.fasta'
    seed.write_text('>A\nAC-G\n>B\nACTG\n>C\nACT-\n')
    refs, n = _read_reference(seed)
    stockholm = tmp_path / 'hmmalign.sto'
    stockholm.write_text('# STOCKHOLM 1.0\nA AG--G\nB ACT-G\nC ACT--\nQ ACW-G\n//\n')
    with pytest.raises(ValueError, match='changed reference'):
        _alignment_queries(stockholm, refs, n, tmp_path/'out.faa', {'Q'})


def test_stockholm_interleaved(tmp_path):
    f = tmp_path/'ali.sto'
    f.write_text('# STOCKHOLM 1.0\nA AC\nB AG\n\nA TT\nB TC\n//\n')
    assert _read_stockholm(f) == {'A':'ACTT','B':'AGTC'}