
import json
from pathlib import Path

import pytest

from hmg.placement.tree_annotation import annotate_tree
from hmg.placement.placement_classification import classify_placements

@pytest.fixture
def example(tmp_path):
    # Root separates the two monophyletic groups; root-adjacent edges 2/5 are boundary.
    tree = '((A:0.1{0},B:0.1{1}):0.2{2},(C:0.1{3},D:0.1{4}):0.2{5});'
    data = {"tree": tree, "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
            "placements": [
                {"n": ["read_in"], "p": [[3 ,-2 ,0.92 ,.01 ,.01] ,[0 ,-3 ,0.08 ,.01 ,.01]]},
                {"n": ["read_out"], "p": [[0 ,-2 ,0.90 ,.01 ,.01] ,[3 ,-3 ,0.10 ,.01 ,.01]]},
                {"nm": [["ambiguous" ,2]], "p": [[2 ,-2 ,0.55 ,.01 ,.01] ,[3 ,-3 ,0.30 ,.01 ,.01]]},
            ], "version": 3}
    jp = tmp_path / 'epa_result.jplace'; jp.write_text(json.dumps(data))
    groups = tmp_path / 'groups.tsv'
    groups.write_text \
        ('reference_id\tgroup\tgroup_name\nA\tOUTGROUP\tOther\nB\tOUTGROUP\tOther\nC\tINGROUP\tSoxB\nD\tINGROUP\tSoxB\n')
    return jp, groups


def test_edges(example):
    jp, groups = example
    result = annotate_tree(jp, groups)
    assert result[0] == {"group": "OUTGROUP", "name": "Other"}
    assert result[3] == {"group": "INGROUP", "name": "SoxB"}
    assert result[2]["group"] == "BOUNDARY"
    assert result[5]["group"] == "BOUNDARY"


def test_placements(example):
    jp, groups = example
    result = classify_placements(jp, groups_file=groups, min_lwr=.75)
    assert result['read_in']['group'] == 'INGROUP'
    assert result['read_out']['group_name'] == 'Other'
    assert result['ambiguous']['group'] == 'UNCERTAIN'
    assert result['ambiguous']['boundary_lwr'] == pytest.approx(.55)
    assert result['ambiguous']['unreported_lwr'] == pytest.approx(.15)


def test_mismatched_groups(example):
    jp, groups = example
    groups.write_text('reference_id\tgroup\tgroup_name\nA\tINGROUP\tSoxB\n')
    with pytest.raises(ValueError, match='mismatch'):
        annotate_tree(jp, groups)

