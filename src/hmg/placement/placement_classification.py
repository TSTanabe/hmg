
"""Classify JPLACE query placements using weighted edge annotations."""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

from .tree_annotation import annotate_tree


def classify_placements(
        jplace_file: str | Path,
        edge_annotations: dict[int, dict] | None = None,
        *,
        groups_file: str | Path | None = None,
        min_lwr: float = 0.75,
) -> dict[str, dict]:
    """Classify queries into INGROUP, named OUTGROUP, or UNCERTAIN.

    `min_lwr` tests the absolute sum of stored LWRs for a group, without
    renormalizing missing placement mass. Boundary mass is kept separate.
    When several named outgroups exist, each must independently exceed the
    threshold. Does not infer organism taxonomy.
    """
    if not (0 < min_lwr <= 1):
        raise ValueError("min_lwr must lie within (0, 1]")
    if edge_annotations is None:
        if groups_file is None:
            raise ValueError("Provide edge_annotations or groups_file")
        edge_annotations = annotate_tree(jplace_file, groups_file)
    with Path(jplace_file).open(encoding="utf-8") as handle:
        data = json.load(handle)
    fields = data.get("fields", [])
    required = {"edge_num", "like_weight_ratio"}
    if len(fields) != len(set(fields)) or not required.issubset(fields):
        raise ValueError(f"JPLACE fields must uniquely include {sorted(required)}")
    edge_i, lwr_i = (fields.index(key) for key in ("edge_num", "like_weight_ratio"))
    output: dict[str, dict] = {}
    for pquery in data.get("placements", []):
        if ("n" in pquery) == ("nm" in pquery):
            raise ValueError("Each pquery must contain exactly one of 'n' or 'nm'")
        ids = pquery["n"] if "n" in pquery else [entry[0] for entry in pquery["nm"]]
        weight_by_name: dict[str, float] = defaultdict(float)
        boundary_weight = 0.0
        total = 0.0
        best_edge = None
        best_weight = -1.0
        for row in pquery.get("p", []):
            if len(row) != len(fields):
                raise ValueError("Placement row length does not match JPLACE fields")
            edge = int(row[edge_i])
            weight = float(row[lwr_i])
            if not math.isfinite(weight) or weight < 0 or weight > 1:
                raise ValueError(f"Invalid LWR: {weight}")
            if edge not in edge_annotations:
                raise ValueError(f"Unknown edge {edge} in annotations")
            total += weight
            if weight > best_weight:
                best_weight, best_edge = weight, edge
            annotation = edge_annotations[edge]
            group = annotation["group"]
            if group == "BOUNDARY":
                boundary_weight += weight
            elif group in {"INGROUP", "OUTGROUP"}:
                name = annotation["name"]
                if not name:
                    raise ValueError(f"Missing group name for edge {edge}")
                weight_by_name[(group, name)] += weight
            else:
                raise ValueError(f"Unknown group label: {group}")
        if total > 1 + 1e-5:
            raise ValueError(f"Stored LWRs sum to >1: {total}")
        ranked = sorted(weight_by_name.items(), key=lambda kv: (-kv[1], kv[0]))
        selected = ranked[0] if ranked and ranked[0][1] >= min_lwr else None
        # Distinct ingroup names should be treated separately rather than pooled.
        classification = selected[0][0] if selected else "UNCERTAIN"
        name = selected[0][1] if selected else None
        record = {
            "group": classification,
            "group_name": name,
            "ingroup_lwr": sum(v for (kind, _), v in weight_by_name.items() if kind == "INGROUP"),
            "outgroup_lwr": sum(v for (kind, _), v in weight_by_name.items() if kind == "OUTGROUP"),
            "outgroup_lwr_by_name": {n: v for (k, n), v in weight_by_name.items() if k == "OUTGROUP"},
            "boundary_lwr": boundary_weight,
            "stored_lwr": total,
            "unreported_lwr": max(0.0, 1.0 - total),
            "best_edge": best_edge,
            "best_edge_lwr": best_weight if best_edge is not None else None,
        }
        for query_id in ids:
            if query_id in output:
                raise ValueError(f"Duplicate query identifier: {query_id}")
            output[query_id] = dict(record)
    return output

