"""Annotate JPLACE tree edges by ingroup/outgroup, without Biopython.

The input tree must be rooted at the INGROUP/OUTGROUP split.
Root-adjacent edges are conservatively classified as BOUNDARY.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class _Node:
    name: str | None = None
    edge: int | None = None
    children: list['_Node'] = field(default_factory=list)
    descendants: frozenset[str] = frozenset()


class _NewickParser:
    """Parse the subset of Newick used by EPA-ng/JPLACE (including {edge})."""

    def __init__(self, source: str):
        self.s = source
        self.i = 0

    def _ws(self):
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    def _token(self, stops: str) -> str:
        start = self.i
        while self.i < len(self.s) and self.s[self.i] not in stops:
            self.i += 1
        return self.s[start:self.i].strip()

    def _label(self) -> str | None:
        self._ws()
        if self.i < len(self.s) and self.s[self.i] == "'":
            self.i += 1
            parts = []
            while self.i < len(self.s):
                if self.s[self.i] == "'":
                    self.i += 1
                    if self.i < len(self.s) and self.s[self.i] == "'":
                        parts.append("'")
                        self.i += 1
                    else:
                        return ''.join(parts)
                else:
                    parts.append(self.s[self.i])
                    self.i += 1
            raise ValueError('Unterminated quoted Newick label')
        result = self._token(':,();{}[]')
        return result or None

    def _comment(self):
        # Standard Newick comments aren't edge identifiers.
        if self.s[self.i] != '[':
            return
        self.i += 1
        level = 1
        while self.i < len(self.s) and level:
            if self.s[self.i] == '[':
                level += 1
            elif self.s[self.i] == ']':
                level -= 1
            self.i += 1
        if level:
            raise ValueError('Unterminated Newick comment')

    def _node(self) -> _Node:
        self._ws()
        children = []
        if self.i < len(self.s) and self.s[self.i] == '(':
            self.i += 1
            while True:
                children.append(self._node())
                self._ws()
                if self.i >= len(self.s):
                    raise ValueError('Unexpected end of Newick')
                sep = self.s[self.i]
                self.i += 1
                if sep == ')':
                    break
                if sep != ',':
                    raise ValueError(f'Expected comma or closing parenthesis, got {sep!r}')
        name = self._label()
        edge = None
        self._ws()
        if self.i < len(self.s) and self.s[self.i] == ':':
            self.i += 1
            length = self._token(',();{}[]')
            try:
                float(length)
            except ValueError as exc:
                raise ValueError(f'Invalid branch length: {length!r}') from exc
        self._ws()
        while self.i < len(self.s) and self.s[self.i] in '{[':
            if self.s[self.i] == '[':
                self._comment()
            else:
                self.i += 1
                number = self._token('}')
                if self.i >= len(self.s) or self.s[self.i] != '}':
                    raise ValueError('Unterminated JPLACE edge identifier')
                self.i += 1
                if edge is not None or not number.isdecimal():
                    raise ValueError(f'Invalid or duplicated edge annotation: {number!r}')
                edge = int(number)
            self._ws()
        if not children and not name:
            raise ValueError('Unnamed leaf in Newick tree')
        return _Node(name, edge, children)

    def parse(self) -> _Node:
        root = self._node()
        self._ws()
        if self.i >= len(self.s) or self.s[self.i] != ';':
            raise ValueError('Newick must end with a semicolon')
        self.i += 1
        self._ws()
        if self.i != len(self.s):
            raise ValueError('Trailing content after Newick tree')
        return root


def _load_groups(path: str | Path) -> dict[str, tuple[str, str]]:
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'reference_id', 'group', 'group_name'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'groups.tsv must contain columns: {sorted(required)}')
        groups = {}
        for row in reader:
            ref = row['reference_id'].strip()
            kind = row['group'].strip().upper()
            name = row['group_name'].strip()
            if not ref or kind not in {'INGROUP', 'OUTGROUP'} or not name:
                raise ValueError(f'Invalid group annotation: {row}')
            if ref in groups:
                raise ValueError(f'Duplicate reference_id: {ref}')
            groups[ref] = (kind, name)
    if not groups:
        raise ValueError('No reference annotations found')
    return groups


def annotate_tree(jplace_file: str | Path, groups_file: str | Path) -> dict[int, dict]:
    """Return edge_id -> {'group': INGROUP|OUTGROUP|BOUNDARY, 'name': str|None}.

    Requires a rooted, bifurcating INGROUP/OUTGROUP split and monophyletic
    named groups. Root-adjacent edges are always marked BOUNDARY.
    """
    with Path(jplace_file).open(encoding='utf-8') as handle:
        data = json.load(handle)
    newick = data.get('tree')
    if not isinstance(newick, str) or not newick.strip():
        raise ValueError("JPLACE file lacks a valid 'tree' string")
    root = _NewickParser(newick).parse()
    groups = _load_groups(groups_file)

    if len(root.children) != 2:
        raise ValueError('Expected a bifurcating rooted tree with two root children')

    nodes = []
    tips = []

    def visit(node: _Node) -> frozenset[str]:
        nodes.append(node)
        if not node.children:
            tips.append(node.name)
            node.descendants = frozenset({node.name})
        else:
            desc = set()
            for child in node.children:
                desc.update(visit(child))
            node.descendants = frozenset(desc)
        return node.descendants

    visit(root)
    if len(tips) != len(set(tips)):
        raise ValueError('Tree tip names must be unique')
    if set(tips) != set(groups):
        raise ValueError(
            f'Tip/group mismatch: missing annotations={sorted(set(tips)-set(groups))[:10]}, '
            f'extra annotations={sorted(set(groups)-set(tips))[:10]}'
        )
    split = [set(groups[t][0] for t in node.descendants) for node in root.children]
    if {frozenset(s) for s in split} != {frozenset({'INGROUP'}), frozenset({'OUTGROUP'})}:
        raise ValueError('Root must separate all INGROUP tips from all OUTGROUP tips')

    for kind, name in set(groups.values()):
        designated = frozenset(t for t, value in groups.items() if value == (kind, name))
        if not any(node.descendants == designated for node in nodes):
            raise ValueError(f'Group {kind}/{name} is not monophyletic')

    result = {}
    for node in nodes:
        if node.edge is None:
            if node is root:
                continue
            raise ValueError('JPLACE tree contains an unnumbered edge')
        if node.edge in result:
            raise ValueError(f'Duplicate JPLACE edge number: {node.edge}')
        annotations = {groups[t] for t in node.descendants}
        if node is root or node in root.children or len(annotations) != 1:
            result[node.edge] = {'group': 'BOUNDARY', 'name': None}
        else:
            kind, name = next(iter(annotations))
            result[node.edge] = {'group': kind, 'name': name}
    if not result:
        raise ValueError('No numbered JPLACE edges found')
    return result