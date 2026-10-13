"""Sequence comparisons for a curated multi-protein dataset.

Group summaries describe the selected sequences. They do not establish thermal
adaptation independently of phylogeny, ecology or classification uncertainty.
"""

import csv
import hashlib
import itertools
import json
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import Bio
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParam import ProteinAnalysis

AA = 'ACDEFGHIKLMNPQRSTVWY'
FOCUS = 'GPR'
MATRIX = substitution_matrices.load('BLOSUM62')


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_dataset(metadata_path):
    with metadata_path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('Empty metadata file')
    project_root = Path(__file__).resolve().parent
    seen_ids, seen_sequences = set(), set()
    proteins = []
    for row in rows:
        accession = row['accession']
        if accession in seen_ids or row['group'] not in {'cold', 'meso'}:
            raise ValueError(f'{accession}: duplicate identifier or invalid group')
        path = Path(row['fasta_path'])
        path = path if path.is_absolute() else project_root / path
        records = list(SeqIO.parse(path, 'fasta'))
        if len(records) != 1:
            raise ValueError(f'{accession}: expected one protein')
        sequence = str(records[0].seq).upper()
        sha = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if not sequence or set(sequence) - set(AA):
            raise ValueError(f'{accession}: empty/nonstandard sequence')
        if len(sequence) != int(row['sequence_length']) or sha != row['sequence_sha256']:
            raise ValueError(f'{accession}: sequence differs from validated metadata')
        if sha in seen_sequences:
            raise ValueError(f'{accession}: duplicate sequence')
        seen_ids.add(accession)
        seen_sequences.add(sha)
        proteins.append({**row, 'sequence': sequence, 'input_path': path})
    if any(not any(p['group'] == group for p in proteins) for group in ['cold', 'meso']):
        raise ValueError('Both cold and meso groups are required')
    return proteins


def sequence_features(protein):
    sequence = protein['sequence']
    counts = Counter(sequence)
    analysis = ProteinAnalysis(sequence)
    row = {key: protein.get(key, '') for key in ['accession', 'group', 'organism', 'strain', 'evidence_level', 'thermal_class']}
    row.update(length=len(sequence), estimated_molecular_weight_da=analysis.molecular_weight(), theoretical_pI=analysis.isoelectric_point(), GRAVY=analysis.gravy())
    for aa in AA:
        row[aa + '_count'] = counts[aa]
        row[aa + '_percent'] = 100 * counts[aa] / len(sequence)
    assert sum(row[a + '_count'] for a in AA) == len(sequence)
    return row


def compare_sequences(a, b, window_size=50):
    aligner = Align.PairwiseAligner(mode='global')
    aligner.substitution_matrix = MATRIX
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(a, b)[0]
    idx_a, idx_b = alignment.indices
    columns, variants, windows, focus_changes = [], [], [], Counter()
    counts = Counter()
    gap_runs = 0
    previous_gap_side = ''
    for col, (ia, ib) in enumerate(zip(idx_a, idx_b), start=1):
        aa = a[int(ia)] if ia >= 0 else '-'
        bb = b[int(ib)] if ib >= 0 else '-'
        if aa == '-':
            kind, gap_side = 'gap_in_a', 'a'
        elif bb == '-':
            kind, gap_side = 'gap_in_b', 'b'
        elif aa == bb:
            kind, gap_side = 'match', ''
        else:
            kind, gap_side = 'substitution', ''
        if gap_side and gap_side != previous_gap_side:
            gap_runs += 1
        previous_gap_side = gap_side
        counts[kind] += 1
        score = float(MATRIX[aa, bb]) if not gap_side else None
        focus = kind == 'substitution' and (aa in FOCUS or bb in FOCUS)
        conservative = kind == 'substitution' and score > 0
        if kind == 'substitution':
            counts['positive_score_substitutions' if conservative else 'nonpositive_score_substitutions'] += 1
        if focus:
            counts['GPR_substitutions'] += 1
            focus_changes[(aa, bb)] += 1
        columns.append((aa, bb, kind))
        if kind != 'match':
            variants.append({'alignment_column': col, 'a_position': int(ia) + 1 if ia >= 0 else '', 'b_position': int(ib) + 1 if ib >= 0 else '', 'a_aa': aa, 'b_aa': bb, 'type': kind, 'blosum62_score': score if score is not None else '', 'positive_score_substitution': conservative, 'focus_GPR': focus})
    nongap = counts['match'] + counts['substitution']
    gaps = counts['gap_in_a'] + counts['gap_in_b']
    assert nongap + gaps == len(columns)
    assert sum(x != '-' for x, _, _ in columns) == len(a)
    assert sum(y != '-' for _, y, _ in columns) == len(b)
    for start in range(0, len(columns), window_size):
        block = columns[start:start + window_size]
        paired = [row for row in block if '-' not in row[:2]]
        row = {'start_column': start + 1, 'end_column': start + len(block), 'comparable_sites': len(paired), 'identical_sites': sum(k == 'match' for _, _, k in paired), 'substitutions': sum(k == 'substitution' for _, _, k in paired), 'gap_columns': len(block) - len(paired)}
        row['identity_percent'] = 100 * row['identical_sites'] / len(paired) if paired else ''
        row['substitution_rate'] = row['substitutions'] / len(paired) if paired else ''
        for side, index in [('a', 0), ('b', 1)]:
            residues = [column[index] for column in block if column[index] != '-']
            row[side + '_nongap_residues'] = len(residues)
            for aa in FOCUS:
                row[side + '_' + aa + '_percent'] = 100 * residues.count(aa) / len(residues) if residues else ''
        windows.append(row)
    metrics = {'alignment_score': float(alignment.score), 'alignment_columns': len(columns), 'aligned_non_gap_sites': nongap, 'identical_sites': counts['match'], 'substitutions': counts['substitution'], 'positive_score_substitutions': counts['positive_score_substitutions'], 'nonpositive_score_substitutions': counts['nonpositive_score_substitutions'], 'gap_columns': gaps, 'gap_in_a_columns': counts['gap_in_a'], 'gap_in_b_columns': counts['gap_in_b'], 'gap_runs': gap_runs, 'identity_on_comparable_sites_percent': 100 * counts['match'] / nongap if nongap else None, 'identity_on_all_alignment_columns_percent': 100 * counts['match'] / len(columns), 'GPR_substitution_sites': counts['GPR_substitutions']}
    return metrics, variants, windows, columns, focus_changes


def group_summaries(features, subset):
    rows = []
    for aa in AA:
        groups = {g: [p for p in features if p['group'] == g] for g in ['cold', 'meso']}
        if not all(groups.values()):
            continue
        row = {'subset': subset, 'amino_acid': aa}
        arrays = {}
        for group, members in groups.items():
            values = np.array([m[aa + '_percent'] for m in members])
            arrays[group] = values
            row[group + '_n'] = len(values)
            row[group + '_mean_percent'] = float(values.mean())
            row[group + '_median_percent'] = float(np.median(values))
            row[group + '_sd_percent'] = float(values.std(ddof=1)) if len(values) > 1 else ''
            row[group + '_pooled_residue_percent'] = 100 * sum(m[aa + '_count'] for m in members) / sum(m['length'] for m in members)
        row['difference_pp_cold_minus_meso'] = float(arrays['cold'].mean() - arrays['meso'].mean())
        nc, nm = len(arrays['cold']), len(arrays['meso'])
        if nc > 1 and nm > 1:
            pooled = np.sqrt(((nc - 1) * arrays['cold'].var(ddof=1) + (nm - 1) * arrays['meso'].var(ddof=1)) / (nc + nm - 2))
            row['unadjusted_hedges_g'] = float((1 - 3 / (4 * (nc + nm) - 9)) * row['difference_pp_cold_minus_meso'] / pooled) if pooled else ''
        else:
            row['unadjusted_hedges_g'] = ''
        rows.append(row)
    return rows


def figures(out, features, summaries, identity):
    full = [r for r in summaries if r['subset'] == 'all']
    colors = {'cold': '#3479b7', 'meso': '#d97c39'}
    fig, ax = plt.subplots(figsize=(13, 5))
    x = np.arange(len(AA))
    for group, shift in [('cold', -.2), ('meso', .2)]:
        ax.bar(x + shift, [r[group + '_mean_percent'] for r in full], width=.4, yerr=[r[group + '_sd_percent'] for r in full], capsize=2, label=f'{group} (n={sum(f["group"] == group for f in features)})', color=colors[group])
    ax.set_xticks(x, list(AA))
    ax.set(ylabel='Per-protein amino acid fraction (%)', title='Selected GH2 proteins: group means with SD')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / 'composition_group.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    for ax, aa in zip(axes, FOCUS):
        for pos, group in enumerate(['cold', 'meso']):
            values = [f[aa + '_percent'] for f in features if f['group'] == group]
            offsets = np.linspace(-.12, .12, len(values))
            ax.scatter(pos + offsets, values, color=colors[group], alpha=.85)
            ax.plot([pos - .18, pos + .18], [np.mean(values)] * 2, color='black', lw=2)
        ax.set_xticks([0, 1], ['cold', 'meso'])
        ax.set(title=aa, ylabel='Fraction (%)')
    fig.suptitle('Individual proteins and group means (black lines)')
    fig.tight_layout()
    fig.savefig(out / 'GPR_distribution.png', dpi=180)
    plt.close(fig)
    labels = [f['accession'] for f in features]
    for filename, values, title, color, maximum in [('composition_heatmap.png', np.array([[f[a + '_percent'] for a in AA] for f in features]), 'Amino acid composition per protein (%)', 'YlGnBu', None), ('identity_heatmap.png', identity, 'Pairwise identity at non-gap sites (%)', 'viridis', 100)]:
        fig, ax = plt.subplots(figsize=(11, 10))
        im = ax.imshow(values, aspect='auto', cmap=color, vmin=0, vmax=maximum)
        ax.set_yticks(np.arange(len(labels)), labels)
        ticks = list(AA) if filename.startswith('composition') else labels
        ax.set_xticks(np.arange(len(ticks)), ticks, rotation=90 if filename.startswith('identity') else 0)
        ax.set_title(title)
        split = sum(f['group'] == 'cold' for f in features) - .5
        ax.axhline(split, color='red', lw=1)
        if filename.startswith('identity'):
            ax.axvline(split, color='red', lw=1)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(out / filename, dpi=180)
        plt.close(fig)


def run_dataset_analysis(metadata_path, out, window_size=50):
    sys.stdout.reconfigure(encoding='utf-8')
    proteins = load_dataset(metadata_path)
    proteins.sort(key=lambda p: (p['group'] != 'cold', p['accession']))
    out.mkdir(parents=True, exist_ok=True)
    pair_dir = out / 'pairwise'
    pair_dir.mkdir(exist_ok=True)
    features = [sequence_features(p) for p in proteins]
    write_csv(out / 'features.csv', features)
    summaries = group_summaries(features, 'all')
    gamma_ids = {p['accession'] for p in proteins if 'Gammaproteobacteria' in p.get('phylogenetic_lineage', '')}
    summaries += group_summaries([f for f in features if f['accession'] in gamma_ids], 'gammaproteobacteria_only')
    stronger_ids = {p['accession'] for p in proteins if p['evidence_level'] == 'strain_literature'}
    summaries += group_summaries([f for f in features if f['accession'] in stronger_ids], 'strain_literature_only')
    write_csv(out / 'group_summary.csv', summaries)
    metrics_rows, window_rows, substitution_rows = [], [], []
    identity = np.eye(len(proteins)) * 100
    total = len(proteins) * (len(proteins) - 1) // 2
    for completed, ((i, a), (j, b)) in enumerate(itertools.combinations(enumerate(proteins), 2), start=1):
        pair = a['accession'] + '__' + b['accession']
        pair_type = 'cross_group' if a['group'] != b['group'] else 'within_' + a['group']
        prefix = {'pair_id': pair, 'pair_type': pair_type, 'a_accession': a['accession'], 'b_accession': b['accession'], 'a_group': a['group'], 'b_group': b['group']}
        metrics, variants, windows, columns, changes = compare_sequences(a['sequence'], b['sequence'], window_size)
        metrics_rows.append({**prefix, **metrics})
        identity[i, j] = identity[j, i] = metrics['identity_on_comparable_sites_percent']
        write_csv(pair_dir / (pair + '_variants.csv'), variants, ['alignment_column', 'a_position', 'b_position', 'a_aa', 'b_aa', 'type', 'blosum62_score', 'positive_score_substitution', 'focus_GPR'])
        write_csv(pair_dir / (pair + '_GPR_substitutions.csv'), [v for v in variants if v['focus_GPR']], ['alignment_column', 'a_position', 'b_position', 'a_aa', 'b_aa', 'type', 'blosum62_score', 'positive_score_substitution', 'focus_GPR'])
        aligned_a = ''.join(x for x, _, _ in columns)
        aligned_b = ''.join(y for _, y, _ in columns)
        (pair_dir / (pair + '_alignment.fasta')).write_text(f'>{a["accession"]}\n{aligned_a}\n>{b["accession"]}\n{aligned_b}\n', encoding='utf-8')
        window_rows.extend({**prefix, **row} for row in windows)
        substitution_rows.extend({**prefix, 'a_aa': aa, 'b_aa': bb, 'count': count} for (aa, bb), count in sorted(changes.items()))
        if completed % 20 == 0 or completed == total:
            print(f'Aligned {completed}/{total} pairs', flush=True)
    write_csv(out / 'pairwise_summary.csv', metrics_rows)
    write_csv(out / 'windows.csv', window_rows)
    write_csv(out / 'GPR_substitution_counts.csv', substitution_rows)
    identity_rows = [{'accession': p['accession'], **{q['accession']: float(identity[i, j]) for j, q in enumerate(proteins)}} for i, p in enumerate(proteins)]
    write_csv(out / 'identity_matrix.csv', identity_rows)
    figures(out, features, summaries, identity)
    pair_counts = Counter(row['pair_type'] for row in metrics_rows)
    full = [row for row in summaries if row['subset'] == 'all' and row['amino_acid'] in FOCUS]
    report = {'executed_at_utc': datetime.now(timezone.utc).isoformat(), 'group_counts': dict(Counter(p['group'] for p in proteins)), 'pair_counts': dict(pair_counts), 'total_pairs': total, 'focus_comparison': full, 'parameters': {'alignment_mode': 'global', 'substitution_matrix': 'BLOSUM62', 'open_gap_score': -10, 'extend_gap_score': -.5, 'window_columns': window_size, 'conservative_substitution_definition': 'nonidentical pair with BLOSUM62 score > 0', 'alignment_selection': 'first optimal alignment; local positions may depend on ties'}, 'versions': {'python': platform.python_version(), 'biopython': Bio.__version__, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__}, 'inputs': [{'accession': p['accession'], 'sequence_sha256': p['sequence_sha256'], 'path': str(p['input_path'])} for p in proteins], 'metadata_sha256': hashlib.sha256(metadata_path.read_bytes()).hexdigest(), 'warnings': ['Exploratory selected dataset; phylogeny not statistically corrected.', 'Cold group includes psychrotolerants; multiple Arthrobacter isolates are not ten independent species.', 'Culture temperatures are not growth optima; some classifications have weaker evidence.', 'Sequence-derived mass, pI and GRAVY are estimates.', 'Positive BLOSUM62 scores do not prove preservation of function.', 'Pairwise substitutions have no inferred evolutionary direction.', 'No group significance or causal cold adaptation is established.']}
    (out / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'group_counts': report['group_counts'], 'pair_counts': report['pair_counts'], 'total_pairs': total, 'focus_comparison': full}, ensure_ascii=False, indent=2))

