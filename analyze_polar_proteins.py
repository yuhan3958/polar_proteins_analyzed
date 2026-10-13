import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

AA = 'ACDEFGHIKLMNPQRSTVWY'
FOCUS = 'GPR'

def read_one(path):
    records = list(SeqIO.parse(str(path), 'fasta'))
    if len(records) != 1:
        raise ValueError(f'{path}: FASTA 파일에는 단백질 한 개만 있어야 합니다 (현재 {len(records)}개).')
    rec = records[0]
    seq = str(rec.seq).upper().rstrip('*')
    bad = sorted(set(seq) - set(AA))
    if bad:
        raise ValueError(f'{path}: 표준 아미노산이 아닌 문자 발견: {bad}')
    if not seq:
        raise ValueError(f'{path}: 서열이 비어 있습니다.')
    return rec, seq

def save_csv(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def main():
    p = argparse.ArgumentParser(description='남극 저온성 vs 대장균 β-galactosidase 분석')
    p.add_argument('--metadata', help='다중 단백질 분석에 사용할 metadata.csv')
    p.add_argument('--cold', default='P81650.fasta')
    p.add_argument('--meso', default='P00722.fasta')
    p.add_argument('--out', default=None)
    p.add_argument('--window', type=int, default=50, help='서열 정렬 창 크기 (열 단위)')
    args = p.parse_args()
    if args.window < 5: p.error('--window는 5 이상이어야 합니다.')
    if args.metadata:
        from dataset_analysis import run_dataset_analysis
        run_dataset_analysis(Path(args.metadata), Path(args.out or 'dataset_results'), args.window)
        return
    args.out = args.out or 'polar_results'
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    cold_rec, cold = read_one(Path(args.cold))
    meso_rec, meso = read_one(Path(args.meso))
    seqs = {'P81650_cold': cold, 'P00722_meso': meso}
    cnt = {k: Counter(s) for k, s in seqs.items()}
    pct = {k: {a: 100 * cnt[k][a] / len(s) for a in AA} for k, s in seqs.items()}

    comp = [{'amino_acid': a, 'cold_count': cnt['P81650_cold'][a],
             'meso_count': cnt['P00722_meso'][a],
             'cold_percent': round(pct['P81650_cold'][a], 4),
             'meso_percent': round(pct['P00722_meso'][a], 4),
             'difference_pp': round(pct['P81650_cold'][a]-pct['P00722_meso'][a], 4)} for a in AA]
    save_csv(out/'composition.csv', list(comp[0]), comp)

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(20)
    ax.bar(x-0.2, [pct['P81650_cold'][a] for a in AA], width=.4, label='P81650 (cold)')
    ax.bar(x+0.2, [pct['P00722_meso'][a] for a in AA], width=.4, label='P00722 (meso)')
    ax.set_xticks(x, list(AA)); ax.set_ylabel('Fraction of residues (%)')
    ax.set_title('Amino acid composition'); ax.legend(); fig.tight_layout()
    fig.savefig(out/'composition.png', dpi=180); plt.close(fig)

    aligner = Align.PairwiseAligner(mode='global')
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(cold, meso)[0]
    idx_a, idx_b = alignment.indices
    matches = substitutions = gap_columns = compared = 0
    variants = []
    aligned_lines = []
    window_rows = []
    align_chars_a, align_chars_b, marks = [], [], []
    # 각 대응 위치를 실제 alignment column에서 검토한다.
    for col, (ia, ib) in enumerate(zip(idx_a, idx_b), start=1):
        a = cold[int(ia)] if ia >= 0 else '-'
        b = meso[int(ib)] if ib >= 0 else '-'
        align_chars_a.append(a); align_chars_b.append(b)
        marks.append('|' if a == b else (' ' if '-' in (a,b) else '.'))
        if '-' in (a,b):
            gap_columns += 1
            kind = 'gap'
        else:
            compared += 1
            if a == b:
                matches += 1
                kind = 'match'
            else:
                substitutions += 1
                kind = 'substitution'
        if kind != 'match':
            variants.append({'alignment_column': col,
                             'cold_position': int(ia)+1 if ia >= 0 else '',
                             'meso_position': int(ib)+1 if ib >= 0 else '',
                             'cold_aa': a, 'meso_aa': b, 'type': kind,
                             'focus_GPR': (a in FOCUS or b in FOCUS) and kind == 'substitution'})
    save_csv(out/'variants.csv', ['alignment_column','cold_position','meso_position',
                                 'cold_aa','meso_aa','type','focus_GPR'], variants)
    # 정렬은 1차 서열의 상동 위치 추정이며 구조와 기능을 증명하지 않는다.
    with (out/'alignment.txt').open('w', encoding='utf-8') as f:
        for i in range(0,len(align_chars_a),80):
            f.write(f'COLD  {i+1:5}: {"".join(align_chars_a[i:i+80])}\n')
            f.write(f'             {"".join(marks[i:i+80])}\n')
            f.write(f'MESO  {i+1:5}: {"".join(align_chars_b[i:i+80])}\n\n')

    for start in range(0,len(idx_a),args.window):
        suba = align_chars_a[start:start+args.window]
        subb = align_chars_b[start:start+args.window]
        paired = [(a,b) for a,b in zip(suba,subb) if a!='-' and b!='-']
        differences = sum(a!=b for a,b in paired)
        window_rows.append({'start_column':start+1,
                            'end_column': min(start+args.window,len(idx_a)),
                            'comparable_sites':len(paired), 'substitutions': differences,
                            'substitution_rate':round(differences/len(paired),5) if paired else ''})
    save_csv(out/'windows.csv', list(window_rows[0]), window_rows)
    fig, ax = plt.subplots(figsize=(12,4.5))
    valid = [r for r in window_rows if r['substitution_rate'] != '']
    ax.plot([(r['start_column']+r['end_column'])/2 for r in valid],
            [r['substitution_rate'] for r in valid], marker='o', markersize=3)
    ax.set(xlabel='Alignment column', ylabel='Substitution fraction (non-gap sites)',
           title=f'Sequence differences by {args.window}-column window', ylim=(0,1))
    fig.tight_layout(); fig.savefig(out/'windows.png', dpi=180); plt.close(fig)

    focus = [v for v in variants if v['focus_GPR']]
    save_csv(out/'GPR_substitutions.csv', list(variants[0]) if variants else
             ['alignment_column','cold_position','meso_position','cold_aa','meso_aa','type','focus_GPR'], focus)
    stats = {
        'cold': {'id': cold_rec.id, 'length':len(cold)},
        'meso': {'id': meso_rec.id, 'length':len(meso)},
        'alignment': {'columns':len(idx_a), 'aligned_non_gap_sites':compared,
                      'identical_sites':matches,'substitutions':substitutions,
                      'gap_columns':gap_columns,
                      'identity_on_comparable_sites_percent':round(100*matches/compared,3) if compared else None,
                      'identity_on_all_alignment_columns_percent':round(100*matches/len(idx_a),3)},
        'focus_percent': {a:{'cold':round(pct['P81650_cold'][a],4),
                             'meso':round(pct['P00722_meso'][a],4),
                             'difference_pp':round(pct['P81650_cold'][a]-pct['P00722_meso'][a],4)} for a in FOCUS},
        'focus_substitution_sites':len(focus),
        'warnings':['두 단백질만으로 통계적 유의성 또는 저온 적응의 인과성을 판정할 수 없습니다.',
                    '서열 정렬은 구조적 활성 부위 위치를 직접 검증하지 않습니다.',
                    'BLOSUM62 전역 정렬은 하나의 가정이며 결과는 갭 벌점에 영향을 받습니다.']
    }
    (out/'summary.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(stats,ensure_ascii=False,indent=2))
    print(f'\n결과 폴더: {out.resolve()}')

if __name__ == '__main__':
    main()
