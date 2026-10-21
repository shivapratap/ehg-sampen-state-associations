"""Render the manuscript display items (Fig. 2-4, S1, Table 1 inputs) from saved CSVs.

    python scripts/build_figure_inputs.py          # 1) build the small adapter tables
    python scripts/make_manuscript_figures.py --overwrite   # 2) render PDF/PNG/LaTeX

Render reviewed frozen CSVs. No fitting, FDR recalculation or data selection by p/q.

Input paths resolve relative to the JSON configuration's directory.
All enabled items validate before any publication output is written.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd

COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#6A6A6A']
MARKERS = ['o', 's', '^', 'D', 'v']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def style():
    font_manager.findfont('DejaVu Sans', fallback_to_default=False)
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.labelsize': 9, 'axes.titlesize': 10,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.linewidth': 0.7,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none',
        'axes.prop_cycle': plt.cycler(color=COLORS),
    })


def numeric(df, columns):
    for col in columns:
        require(col in df, f'missing numeric column: {col}')
        df[col] = pd.to_numeric(df[col], errors='raise')
        require(np.isfinite(df[col]).all(), f'nonfinite {col}; report/resolve undefined cases explicitly')


def prepare(item, config_dir):
    require(item['kind'] in ['forest', 'gain', 'agreement', 'table'], 'unknown kind')
    require(item.get('caption'), 'caption is required')
    require('DRAFT:' not in item['caption'], 'replace the draft caption with a verified caption before enabling')
    require('REPLACE_' not in json.dumps(item), 'resolve configuration placeholders before enabling')
    src = (config_dir / item['source']).resolve()
    require(src.is_file(), f'missing source: {src}')
    df = pd.read_csv(src)
    source_rows = len(df)
    df = df.rename(columns=item.get('rename', {}))
    require(not df.columns.duplicated().any(), 'duplicate columns after rename')
    for col, mapping in item.get('value_maps', {}).items():
        require(col in df, f'value-map column absent: {col}')
        df[col] = df[col].map(lambda x: mapping.get(str(x), x))
    for col, vals in item.get('filters', {}).items():
        require(col in df, f'filter column absent: {col}')
        df = df[df[col].isin(vals)].copy()
    require(len(df) > 0, 'no rows after declared filters')
    for col, vals in item.get('expected', {}).items():
        require(col in df, f'expected metadata column absent: {col}')
        require(not df[col].isna().any(), f'empty metadata: {col}')
        require(set(df[col]) == set(vals), f'{col} values differ from declared expected values')
    if 'expected_rows' in item:
        require(len(df) == item['expected_rows'], 'row count differs from expected_rows')
    kind = item['kind']
    if kind == 'forest':
        cols = ['label', 'group', 'estimate', 'ci_low', 'ci_high', 'contrast',
                'scale', 'interval_method', 'family_id', 'family_size', 'n_women', 'n_windows']
        require(set(cols) <= set(df), f'forest missing: {set(cols)-set(df)}')
        require(not df[cols].isna().any().any(), 'empty required forest fields')
        for col in ['contrast', 'scale', 'interval_method', 'family_id', 'family_size']:
            require(col in item.get('expected', {}), f'declare expected {col} to prevent silent mixing')
            require(df[col].nunique() == 1, f'mixed {col}: use separate panels/config items')
        numeric(df, ['estimate', 'ci_low', 'ci_high', 'n_women', 'n_windows', 'family_size'])
        require((df.n_women > 0).all() and (df.n_windows >= df.n_women).all(), 'invalid denominators')
        require(all((df[c] % 1 == 0).all() for c in ['n_women', 'n_windows', 'family_size']), 'counts must be integers')
        require(((df.ci_low <= df.estimate) & (df.estimate <= df.ci_high)).all(), 'interval does not bracket estimate')
        require(not df.duplicated(['label', 'group']).any(), 'duplicate forest label/group')
        require(set(item['order']) == set(df.label), 'order must list every retained label exactly once')
        require(len(item['order']) == len(set(item['order'])), 'duplicate label order')
        require(set(item['groups']) == set(df.group), 'groups must enumerate all groups')
        require(all(set(g.label) == set(item['order']) for _, g in df.groupby('group')), 'missing feature/group cell; separate or document unavailable comparisons')
        require(item.get('xlabel'), 'forest xlabel must specify units/contrast')
    elif kind == 'gain':
        require({'window_id', 'group', 'gain', 'value'} <= set(df), 'gain needs window_id,group,gain,value')
        require(not df[['window_id', 'group']].isna().any().any(), 'missing gain identifier')
        numeric(df, ['gain', 'value'])
        require((df.gain > 0).all(), 'gain must be positive')
        require(not df.duplicated(['window_id', 'group', 'gain']).any(), 'duplicate gain cell')
        gains = item['gains']
        require(1 in gains, 'declare unit gain in gains')
        require(set(item['groups']) == set(df.group), 'declare every gain group')
        require(all(set(g.gain) == set(gains) for _, g in df.groupby(['window_id', 'group'])), 'incomplete gain trajectories')
        require(item.get('mode') in ['relative_percent', 'absolute'], 'declare gain mode')
        ref = df[df.gain.eq(1)][['window_id', 'group', 'value']].rename(columns={'value': '_reference'})
        df = df.merge(ref, on=['window_id', 'group'], validate='many_to_one')
        if item['mode'] == 'relative_percent':
            require((df._reference.abs() > item.get('baseline_tolerance', 1e-12)).all(), 'near-zero reference: choose absolute change and label units')
            df['change'] = 100 * (df.value / df._reference - 1)
        else:
            require(item.get('value_units'), 'absolute gain change requires value_units')
            df['change'] = df.value - df._reference
    elif kind == 'agreement':
        require({'reference', 'candidate'} <= set(df), 'agreement needs reference,candidate')
        numeric(df, ['reference', 'candidate'])
        require(item.get('xlabel') and item.get('ylabel'), 'agreement labels required')
    else:
        require(item.get('columns'), 'declare table columns')
        require(set(item['columns']) <= set(df), 'missing requested table columns')
    for col in ['p_value', 'q_value']:
        if col in df:
            numeric(df, [col])
            require(df[col].between(0, 1).all(), f'invalid {col}')
    audit = {'source': str(src), 'sha256': digest(src), 'source_rows': source_rows,
             'retained_rows': len(df), 'configuration': item}
    return df, audit


def encoding(item, group, i):
    spec = item.get('series_styles', {}).get(group, {})
    default_color = '#6A6A6A' if group.startswith('Original') else '#0072B2' if group.startswith('GI') else COLORS[i % len(COLORS)]
    default_marker = 'o' if group.startswith('Original') else 's' if group.startswith('GI') else MARKERS[i % len(MARKERS)]
    return spec.get('color', default_color), spec.get('marker', default_marker), spec.get('linestyle', '-')


def latex_table(df, path):
    """Escaped tabular without an optional pandas/Jinja dependency."""
    escapes = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
               '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}',
               '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    def cell(value):
        if pd.isna(value):
            return '--'
        text = f'{value:.5g}' if isinstance(value, (float, np.floating)) else str(value)
        return ''.join(escapes.get(c, c) for c in text)
    lines = [r'\begin{tabular}{' + 'l' * len(df.columns) + '}', r'\toprule',
             ' & '.join(cell(c) for c in df.columns) + r' \\', r'\midrule']
    lines.extend(' & '.join(cell(v) for v in row) + r' \\' for row in df.itertuples(index=False, name=None))
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    path.write_text('\n'.join(lines) + '\n')


def forest(df, item):
    labels, groups = item['order'], item['groups']
    fig, ax = plt.subplots(figsize=(item.get('width', 7.1), max(2.8, 1.0 + 0.27 * len(labels))), layout='constrained')
    for i, group in enumerate(groups):
        sub = df[df.group.eq(group)].set_index('label').loc[labels]
        y = np.arange(len(labels)) + (i-(len(groups)-1)/2)*0.65/len(groups)
        color, marker, _ = encoding(item, group, i)
        ax.errorbar(sub.estimate, y, xerr=[sub.estimate-sub.ci_low, sub.ci_high-sub.estimate],
                    fmt=marker, color=color, ecolor=color, markersize=3.6, capsize=2,
                    linewidth=0.8, label=group)
    ax.axvline(0, color='#555555', linestyle='--', linewidth=0.8, zorder=0)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlabel(item['xlabel'])
    ax.grid(axis='x', color='#E6E6E6', linewidth=0.5)
    if len(groups) > 1:
        ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=min(3, len(groups)))
    return fig, None


def gain(df, item):
    fig, ax = plt.subplots(figsize=(item.get('width', 7.1), 3.6), layout='constrained')
    summaries = []
    for i, group in enumerate(item['groups']):
        sub = df[df.group.eq(group)]
        summary = sub.groupby('gain').change.agg(median='median', q25=lambda s: s.quantile(.25), q75=lambda s: s.quantile(.75), n='size').reset_index()
        summary['group'] = group
        summaries.append(summary)
        color, marker, line = encoding(item, group, i)
        ax.plot(summary.gain, summary['median'], marker=marker, linestyle=line, color=color, markersize=4, linewidth=1, label=group)
        ax.fill_between(summary.gain, summary.q25, summary.q75, color=color, alpha=.13)
    ax.set_xscale('log', base=2)
    ax.set_xticks(item['gains'])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.axhline(0, color='#666666', linestyle='--', linewidth=.7)
    if item.get('symlog_linthresh'):
        ax.set_yscale('symlog', linthresh=item['symlog_linthresh'])
    ax.set_xlabel('Positive constant gain (×)')
    ax.set_ylabel('Change from unit gain (%)' if item['mode']=='relative_percent' else f"Change from unit gain ({item['value_units']})")
    ax.legend(frameon=False)
    ax.grid(axis='y', color='#E6E6E6', linewidth=.5)
    return fig, pd.concat(summaries, ignore_index=True)


def agreement(df, item):
    fig, ax = plt.subplots(figsize=(3.5, 3.5), layout='constrained')
    ax.scatter(df.reference, df.candidate, s=16, color=COLORS[0], marker='o', alpha=.75)
    lo = min(df.reference.min(), df.candidate.min())
    hi = max(df.reference.max(), df.candidate.max())
    pad = .05 * (hi-lo) if hi>lo else max(abs(lo)*.05, 1e-6)
    ax.plot([lo-pad, hi+pad], [lo-pad, hi+pad], color='#666666', linestyle='--', linewidth=.8)
    ax.set(xlim=(lo-pad, hi+pad), ylim=(lo-pad, hi+pad), xlabel=item['xlabel'], ylabel=item['ylabel'])
    ax.set_aspect('equal', adjustable='box')
    ax.text(.03, .97, f"n = {len(df)}\nMaximum absolute difference = {np.abs(df.reference-df.candidate).max():.3g}", transform=ax.transAxes, va='top', fontsize=8)
    return fig, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--config', type=Path, default=root / 'configs' / 'figures.json')
    parser.add_argument('--out', type=Path, default=root / 'results' / 'manuscript' / 'figures')
    parser.add_argument('--overwrite', action='store_true', help='replace an existing output directory')
    args = parser.parse_args()
    try:
        if args.overwrite and args.out.exists():
            import shutil
            shutil.rmtree(args.out)
        require(not args.out.exists(), 'output directory exists: pass --overwrite or use a new --out')
        config = json.loads(args.config.read_text())
        items = config['items']
        ids = [i['id'] for i in items]
        require(len(ids) == len(set(ids)), 'duplicate display ID')
        require(all(i and all(c.isalnum() or c in '_-' for c in i) for i in ids), 'unsafe display ID')
        ready = []
        for item in items:
            if item.get('enabled', False):
                try:
                    df, audit = prepare(item, args.config.resolve().parent)
                    ready.append((item, df, audit))
                except Exception as exc:
                    raise ValueError(f"{item['id']}: {exc}") from exc
        require(ready, 'no enabled displays: complete Pass 1 configuration first')
        style()
        args.out.mkdir(parents=True)
        manifest = {'status': 'rendering', 'config_sha256': digest(args.config),
                    'script_sha256': digest(__file__), 'python': platform.python_version(),
                    'packages': {'pandas': pd.__version__, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__},
                    'disabled': [i['id'] for i in items if not i.get('enabled', False)], 'outputs': []}
        (args.out/'BUILD_STATUS.json').write_text(json.dumps(manifest, indent=2))
        for item, df, audit in ready:
            stem = args.out/item['id']
            df.to_csv(stem.with_suffix('.csv'), index=False)
            if item['kind']=='table':
                latex_table(df[item['columns']], stem.with_suffix('.tex'))
            else:
                fig, summary = {'forest':forest, 'gain':gain, 'agreement':agreement}[item['kind']](df, item)
                fig.savefig(stem.with_suffix('.pdf'))
                fig.savefig(stem.with_suffix('.png'), dpi=600)
                plt.close(fig)
                if summary is not None:
                    summary.to_csv(args.out/(item['id']+'_descriptive_summary.csv'), index=False)
            stem.with_suffix('.caption.md').write_text(item['caption']+'\n')
            audit['id'] = item['id']
            manifest['outputs'].append(audit)
        manifest['status'] = 'rendered_pending_visual_review'
        manifest['output_hashes'] = {p.name:digest(p) for p in sorted(args.out.iterdir()) if p.name!='BUILD_STATUS.json'}
        (args.out/'BUILD_STATUS.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(f'Rendered {len(ready)} displays. Visual review required: {args.out}')
    except Exception as exc:
        print(f'BUILD FAILED: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
