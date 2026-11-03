#!/usr/bin/env python3
"""Aggregate partial high-divalent MD results with paired Student-t contrasts."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'ComponentEntropyStudy_20261030'
STATES = {'B3': 'baseline_3C', 'B5': 'mono_f010', 'U4': 'weak_f030_e025',
          'L5': 'dispersion_split_d000', 'M5': 'dispersion_split_d025',
          'H5': 'dispersion_split_d040', 'U6': 'mono_weak_uniform6C',
          'H7': 'mono_weak_split7C'}
CONTRASTS = ['H5/U4', 'M5/U4', 'L5/U4', 'H7/U6', 'H5/B3',
             'H7/B5', 'B5/B3', 'U4/B3']
METRICS = ['D_net', 'D_probe', 'D_+2', 'D_-2', 'eta', 'sigma', 's_ex', 'Z', 'u',
           'sigma_2', 'sigma_5', 'sigma_10', 'sigma_20', 'sigma_NE', 'inv_Haven']
REFERENCE_CONTRASTS = {'H5/U4': 'high_split', 'H7/U6': 'mono_high_split',
                       'H5/B3': 'high_vs_baseline', 'U4/B3': None}


def estimate(values):
    values = np.asarray(values, float)
    if not len(values):
        return {'n': 0, 'mean': None, 'CI': None}
    mean = float(values.mean())
    if len(values) == 1:
        return {'n': 1, 'mean': mean, 'CI': None}
    half = float(student_t.ppf(.975, len(values)-1) * values.std(ddof=1) / np.sqrt(len(values)))
    return {'n': len(values), 'mean': mean, 'CI': [mean-half, mean+half]}


def ratio_contrast(left, right):
    samples = [np.log(a/b) for a, b in zip(left, right) if a > 0 and b > 0]
    stat = estimate(samples)
    stat['percent'] = 100 * np.expm1(stat['mean']) if stat['mean'] is not None else None
    stat['percent_CI'] = ([float(100*np.expm1(x)) for x in stat['CI']]
                          if stat['CI'] is not None else None)
    return stat


def metric(row, key):
    transport = row['transport']
    if key == 's_ex':
        return row['s_ex']
    if key in ('Z', 'u'):
        return row[key]
    if key == 'D_+2':
        return transport['D_species']['probe_plus2']
    if key == 'D_-2':
        return transport['D_species']['probe_minus2']
    if key == 'eta':
        return transport['eta']['plateau']
    if key == 'sigma':
        return transport['sigma']['plateau']
    if key.startswith('sigma_') and key[6:] in ('2', '5', '10', '20'):
        return transport['sigma']['cutoffs'][key[6:]]
    if key == 'sigma_NE':
        return transport['sigma']['sigma_NE']
    if key == 'inv_Haven':
        return transport['sigma']['inverse_Haven_ratio']
    return transport[key]


def pooled_inverse_haven(by):
    """Pool independent state replicas, omitting the L5 label-only control."""
    return estimate([metric(row, 'inv_Haven') for code, reps in by.items()
                     if code != 'L5' for row in reps.values()])


FAMILIES = {'no_salt': ['B3', 'U4', 'M5', 'H5'], 'mono_salt': ['B5', 'U6', 'H7']}
FIT_METRICS = ['D_net', 'D_probe', 'D_+2', 'D_-2', 'sigma', 'sigma_NE']
LOCAL_SLOPES = ['H5/U4', 'H7/U6', 'M5/U4']


def log_linear_fit(x, y):
    """Unweighted least squares of ln(y) against x; returns slope B, intercept, R^2."""
    x, logy = np.asarray(x, float), np.log(np.asarray(y, float))
    slope, intercept = np.polyfit(x, logy, 1)
    residual = logy - (intercept + slope * x)
    total = np.sum((logy - logy.mean()) ** 2)
    r2 = 1 - np.sum(residual ** 2) / total if total > 0 else 1.0
    return {'B': float(slope), 'intercept': float(intercept), 'R2': float(r2)}


def rosenfeld_fits(systems):
    """Family fits of ln(mean metric) vs mean s_ex; the Rosenfeld factor is a common constant."""
    return {family: {key: log_linear_fit([systems[c]['metrics']['s_ex']['mean'] for c in codes],
                                         [systems[c]['metrics'][key]['mean'] for c in codes])
                     for key in FIT_METRICS}
            for family, codes in FAMILIES.items()}


def local_slopes(contrasts):
    """Controlled slopes ln(1+percent/100)/Delta s_ex for matched pairs."""
    return {label: {key: float(np.log1p(contrasts[label]['metrics'][key]['percent'] / 100)
                               / contrasts[label]['metrics']['s_ex']['mean'])
                    for key in ('D_net', 'D_probe', 'sigma_NE')}
            for label in LOCAL_SLOPES}


def sign_differences(by):
    return {code: estimate([row['transport']['D_species']['probe_plus2']
                            - row['transport']['D_species']['probe_minus2']
                            for row in by.get(code, {}).values()])
            for code in STATES}


def concentration_fluctuations(by):
    states = {}
    for code in STATES:
        values = [v for row in by.get(code, {}).values()
                  for v in row['transport']['concentration_fluctuations'].values()]
        states[code] = {'min': float(min(values)), 'max': float(max(values))} if values else None
    pooled = [v for v in states.values() if v]
    return {'states': states, 'overall': {'min': min(v['min'] for v in pooled),
                                          'max': max(v['max'] for v in pooled)}}


def format_stat(stat, field='mean', ci='CI', digits=5, scale=1):
    if stat is None or stat.get(field) is None:
        return '—'
    value = f"{scale * stat[field]:.{digits}g}"
    interval = stat.get(ci)
    if interval is not None:
        value += f" [{scale * interval[0]:.{digits}g}, {scale * interval[1]:.{digits}g}]"
    return value


def comparison_data(systems, contrasts, old, follow, hnc):
    details = {}
    reference = follow['paired_contrasts']['long_N4000']
    follow_systems = follow['systems']['long_N4000']
    for label, old_label in REFERENCE_CONTRASTS.items():
        a, b = label.split('/')
        previous = reference[label]['20_100']
        old_entropy = old['contrasts'][old_label]['delta_s_ex'] if old_label else None
        hnc_pair = hnc.get('contrasts', {}).get(label.replace('/', '-'))
        details[label] = {
            'x075': {key: contrasts[label]['metrics'][key] for key in
                     ('D_net', 'D_probe', 'eta', 's_ex')},
            'x003': {
                'D_net': {'percent': 100 * previous['net']['ratio_minus_one'],
                          'percent_CI': [100 * x for x in previous['net']['ratio_CI']]},
                'D_probe': {'percent': 100 * previous['probe']['ratio_minus_one'],
                            'percent_CI': [100 * x for x in previous['probe']['ratio_CI']]},
                'eta': {'percent': 100 * (follow_systems[a]['eta']['mean'] /
                                          follow_systems[b]['eta']['mean'] - 1),
                        'method': 'ratio of means'},
                's_ex': old_entropy},
            'hnc_x075_delta_s_ex': hnc_pair['new'] if hnc_pair else None,
        }
    shifts = {}
    for code in ('B3', 'U4', 'H5', 'U6', 'H7'):
        current = systems[code]['metrics']
        prior = follow_systems[code]
        old_system = old['systems'].get(STATES[code], {})
        shifts[code] = {}
        for key in ('s_ex', 'D_probe', 'eta', 'Z', 'u'):
            new = current[key]
            earlier = old_system.get('s_ex') if key == 's_ex' else prior.get(key)
            delta = (new['mean'] - earlier['mean'] if earlier is not None else None)
            shifts[code][key] = {'x075': new, 'x003': earlier,
                                 'difference': (delta if key not in ('D_probe', 'eta')
                                                else 100 * (new['mean'] / earlier['mean'] - 1))
                                 if earlier is not None else None,
                                 'difference_unit': '%' if key in ('D_probe', 'eta') else 'absolute'}
    return {'contrasts': details, 'absolute_shifts': shifts}


def plot_contrasts(contrasts, comparison, out):
    labels = ['H5/U4', 'M5/U4', 'H7/U6', 'H5/B3', 'H7/B5', 'U4/B3', 'B5/B3']
    keys = ('D_probe', 'eta', 'sigma', 'sigma_NE')
    colors = ('#326a99', '#dc8d34', '#4b936c', '#9a5aa1')
    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(labels))
    for j, (key, color) in enumerate(zip(keys, colors)):
        centers = x + (j - 1.5) * .19
        stats = [contrasts[label]['metrics'][key] for label in labels]
        values = np.array([s['percent'] for s in stats])
        errors = np.array([[s['percent'] - s['percent_CI'][0] for s in stats],
                           [s['percent_CI'][1] - s['percent'] for s in stats]])
        ax.bar(centers, values, width=.18, color=color, label=key)
        ax.errorbar(centers, values, yerr=errors, fmt='none', ecolor='#292929',
                    capsize=2, linewidth=.8)
    for i, label in enumerate(labels):
        ref = comparison['contrasts'].get(label)
        if ref:
            stat = ref['x003']['D_probe']
            ax.errorbar(x[i] - .285, stat['percent'],
                        yerr=[[stat['percent'] - stat['percent_CI'][0]],
                              [stat['percent_CI'][1] - stat['percent']]],
                        fmt='D', color='#182b41', markersize=5, capsize=3,
                        label='x=0.03 D_probe' if i == 0 else None)
    ax.axhline(0, color='#444444', linewidth=.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel('Paired contrast (%)')
    ax.set_title('x=0.075 MD transport contrasts (95% Student-t CI)')
    ax.legend(ncol=5, fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.1))
    ax.grid(axis='y', alpha=.2)
    fig.tight_layout()
    fig.savefig(out / 'md_contrasts.pdf', bbox_inches='tight')
    fig.savefig(out / 'md_contrasts.png', dpi=180, bbox_inches='tight')
    plt.close(fig)


def aggregate():
    by = {}
    for path in sorted((ROOT / 'md').glob('*/replica_*/replica_result.json')):
        row = json.loads(path.read_text())
        by.setdefault(row['code'], {})[row['replica']] = row
    systems = {}
    for code in STATES:
        reps = by.get(code, {})
        systems[code] = {'missing_replicas': sorted(set(range(1, 5)) - set(reps)),
                         'metrics': {key: estimate([metric(row, key) for row in reps.values()])
                                     for key in METRICS}}
    contrasts = {}
    for label in CONTRASTS:
        a, b = label.split('/')
        common = sorted(set(by.get(a, {})) & set(by.get(b, {})))
        contrasts[label] = {'paired_replicas': common, 'missing_pairs': sorted(set(range(1, 5))-set(common)),
                            'metrics': {key: (estimate([metric(by[a][i], key)-metric(by[b][i], key)
                                                        for i in common]) if key == 's_ex' else
                                              ratio_contrast([metric(by[a][i], key) for i in common],
                                                             [metric(by[b][i], key) for i in common]))
                                        for key in METRICS if key not in ('Z', 'u')}}
    old = json.loads((OLD / 'results/summary.json').read_text())
    follow = json.loads((OLD / 'followup_20261101/results/summary.json').read_text())
    hnc_path = ROOT / 'results/hnc_summary.json'
    hnc = json.loads(hnc_path.read_text()) if hnc_path.exists() else {}
    cross_comparison = comparison_data(systems, contrasts, old, follow, hnc)
    comparison = {}
    for code, name in STATES.items():
        comparison[code] = {
            'x015_N4000_4x3000tau': systems[code]['metrics'],
            'x003_N4000_3x1000tau_EOS': old['systems'].get(name),
            'x003_N4000_6x5000tau_followup': follow['systems']['long_N4000'].get(code),
        }
    summary = {'completed_replicas': sum(map(len, by.values())), 'expected_replicas': 32,
               'statistics': '95% Student-t CIs across replicas; paired percent = 100*(exp(mean ln ratio)-1)',
               'systems': systems, 'paired_contrasts': contrasts,
               'x003_comparison': comparison,
               'x003_contrast_comparison': cross_comparison,
               'pooled_inverse_Haven': pooled_inverse_haven(by),
               'sign_difference': sign_differences(by),
               'concentration_fluctuations': concentration_fluctuations(by),
               'rosenfeld_fits': rosenfeld_fits(systems),
               'local_slopes': local_slopes(contrasts),
               'source_files': {'EOS': str(OLD / 'results/summary.json'),
                                'followup': str(OLD / 'followup_20261101/results/summary.json')}}
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    lines = [f'# High-divalent MD: {summary["completed_replicas"]}/32 replicas', '',
             '95% Student-t intervals use independent replicas. Paired transport contrasts are',
             '100[exp(mean ln ratio)-1]%; entropy contrasts are direct differences.', '',
             '| State | n | Missing replicas | s_ex | D_net | D_probe | eta | sigma |',
             '|---|---:|---|---:|---:|---:|---:|---:|']
    for code, row in systems.items():
        metrics = row['metrics']
        fmt = lambda key, data=metrics: format_stat(data[key])
        lines.append(f"| {code} | {metrics['s_ex']['n']} | {row['missing_replicas']} | "
                     f"{fmt('s_ex')} | {fmt('D_net')} | {fmt('D_probe')} | {fmt('eta')} | {fmt('sigma')} |")
    lines += ['', '## Paired contrasts', '', '| Contrast | Pairs | D_net % | D_probe % | eta % | sigma % | Δs_ex |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for label, row in contrasts.items():
        fmt = lambda key, data=row['metrics']: format_stat(data[key], 'percent', 'percent_CI', 4)
        lines.append(f"| {label} | {len(row['paired_replicas'])} | {fmt('D_net')} | {fmt('D_probe')} | "
                     f"{fmt('eta')} | {fmt('sigma')} | {format_stat(row['metrics']['s_ex'])} |")
    lines += ['', '## Conductivity convergence and Nernst-Einstein', '',
              'Conductivity cutoffs are running Green-Kubo integrals at lag τ. '
              'The reported sigma is the 5–10 τ plateau mean.', '',
              '| State | sigma_2 | sigma_5 | sigma_10 | sigma_20 | plateau sigma | sigma_NE | inverse Haven |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
    for code, row in systems.items():
        stats = row['metrics']
        lines.append('| ' + code + ' | ' + ' | '.join(format_stat(stats[key]) for key in
                     ('sigma_2', 'sigma_5', 'sigma_10', 'sigma_20', 'sigma', 'sigma_NE', 'inv_Haven')) + ' |')
    lines += ['', f"Pooled inverse Haven, excluding L5 (n={summary['pooled_inverse_Haven']['n']}): "
              f"{format_stat(summary['pooled_inverse_Haven'])}.", '',
              '### Paired conductivity contrasts', '',
              '| Contrast | sigma % | sigma_NE % | inverse Haven % |',
              '|---|---:|---:|---:|']
    for label, row in contrasts.items():
        stats = row['metrics']
        lines.append('| ' + label + ' | ' + ' | '.join(
            format_stat(stats[key], 'percent', 'percent_CI', 4)
            for key in ('sigma', 'sigma_NE', 'inv_Haven')) + ' |')
    lines += ['', '## Comparison with x=0.03', '',
              'Percentages are paired log-ratio contrasts except x=0.03 eta, which is a ratio of means.',
              'x=0.03 transport uses the 20–100 τ window; x=0.03 entropy uses the EOS campaign.', '',
              '| Contrast | x=0.075 D_net % | x=0.03 D_net % | x=0.075 D_probe % | x=0.03 D_probe % | x=0.075 eta % | x=0.03 eta % (ratio of means) | x=0.075 Δs_ex | x=0.03 Δs_ex | HNC x=0.075 Δs_ex |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for label, row in cross_comparison['contrasts'].items():
        new, prior = row['x075'], row['x003']
        cells = [format_stat(new[key], 'percent', 'percent_CI', 4) if key != 's_ex'
                 else format_stat(new[key]) for key in ('D_net', 'D_probe', 'eta', 's_ex')]
        old_cells = [format_stat(prior[key], 'percent', 'percent_CI', 4) if key != 's_ex'
                     else format_stat(prior[key]) for key in ('D_net', 'D_probe', 'eta', 's_ex')]
        lines.append(f"| {label} | {cells[0]} | {old_cells[0]} | {cells[1]} | {old_cells[1]} | "
                     f"{cells[2]} | {old_cells[2]} | {cells[3]} | {old_cells[3]} | "
                     f"{row['hnc_x075_delta_s_ex'] if row['hnc_x075_delta_s_ex'] is not None else '—'} |")
    lines += ['', '### Absolute shift: x=0.075 versus x=0.03', '',
              'Differences are x=0.075 minus x=0.03; D_probe and eta differences are percentages.', '',
              '| State | Metric | x=0.075 | x=0.03 | Difference |',
              '|---|---|---:|---:|---:|']
    for code, metrics in cross_comparison['absolute_shifts'].items():
        for key, entry in metrics.items():
            diff = entry['difference']
            suffix = '%' if entry['difference_unit'] == '%' else ''
            difference = f'{diff:.5g}{suffix}' if diff is not None else '—'
            lines.append(f"| {code} | {key} | {format_stat(entry['x075'])} | "
                         f"{format_stat(entry['x003'])} | {difference} |")
    lines += ['', '## Charge-sign symmetry and homogeneity', '',
              '| State | D_+2 - D_-2 | S_cc range |', '|---|---:|---:|']
    for code in STATES:
        cf = summary['concentration_fluctuations']['states'][code]
        lines.append(f"| {code} | {format_stat(summary['sign_difference'][code])} | "
                     f"{cf['min']:.3f}–{cf['max']:.3f} |")
    overall = summary['concentration_fluctuations']['overall']
    lines += ['', f"Low-k concentration fluctuations overall: {overall['min']:.3f}–{overall['max']:.3f}.",
              '', '## Rosenfeld-like family fits', '',
              'Unweighted fits of ln(mean metric) against mean s_ex (L5 excluded).', '',
              '| Family | Metric | B | R² |', '|---|---|---:|---:|']
    for family, fits in summary['rosenfeld_fits'].items():
        for key, fit in fits.items():
            lines.append(f"| {family} | {key} | {fit['B']:.3f} | {fit['R2']:.4f} |")
    lines += ['', '| Contrast | local slope D_net | D_probe | sigma_NE |', '|---|---:|---:|---:|']
    for label, slopes in summary['local_slopes'].items():
        lines.append(f"| {label} | {slopes['D_net']:.3f} | {slopes['D_probe']:.3f} | {slopes['sigma_NE']:.3f} |")
    lines += ['', '## x=0.03 reference', '',
              'EOS: ComponentEntropyStudy_20261030/results/summary.json (3 replicas, 1000 tau dense).',
              'Transport: ComponentEntropyStudy_20261030/followup_20261101/results/summary.json '
              '(6 replicas, 5000 tau; B3/U4/H5/U6/H7 only).',
              'The x=0.03 and x=0.075 rows are independent campaigns; see summary.json for all values and CIs.', '']
    (out / 'report.md').write_text('\n'.join(lines))
    plot_contrasts(contrasts, cross_comparison, out)
    print(f"Aggregated {summary['completed_replicas']}/32 replicas")


if __name__ == '__main__':
    aggregate()
