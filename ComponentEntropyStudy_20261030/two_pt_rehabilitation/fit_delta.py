#!/usr/bin/env python3
"""Post-hoc global calibration of the revised-2PT partition exponent."""

import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import t

from run_rehabilitation import CACHE, HERE, STUDY, cache_replica, digest, spectrum
from two_pt import entropy_from_moments, evaluate, log_moment, trap

ORDER = [('B3', 'baseline_3C'), ('B5', 'mono_f010'),
         ('U4', 'weak_f030_e025'), ('L5', 'dispersion_split_d000'),
         ('M5', 'dispersion_split_d025'), ('H5', 'dispersion_split_d040'),
         ('U6', 'mono_weak_uniform6C'), ('H7', 'mono_weak_split7C')]
FIT_CODES = [code for code, _ in ORDER if code != 'L5']
SALT = {'B5', 'U6', 'H7'}
PAIRS = [('M5/U4', 'M5', 'U4'), ('H5/U4', 'H5', 'U4'),
         ('H7/U6', 'H7', 'U6'), ('H5/B3', 'H5', 'B3'),
         ('H7/B5', 'H7', 'B5'), ('B5/B3', 'B5', 'B3'),
         ('L5/U4', 'L5', 'U4')]
GRID = np.arange(160, 401, dtype=float) / 200
REPS = (1, 2, 3)


def estimate(values):
    values = np.asarray(values, float)
    mean = float(values.mean())
    half = float(t.ppf(.975, len(values) - 1) * values.std(ddof=1)
                 / np.sqrt(len(values)))
    return {'values': values.tolist(), 'mean': mean,
            'CI': [mean - half, mean + half]}


def compact(values):
    e = estimate(values)
    return {k: e[k] for k in ('mean', 'CI')}


def fmt(e, digits=5):
    return (f"{e['mean']:+.{digits}f} "
            f"[{e['CI'][0]:+.{digits}f}, {e['CI'][1]:+.{digits}f}]")


def main():
    data = {}
    for code, name in ORDER:
        data[code] = {}
        for rep in REPS:
            folder = STUDY / 'md' / name / f'replica_{rep}'
            meta, arrays = cache_replica(folder)
            assert meta['name'] == name and meta['replica'] == rep
            stored = json.loads((HERE / 'replicas' / f'{name}_r{rep}.json').read_text())
            ti = json.loads((folder / 'replica_result.json').read_text())['s_ex']
            assert ti == stored['TI_s_ex']
            groups = []
            for i, group in enumerate(meta['groups']):
                if not group['primary']:
                    continue
                nu, G, _ = spectrum(arrays[f'g{i}_full'], meta['dt'], 10, 2.)
                moments = (trap(G, nu), log_moment(nu, G), float(G[0]))
                groups.append((group, moments, nu, G))
            assert abs(sum(g['x'] for g, *_ in groups) - 1) < 1e-12

            def evaluate_rep(delta, details=False, groups=groups):
                mixture = 0.
                species = {}
                for group, (modes, logmoment, s0), _, _ in groups:
                    result = entropy_from_moments(modes, logmoment, s0, 2., .7,
                                                  .7, exponent=float(delta))
                    mixture += group['x'] * result['s_ex']
                    if details:
                        species[group['label']] = {
                            's_ex': result['s_ex'], 'f': result['f'],
                            'x': group['x'], 'probe': group['probe'],
                            'primary_solvent': group['key'][0] == 0
                            and group['key'][1] == 1.}
                return (mixture, species) if details else mixture

            for delta in (1., 1.5):
                direct = sum(g['x'] * evaluate(nu, G, exponent=delta,
                                              diagnostics=False)['s_ex']
                             for g, _, nu, G in groups)
                saved = stored['variants']['full'][f'tau10_d{delta:g}']['s_ex']
                assert abs(direct - saved) <= 1e-10, (code, rep, delta, direct, saved)
                assert abs(evaluate_rep(delta) - saved) <= 1e-10
            curves = [evaluate_rep(delta, details=True) for delta in GRID]
            data[code][rep] = {'TI': ti, 'curve': np.array([v[0] for v in curves]),
                               'species': {label: {'s_ex': [v[1][label]['s_ex'] for v in curves],
                                                   'f': [v[1][label]['f'] for v in curves],
                                                   'x': curves[0][1][label]['x'],
                                                   'probe': curves[0][1][label]['probe'],
                                                   'primary_solvent': curves[0][1][label]['primary_solvent']}
                                           for label in curves[0][1]},
                               'eval': evaluate_rep,
                               'cache': str((CACHE / f'{name}_r{rep}.npz').relative_to(HERE))}
        print(f'{code}: checked three stored replicas at delta=1 and 1.5', flush=True)

    def value(code, rep, delta):
        return data[code][rep]['eval'](delta)

    def bias(code, delta, reps=REPS):
        return np.mean([value(code, rep, delta) - data[code][rep]['TI']
                        for rep in reps])

    def fit(codes, reps=REPS, weights=None):
        def loss(delta):
            b = np.array([bias(code, delta, reps) for code in codes])
            return float(np.sum(b * b if weights is None else weights * b * b))
        result = minimize_scalar(loss, bounds=(.8, 2.), method='bounded',
                                 options={'xatol': 1e-12})
        candidates = [(result.fun, result.x), (loss(.8), .8), (loss(2.), 2.)]
        return float(min(candidates)[1])

    optimum = fit(FIT_CODES)
    replica_fits = estimate([fit(FIT_CODES, (rep,)) for rep in REPS])
    variances = np.array([np.var([value(code, rep, optimum)
                                   - data[code][rep]['TI'] for rep in REPS], ddof=1) / 3
                          for code in FIT_CODES])
    assert np.all(variances > 0)
    weighted = fit(FIT_CODES, weights=1 / variances)
    methods = {'1.0': 1., '1.25': 1.25, '1.5': 1.5, 'fitted': optimum}
    state_results = {}
    for code, _ in ORDER:
        bcurve = np.mean([data[code][r]['curve'] - data[code][r]['TI']
                          for r in REPS], axis=0)
        def root(reps, code=code):
            a, b = bias(code, .8, reps), bias(code, 2., reps)
            return float(brentq(lambda delta, code=code, reps=reps: bias(code, delta, reps),
                                .8, 2.)) if a * b <= 0 else None
        matches = [root((rep,)) for rep in REPS]
        state_results[code] = {
            'name': dict(ORDER)[code],
            'TI': estimate([data[code][r]['TI'] for r in REPS]),
            'delta_state': root(REPS),
            'replica_matching_exponents': matches,
            'replica_matching_spread': compact(matches) if all(v is not None for v in matches) else None,
            'biases': {key: estimate([value(code, r, d) - data[code][r]['TI']
                                      for r in REPS]) for key, d in methods.items()},
            'two_pt': {key: estimate([value(code, r, d) for r in REPS])
                       for key, d in methods.items()},
            'curve_mean_bias': bcurve.tolist(),
            'curve_replica_bias': [(data[code][r]['curve'] - data[code][r]['TI']).tolist()
                                   for r in REPS],
            'species': {str(r): data[code][r]['species'] for r in REPS},
        }
    metrics = {}
    for key in methods:
        b = np.array([state_results[c]['biases'][key]['mean'] for c in FIT_CODES])
        metrics[key] = {'RMS_absolute_state_bias': float(np.sqrt(np.mean(b**2))),
                        'max_absolute_state_bias': float(np.max(abs(b)))}
    increments = {}
    for label, a, b in PAIRS:
        ti = [data[a][r]['TI'] - data[b][r]['TI'] for r in REPS]
        increments[label] = {
            'null_check': label == 'L5/U4', 'TI': estimate(ti),
            'methods': {key: {'delta': estimate([value(a, r, d) - value(b, r, d)
                                                  for r in REPS]),
                              'bias': estimate([value(a, r, d) - value(b, r, d) - ti[r-1]
                                                for r in REPS])}
                        for key, d in methods.items()}}
    loso = {}
    for code in FIT_CODES:
        delta = fit([c for c in FIT_CODES if c != code])
        loso[code] = {'training_delta': delta, 'held_out_bias': float(bias(code, delta))}
    loso_rms = float(np.sqrt(np.mean([v['held_out_bias']**2 for v in loso.values()])))
    no_salt = [state_results[c]['delta_state'] for c in FIT_CODES if c not in SALT]
    salt = [state_results[c]['delta_state'] for c in FIT_CODES if c in SALT]
    output = {
        'provenance': {'note': 'Post-hoc calibration to TI; not prespecified in protocol.json',
                       'settings': {'segment': 'full', 'taper_tau': 10, 'temperature': 2.,
                                    'density': .7, 'grid': {'start': .8, 'stop': 2., 'step': .005},
                                    'independent_states': FIT_CODES, 'excluded_duplicate': 'L5'},
                       'code_sha256': {n: digest(HERE / n) for n in
                                       ('fit_delta.py', 'run_rehabilitation.py', 'two_pt.py')},
                       'protocol_sha256': digest(HERE / 'protocol.json'),
                       'cache_files': [data[c][r]['cache'] for c, _ in ORDER for r in REPS],
                       'stored_reproduction_tolerance': 1e-10,
                       'weighted_fit': 'Inverse variance of three-replica mean bias, fixed at unweighted optimum'},
        'grid': GRID.tolist(), 'global_fit': {'delta_star': optimum,
            'replica_index_fits': replica_fits, 'inverse_variance_weighted_delta': weighted,
            'inverse_variance_weights': dict(zip(FIT_CODES, (1 / variances).tolist()))},
        'methods': methods, 'states': state_results, 'metrics': metrics,
        'increments': increments, 'loso': {'states': loso, 'RMS_held_out_bias': loso_rms},
        'matching_exponent_groups': {'without_monovalent_salt': no_salt,
                                     'with_monovalent_salt': salt,
                                     'without_mean': float(np.mean(no_salt)),
                                     'with_mean': float(np.mean(salt))},
    }
    (HERE / 'delta_fit.json').write_text(json.dumps(output, indent=2) + '\n')
    write_report(output)
    write_figures(output)
    print(f'Global delta={optimum:.8f}; LOSO RMS={loso_rms:.6f}; output written.', flush=True)


def write_report(o):
    fit = o['global_fit']; lines = [
        '# Post-hoc global partition-exponent calibration', '',
        'This delta fit uses TI results from the same seven independent states, so it is exploratory and cannot be treated as prespecified validation. L5 is the exact mechanical duplicate of U4 and is a null check only. The published protocol used delta=1 and an unfitted delta=1.5 sensitivity.', '',
        f"Unweighted global delta* = **{fit['delta_star']:.7f}**. Replica-index fits: {fmt(fit['replica_index_fits'], 6)} (nominal 95% Student-t CI, df=2). Inverse-variance weighted fit: {fit['inverse_variance_weighted_delta']:.7f}; weights are fixed using the replica-mean bias variance at the unweighted optimum. Leave-one-state-out RMS held-out bias: {o['loso']['RMS_held_out_bias']:.5f} kB.", '',
        '## Absolute-state bias (2PT − TI), kB per particle', '',
        '| State | Matching delta | Replica-root 95% CI | delta=1 | delta=1.25 | delta=1.5 | delta* |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for code, _ in ORDER:
        s = o['states'][code]
        spread = s['replica_matching_spread']
        ci = (f"[{spread['CI'][0]:.5f}, {spread['CI'][1]:.5f}]"
              if spread is not None else 'unavailable')
        lines.append('| ' + code + ' | ' + (f"{s['delta_state']:.5f}" if s['delta_state'] is not None else 'outside grid') + ' | '
                     + ci + ' | '
                     + ' | '.join(fmt(s['biases'][k]) for k in o['methods']) + ' |')
    lines += ['', '| Exponent | RMS state bias | Maximum absolute state bias |',
              '|---|---:|---:|']
    for key, m in o['metrics'].items():
        lines.append(f"| {key} | {m['RMS_absolute_state_bias']:.5f} | {m['max_absolute_state_bias']:.5f} |")
    lines += ['', '## Paired increments, kB per particle', '',
              '| Comparison | TI | 2PT delta=1.25 | Bias [95% CI] | 2PT delta* | Bias [95% CI] |',
              '|---|---:|---:|---:|---:|---:|']
    for label, p in o['increments'].items():
        a, b = p['methods']['1.25'], p['methods']['fitted']
        lines.append(f"| {label}{' (null)' if p['null_check'] else ''} | {p['TI']['mean']:+.5f} | {a['delta']['mean']:+.5f} | {fmt(a['bias'])} | {b['delta']['mean']:+.5f} | {fmt(b['bias'])} |")
    lines += ['', '## Leave-one-state-out prediction', '',
              '| Held-out state | Training delta | Held-out bias |', '|---|---:|---:|']
    for code, v in o['loso']['states'].items():
        lines.append(f"| {code} | {v['training_delta']:.5f} | {v['held_out_bias']:+.5f} |")
    g = o['matching_exponent_groups']
    lines += ['', f"Mean state-matching delta without monovalent salt (B3, U4, M5, H5): {g['without_mean']:.5f}; with monovalent salt (B5, U6, H7): {g['with_mean']:.5f}. The state roots vary from 1.065 (H7) to 1.321 (B5), while the two salt-group means are nearly equal. This is state variation without a simple salt/no-salt separation; the group means are descriptive, not an independent significance test.", '',
              'Replica bands and all species entropy/fluidicity curves are in [delta_fit.json](delta_fit.json). Figures: [bias curves](delta_fit.pdf) and [absolute entropy comparison](entropy_validation_delta_fit.pdf). Intervals quantify three-replica sampling only; calibration and cross-validation use the same finite set of states.', '']
    (HERE / 'delta_fit.md').write_text('\n'.join(lines))


def write_figures(o):
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    grid = np.array(o['grid'])
    fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
    for code in FIT_CODES:
        s = o['states'][code]; band = np.array(s['curve_replica_bias'])
        color = '#b2182b' if code in SALT else '#2166ac'
        mean = band.mean(axis=0)
        half = t.ppf(.975, 2) * band.std(axis=0, ddof=1) / np.sqrt(3)
        ax.plot(grid, mean, color=color, lw=1.2, label=code)
        ax.fill_between(grid, mean-half, mean+half, color=color, alpha=.045)
    ax.axhline(0, color='black', lw=.8)
    for key, d, style in [('1.0', 1., ':'), ('1.25', 1.25, '--'),
                          ('1.5', 1.5, ':'), ('fitted', o['global_fit']['delta_star'], '-.')]:
        ax.axvline(d, color='#555555', lw=1, ls=style)
        ax.text(d, ax.get_ylim()[1], key, rotation=90, va='top', ha='right', fontsize=8)
    ax.set(xlabel=r'Partition exponent $\delta$', ylabel=r'2PT − TI $s_{ex}/k_B$',
           title='Post-hoc exponent sensitivity; nominal 95% replica bands')
    ax.grid(alpha=.2); ax.legend(ncol=4, frameon=False, fontsize=8)
    for ext in ('pdf', 'png'):
        fig.savefig(HERE / f'delta_fit.{ext}', dpi=210)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout='constrained')
    for ax, key, title in zip(axes, ('1.0', '1.25', 'fitted'),
                              ('delta=1', 'delta=1.25', f"fitted delta={o['global_fit']['delta_star']:.3f}")):
        limits = []
        for code in FIT_CODES:
            s = o['states'][code]; x = s['TI']; y = s['two_pt'][key]
            color = '#b2182b' if code in SALT else '#2166ac'
            ax.errorbar(x['mean'], y['mean'], xerr=(x['CI'][1]-x['CI'][0])/2,
                        yerr=(y['CI'][1]-y['CI'][0])/2, fmt='o', color=color, capsize=3)
            ax.annotate(code, (x['mean'], y['mean']), xytext=(4, 5),
                        textcoords='offset points', fontsize=8)
            limits.extend(x['CI'] + y['CI'])
        lo, hi = min(limits)-.025, max(limits)+.025
        ax.plot([lo, hi], [lo, hi], 'k--', lw=.8)
        ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=r'TI $s_{ex}/k_B$',
               ylabel=r'R2PT $s_{ex}/k_B$', title=title)
        ax.grid(alpha=.2)
    fig.suptitle('Nominal 95% replica intervals; dashed line: equality', fontsize=10)
    for ext in ('pdf', 'png'):
        fig.savefig(HERE / f'entropy_validation_delta_fit.{ext}', dpi=210)
    plt.close(fig)


if __name__ == '__main__':
    main()
