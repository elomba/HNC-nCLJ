#!/usr/bin/env python3
"""Package completed dense MD trajectories for external Voronoi analysis."""
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parent.parent
STUDY = ROOT / 'ComponentEntropyStudy_20261030'
PAIR = ROOT / 'MDPairEntropyAssessment'
OUT = ROOT / 'VoronoiEntropyData'
ORDER = [('B3', 'baseline_3C'), ('B5', 'mono_f010'),
         ('U4', 'weak_f030_e025'), ('L5', 'dispersion_split_d000'),
         ('M5', 'dispersion_split_d025'), ('H5', 'dispersion_split_d040'),
         ('U6', 'mono_weak_uniform6C'), ('H7', 'mono_weak_split7C')]
COMPARISONS = {
    'B3': 'Original baseline; compare B5 and B3/U4/H5.',
    'B5': 'Monovalent baseline; compare B3 and B5/U6/H7.',
    'U4': 'Label null with L5; matched-amplitude series U4/M5/H5.',
    'L5': 'Label-only null with U4; two epsilon=0.25 solvents are identical.',
    'M5': 'Matched-amplitude physical split in U4/M5/H5.',
    'H5': 'Matched-amplitude physical split in U4/M5/H5; compare B3.',
    'U6': 'Uniform solvent with monovalent background; compare H7 and B5.',
    'H7': 'Diverse solvent with monovalent background; compare U6 and B5.',
}
DESCRIPTIONS = {
    'B3': '94% LJ solvent with +2/-2 probes (baseline electrolyte)',
    'B5': 'B3 plus 5% of each monovalent ion sign (salt-addition baseline)',
    'U4': '64% LJ solvent + 30% uniform weak solvent (epsilon=0.25)',
    'L5': 'U4 with the weak solvent split into two identical labels (label-only null)',
    'M5': 'U4 weak solvent split into epsilon=0.0625/0.5625 (moderate diversification)',
    'H5': 'U4 weak solvent split into epsilon=0.01/0.81 (strong diversification)',
    'U6': '54% LJ solvent + 30% uniform weak solvent + 5% of each monovalent sign',
    'H7': 'U6 weak solvent split into epsilon=0.01/0.81 (strong diversification with salt)',
}
LABELS = ['solvent', 'weak_solvent', 'weak_A', 'weak_B',
          'background_plus1', 'background_minus1', 'probe_plus2', 'probe_minus2']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def create(path, content):
    """Create text output, refusing to replace non-identical content."""
    data = content.encode() if isinstance(content, str) else content
    if path.exists():
        if path.read_bytes() != data:
            raise FileExistsError(f'Existing generated file differs: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def copy_checked(src, dest):
    source_hash = digest(src)
    if dest.exists():
        if digest(dest) != source_hash:
            raise FileExistsError(f'Existing copy differs: {dest}')
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        if digest(dest) != source_hash:
            raise OSError(f'Copy verification failed: {dest}')
    return source_hash


def csv_text(rows, columns):
    import io
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def mean_ci(values):
    a = np.asarray(values, dtype=float)
    return float(a.mean()), float(t.ppf(.975, len(a)-1) * a.std(ddof=1) / np.sqrt(len(a)))


def per_state_readme(code, item, rows):
    species = item['species']
    lines = [f'# {code}: {DESCRIPTIONS[code]}', '', COMPARISONS[code], '',
             '## Model and composition', '',
             'T*=2, rho*=0.7, N=4000; periodic cubic box with '
             f'L=(4000/0.7)^(1/3)={((4000/.7)**(1/3)):.9f} sigma. '
             'LJ reduced units: sigma=1, mass=1, primary epsilon=1, '
             'tau=sigma sqrt(m/epsilon_1).', '',
             '| LAMMPS type | Label | q | x | N_i | epsilon_i | sigma | mass |',
             '|---:|---|---:|---:|---:|---:|---:|---:|']
    for j, sp in enumerate(species, 1):
        lines.append(f'| {j} | {sp["label"]} | {sp["q"]} | {sp["x"]} | '
                     f'{round(4000*sp["x"])} | {sp["eps"]:.6g} | {sp["sigma"]} | {sp["mass"]} |')
    lines += ['', f'I = {item["I"]:.2f}; ideal mixing entropy '
              f'-sum(x_i ln x_i) = {item["mixing_entropy"]:.8f} kB/particle '
              '(excluded from excess entropy).', '',
              'Pair style `lj/cut/coul/long`; shifted LJ cutoff 2.5 sigma, no tail correction. '
              'Geometric epsilon_ij=sqrt(epsilon_i epsilon_j); all sigma_ij=1. '
              'Coulomb prefactor 1 gives Bjerrum length 0.5 at T=2. '
              'PPPM accuracy 1e-5 and real-space Coulomb cutoff 5 sigma.', '',
              '| type i | type j | epsilon_ij |', '|---:|---:|---:|']
    for i, a in enumerate(species, 1):
        for j, b in enumerate(species[i-1:], i):
            lines.append(f'| {i} | {j} | {math.sqrt(a["eps"]*b["eps"]):.8g} |')
    lines += ['', 'NVT Nose–Hoover thermostat (Tdamp=0.1 tau), dt=0.002 tau; '
              '100,000 equilibration and 500,000 production steps. '
              '2001 positions frames every 250 steps = 0.5 tau over 1000 tau. '
              'Coordinates are unwrapped; wrap by the per-frame cell lengths before periodic Voronoi analysis.',
              '', '## Replicas and references', '',
              '| replica | seed | TI s_ex | HNC s_ex | MD s2 | D_net | D_probe |',
              '|---:|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        lines.append('| ' + ' | '.join((f'{row[x]:.6g}' if isinstance(row[x], float) else str(row[x])) for x in
                     ['replica', 'seed', 's_ex_TI', 's_ex_HNC', 's2_MD', 'D_net', 'D_probe']) + ' |')
    lines += ['', 'The transport coefficients use a common whole-system COM frame; '
              '20–100 tau is the primary MSD fitting window. See `reference.json` '
              'and the top-level CSV for all species and alternate probe windows.', '']
    if code == 'L5':
        lines += ['The two epsilon=0.25 neutral labels are mechanically identical. '
                  'Merge them for any label-sensitive structural estimator and report '
                  'both merged and unmerged results.', '']
    return '\n'.join(lines)


def top_readme(items, rows, means):
    lines = ['# Voronoi entropy analysis input', '',
             'Local-only, checksum-verified copies of the completed dense MD trajectories. '
             'Use these for independent Voronoi topology and cluster-entropy analysis; '
             'no structural entropy estimates are asserted here. See '
             '`../doc/voronoi_entropy_work_plan.tex` for the scientific plan.', '',
             '## States', '',
             '| code | Description | labels / distinct species | I | HNC s_ex | MD TI s_ex mean |',
             '|---|---|---:|---:|---:|---:|']
    by_code = {row['code']: row for row in means}
    for code, name in ORDER:
        item = items[name]
        lines.append(f'| {code} | {DESCRIPTIONS[code]} | '
                     f'{len(item["species"])} / {item["mechanically_distinct_components"]} | '
                     f'{item["I"]:.2f} | {rows[code][0]["s_ex_HNC"]:.5f} | '
                     f'{by_code[code]["s_ex_TI_mean"]:.5f} |')
    lines += ['', 'B3/B5 are baselines; U4/L5 is the label null; U4/M5/H5 have '
              'matched mean bare LJ amplitude; U6/H7 repeat the uniform/diverse '
              'comparison with ±1 background.', '',
              '## Layout and trajectory', '',
              '`<CODE>_<name>/replica_<k>/` holds copied `positions.nc`, exact '
              '`in.run`, initial `data.mixture` (type and charge by atom id), '
              '`log.lammps`, and derived `reference.json`. State READMEs describe '
              'all species and pair amplitudes. `reference_values.csv` contains 24 '
              'replica rows; `reference_state_means.csv` contains means and nominal '
              '95% Student-t half-widths (3 independent replicas). '
              '`MANIFEST.json` and `checksums.sha256` map each copied file to its '
              'source and SHA256.', '',
              'LAMMPS `dump netcdf` writes AMBER-convention NetCDF: 2001 frames '
              'at 0.5 tau intervals over 1000 tau, all 4000 atoms, variables '
              '`id`, `type`, `unwrapped_coordinates`, `cell_lengths`, and `time`. '
              'Coordinates are **unwrapped**; wrap each frame into its periodic '
              'cubic cell for Voronoi analysis. Match atom types to the state '
              'README or `data.mixture`. For example:', '',
              '```python', 'from netCDF4 import Dataset',
              'with Dataset("B3_baseline_3C/replica_1/positions.nc") as nc:',
              '    xyz = nc.variables["unwrapped_coordinates"][0]  # atom, xyz',
              '    box = nc.variables["cell_lengths"][0]',
              '    wrapped = xyz % box', '```', '',
              'MDAnalysis, ASE, and OVITO can also inspect/convert trajectories; '
              'confirm their AMBER NetCDF atom ordering, units, and periodic cell '
              'handling before computing topology.', '',
              'LJ reduced units: sigma=1, mass=1, primary solvent epsilon_1=1; '
              'time tau=sigma sqrt(m/epsilon_1), T*=2, rho*=0.7, '
              'diffusion sigma²/tau, entropy kB per particle. '
              'The box is periodic and cubic with L=(4000/0.7)^(1/3). '
              'LAMMPS type ids start at 1 in selected.json species order. '
              'L5 splits the epsilon=0.25 solvent into two mechanically identical '
              'labels. Merge mechanically identical labels when the estimator is '
              'label-sensitive; report merged and unmerged estimates.', '',
              '## Reference values and plot', '',
              '`s_ex_TI` is MD thermodynamic-integration excess entropy per '
              'particle (19 density points, configurational-pressure ideal '
              'reference, screened Debye–Hückel dilute endpoint); it excludes '
              'ideal mixing entropy. `D_net`, `D_probe`, and species D use '
              'whole-system COM displacement in the MSD, with 20–100 tau as '
              'the primary fit and 10–50 / 50–200 tau probe alternatives. '
              '`s2_MD` is all-species MD pair entropy with the fine-core '
              'resolved 20-block jackknife hybrid at width 0.005, exactly the '
              'field used by `doc/component_study_pair_entropy_states.tex`; '
              'it is an approximation to TI excess entropy. `s_ex_HNC` is the '
              'refined HNC prediction. Sources: study `manifest.json` and '
              '`replica_result.json`, `MDPairEntropyAssessment/refined/` '
              'and `summarize.py`; see the source paths in each reference.json.', '',
              'To plot external estimates in a CSV with `code,replica,s_est` '
              '(replica optional):', '',
              '```bash',
              'python3 plot_diffusivity_vs_entropy.py --entropy external.csv',
              'python3 plot_diffusivity_vs_entropy.py --entropy reference_values.csv --column s2_MD',
              '```', '',
              'The plots compare ln D with the supplied estimate and the TI '
              'reference and print fitted descriptive Rosenfeld slopes. '
              'L5 is a hollow label-null duplicate. Confidence bars reflect '
              'only three independent replicas, without finite-size or estimator '
              'bias corrections. Extended N=4000 and N=8000 trajectories in '
              '`ComponentEntropyStudy_20261030/followup_20261101` are in '
              'progress and can be added later.', '']
    return '\n'.join(lines)


def main():
    selected = {x['name']: x for x in json.loads((STUDY/'selected.json').read_text())}
    OUT.mkdir(exist_ok=True)
    copied = []
    rows = []
    grouped = {}
    for code, name in ORDER:
        item = selected[name]
        grouped[code] = []
        for rep in (1, 2, 3):
            src = STUDY/'md'/name/f'replica_{rep}'
            dest = OUT/f'{code}_{name}'/f'replica_{rep}'
            manifest = json.loads((src/'manifest.json').read_text())
            result = json.loads((src/'replica_result.json').read_text())
            pair_path = PAIR/'refined'/f'{code}_r{rep}.json'
            pair = json.loads(pair_path.read_text())
            assert pair['original_core_histograms_match_exactly']
            assert manifest['composition'] == item and manifest['N'] == 4000
            assert result['replica'] == rep and result['seed'] == manifest['seed']
            transport = result['transport']
            labels = [s['label'] for s in item['species']]
            windows = transport['fitting_windows']
            probe_indices = [labels.index(x) + 1 for x in ('probe_plus2', 'probe_minus2')]
            row = dict(code=code, name=name, replica=rep, seed=manifest['seed'],
                       N=4000, n_labels=len(labels),
                       n_distinct=item['mechanically_distinct_components'],
                       I=item['I'], s_ex_TI=result['s_ex'],
                       s_ex_HNC=result['HNC_s_ex'],
                       s2_MD=pair['jackknife20_hybrid']['0.005']['s2'],
                       D_net=transport['D_net'], D_probe=transport['D_probe'],
                       D_probe_10_50=float(np.mean(np.asarray(windows['10_50'])[probe_indices])),
                       D_probe_50_200=float(np.mean(np.asarray(windows['50_200'])[probe_indices])))
            for label in LABELS:
                row[f'D_{label}'] = transport['D_species'].get(label, '')
            assert all(np.isclose(windows['20_100'][j+1], transport['D_species'][label])
                       for j, label in enumerate(labels))
            sources = {}
            for filename in ('positions.nc', 'in.run', 'data.mixture', 'log.lammps'):
                source = src/'transport'/filename
                target = dest/filename
                sha = copy_checked(source, target)
                source_rel = str(source.relative_to(ROOT))
                sources[filename] = {'path': source_rel, 'sha256': sha}
                copied.append({'path': str(target.relative_to(OUT)),
                               'source': source_rel, 'sha256': sha})
            reference = dict(row)
            reference['sources'] = sources
            reference['source_manifest'] = str((src/'manifest.json').relative_to(ROOT))
            reference['source_result'] = str((src/'replica_result.json').relative_to(ROOT))
            reference['source_pair_entropy'] = str(pair_path.relative_to(ROOT))
            create(dest/'reference.json', json.dumps(reference, indent=2, allow_nan=False)+'\n')
            rows.append(row)
            grouped[code].append(row)
        create(OUT/f'{code}_{name}'/'README.md', per_state_readme(code, item, grouped[code]))
    columns = ['code', 'name', 'replica', 'seed', 'N', 'n_labels', 'n_distinct', 'I',
               's_ex_TI', 's_ex_HNC', 's2_MD', 'D_net', 'D_probe'] + [f'D_{x}' for x in LABELS] + [
                   'D_probe_10_50', 'D_probe_50_200']
    create(OUT/'reference_values.csv', csv_text(rows, columns))
    means = []
    metrics = ['s_ex_TI', 's2_MD', 'D_net', 'D_probe', 'D_plus2', 'D_minus2']
    for code, name in ORDER:
        line = {'code': code, 'name': name}
        for metric in metrics:
            source_metric = {'D_plus2': 'D_probe_plus2',
                             'D_minus2': 'D_probe_minus2'}.get(metric, metric)
            line[f'{metric}_mean'], line[f'{metric}_ci95_halfwidth'] = mean_ci(
                [r[source_metric] for r in grouped[code]])
        means.append(line)
    mean_cols = ['code', 'name'] + [f'{m}_{suffix}' for m in metrics
                                   for suffix in ('mean', 'ci95_halfwidth')]
    create(OUT/'reference_state_means.csv', csv_text(means, mean_cols))
    create(OUT/'README.md', top_readme(selected, grouped, means))
    create(OUT/'plot_diffusivity_vs_entropy.py', PLOT_SCRIPT)
    create(OUT/'MANIFEST.json', json.dumps({'copied_files': copied}, indent=2)+'\n')
    create(OUT/'checksums.sha256', ''.join(f'{r["sha256"]}  {r["path"]}\n'
                                            for r in copied))
    print(f'Packaged {len(rows)} replicas, {len(copied)} verified copied files in {OUT}')


PLOT_SCRIPT = '''#!/usr/bin/env python3
"""Compare supplied structural entropy estimates with MD diffusion and TI."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t

HERE = Path(__file__).resolve().parent
ORDER = ['B3', 'U4', 'M5', 'H5', 'L5', 'B5', 'U6', 'H7']
RED = {'B5', 'U6', 'H7'}
OBS = [('D_net', 'net'), ('D_probe', 'probe'),
       ('D_probe_plus2', '+2'), ('D_probe_minus2', '-2')]


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def estimate(values):
    x = np.asarray(values, float)
    return x.mean(), t.ppf(.975, len(x)-1)*x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else 0.


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--entropy', type=Path, required=True)
    p.add_argument('--column', default='s_est')
    p.add_argument('--label', default=None)
    args = p.parse_args()
    refs = read_csv(HERE/'reference_values.csv')
    supplied = read_csv(args.entropy)
    if not supplied or not {'code', args.column} <= supplied[0].keys():
        p.error(f'entropy CSV must contain code and {args.column}')
    has_replica = 'replica' in supplied[0]
    by_code = {}
    for code in ORDER:
        rr = [r for r in refs if r['code'] == code]
        ee = [r for r in supplied if r['code'] == code]
        if not ee:
            continue
        if has_replica:
            emap = {int(r['replica']): float(r[args.column]) for r in ee}
            if set(emap) != {1, 2, 3} or len(ee) != 3:
                p.error(f'{code}: expected exactly replicas 1, 2, 3')
            ev = [emap[int(r['replica'])] for r in rr]
        else:
            if len(ee) != 1:
                p.error(f'{code}: expected one state estimate')
            ev = [float(ee[0][args.column])]
        by_code[code] = {'est': estimate(ev), 'ti': estimate([r['s_ex_TI'] for r in rr]),
                         'diff': {key: estimate([np.log(float(r[key])) for r in rr]) for key, _ in OBS}}
    if len(by_code) < 2:
        p.error('need at least two matching states')
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    entropy_label = args.label or args.column
    for ax, (key, label) in zip(axes.flat[:4], OBS):
        for code, d in by_code.items():
            color = 'tab:red' if code in RED else 'tab:blue'
            markerface = 'none' if code == 'L5' else color
            for xkey, marker, alpha in [('est', 'o', 1), ('ti', 's', .55)]:
                x, xerr = d[xkey]
                y, yerr = d['diff'][key]
                ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=marker,
                            markerfacecolor=markerface if xkey == 'est' else 'none',
                            markeredgecolor=color, ecolor=color, alpha=alpha, capsize=2)
                ax.annotate(code, (x, y), xytext=(3, 3), textcoords='offset points', fontsize=7)
        ax.set_xlabel(f'{entropy_label} (circle); TI s_ex (square) / kB')
        ax.set_ylabel(f'ln D_{label}')
        for xkey in ('est', 'ti'):
            fit_codes = [c for c in by_code if c != 'L5']
            x = np.array([by_code[c][xkey][0] for c in fit_codes])
            y = np.array([by_code[c]['diff'][key][0] for c in fit_codes])
            slope = np.polyfit(x, y, 1)[0]
            print(f'{key} vs {xkey}: Rosenfeld B = {slope:.6g} (n={len(x)})')
    ax = axes.flat[4]
    for code, d in by_code.items():
        color = 'tab:red' if code in RED else 'tab:blue'
        ax.errorbar(d['ti'][0], d['est'][0], xerr=d['ti'][1], yerr=d['est'][1],
                    fmt='o', color=color, markerfacecolor='none' if code == 'L5' else color,
                    capsize=2)
        ax.annotate(code, (d['ti'][0], d['est'][0]), xytext=(3, 3),
                    textcoords='offset points', fontsize=7)
    fit_codes = [c for c in by_code if c != 'L5']
    x = np.array([by_code[c]['ti'][0] for c in fit_codes])
    y = np.array([by_code[c]['est'][0] for c in fit_codes])
    b, a = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 100)
    ax.plot(xx, a+b*xx, 'k--', label=f'linear fit: y={a:.3f}+{b:.3f}x')
    ax.legend(fontsize=8)
    ax.set_xlabel('TI s_ex / kB')
    ax.set_ylabel(f'{entropy_label} / kB')
    axes.flat[5].axis('off')
    axes.flat[5].text(.03, .9, 'Blue: no monovalent background\\nRed: with ±1 background\\n'
                      'Hollow L5: label-only duplicate\\nBars: 95% Student-t CI', va='top')
    for ext in ('pdf', 'png'):
        target = HERE/f'diffusivity_vs_entropy.{ext}'
        fig.savefig(target, dpi=200)
        print(target)


if __name__ == '__main__':
    main()
'''

if __name__ == '__main__':
    main()
