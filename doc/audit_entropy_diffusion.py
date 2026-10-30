#!/usr/bin/env python3
"""Audit the reported manuscript values; does not re-estimate trajectory errors.

Run from any directory with Python, NumPy and Matplotlib. Reads the archived
Table II CSV and the Figure 5 solvent constants; Table IV GK values below are
transcribed explicitly. Writes JSON and a replacement Figure 5 in doc/.
"""
from pathlib import Path
import csv
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

DOC = Path(__file__).resolve().parent
ROOT = DOC.parent


def fit(x, d):
    x, y = np.asarray(x), np.log(d)
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope*x + intercept)
    return dict(B=float(slope), intercept=float(intercept),
                R2=float(1 - residual@residual / ((y-y.mean())@(y-y.mean()))), n=len(x))


def main():
    with (ROOT/'lammps_study/results/rosenfeld_diffusion_coefficients.csv').open() as f:
        rows = list(csv.DictReader(f))
    s = np.array([float(r['s_ex_MD']) for r in rows])
    h = np.array([float(r['s_ex_HNC']) for r in rows])
    d = np.array([float(r['D_net']) for r in rows])
    ion = np.array([float(r['D_z2']) for r in rows])
    # Reported in Table IV, in the same system order as the CSV.
    gk_ion = np.array([.2038, .2017, .2057, .2037])
    msd_ion_2pt = np.array([.2027, .1813, .2072, .2036])
    compositions = [[.94,.03,.03], [.89,.05,.03,.03],
                    [.89,.05,.03,.03], [.84,.05,.03,.05,.03]]
    mixing = [-float(np.dot(x, np.log(x))) for x in compositions]
    # Read exactly the constants used in the existing figure generator.
    source = (DOC/'generate_manuscript_figures.py').read_text()
    block = source.split('    solvents = [',1)[1].split('    ]',1)[0]
    solvent = re.findall(r'\("[^"]+",\s*([-\d.]+),\s*([\d.]+) \*', block)
    solvent = np.asarray(solvent, dtype=float)
    results = {
        'provenance': 'Reported values only; no new simulation or uncertainty estimate.',
        'net_MD_three_points': fit(s[:3], d[:3]),
        'net_HNC_three_points': fit(h[:3], d[:3]),
        'net_MD_four_points': fit(s,d),
        'solvent_pair_entropy_six_points': fit(solvent[:,0], solvent[:,1]),
        'ideal_mixing_entropy_per_kB': dict(zip([r['System'] for r in rows], mixing)),
        'five_minus_three_total_entropy_per_kB': float(mixing[3]-mixing[0]+s[3]-s[0]),
        'ion_5C_over_3C_table_II': float(ion[3]/ion[0]),
        'ion_5C_over_3C_table_IV_MSD': float(msd_ion_2pt[3]/msd_ion_2pt[0]),
        'ion_5C_over_3C_table_IV_GK': float(gk_ion[3]/gk_ion[0]),
    }
    (DOC/'entropy_diffusion_audit.json').write_text(json.dumps(results, indent=2)+'\n')
    scale = .7**(1/3)/np.sqrt(2.)
    colors = ['#1f77b4','#2ca02c','#ff7f0e','#d62728']
    labels = ['3C',r'4C ($\epsilon_2=1.25$)',r'4C ($\epsilon_2=0.75$)','5C']
    fig, axes = plt.subplots(1,2,figsize=(11.8,4.8),layout='constrained')
    for i, (c, label) in enumerate(zip(colors,labels)):
        axes[0].scatter(s[i],d[i]*scale,color=c,label=label)
        axes[1].scatter(s[i],ion[i]*scale,color=c,marker='o')
        axes[1].scatter(s[i],gk_ion[i]*scale,color=c,marker='s',facecolors='none')
    p = results['net_MD_three_points']
    xx = np.linspace(s.min()-.002,s.max()+.002,100)
    axes[0].plot(xx,np.exp(p['intercept']+p['B']*xx)*scale,'k--',
                 label=f"3-point fit: B={p['B']:.2f}, $R^2$={p['R2']:.2f}")
    axes[1].scatter([],[],color='black',marker='o',label='MSD, Table II')
    axes[1].scatter([],[],edgecolors='black',facecolors='none',marker='s',label='GK, Table IV')
    for ax in axes:
        ax.set_xlabel(r'Reported MD $s_{\rm ex}/k_B$')
        ax.set_ylabel(r'$D\rho^{1/3}/\sqrt{k_BT/m}$')
        ax.set_yscale('log')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    axes[0].set_title('(a) Mean self-diffusion')
    axes[1].set_title(r'(b) Mean self-diffusion of $|z|=2$ ions')
    fig.savefig(DOC/'figures/rosenfeld_scaling_review.png',dpi=250)
    plt.close(fig)
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    main()
