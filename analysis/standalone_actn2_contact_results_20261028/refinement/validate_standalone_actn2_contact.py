#!/usr/bin/env python3
"""Multi-seed and timestep checks for the standalone experiment (run in Slurm).

Uses independent random streams at each resolution. Six-standard-error
comparisons detect large resolution effects; they are not equivalence tests.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import standalone_actn2_contact as m


def compare(a, b, name, tolerance=0.):
    a, b = np.asarray(a, float), np.asarray(b, float)
    error = float(abs(a.mean()-b.mean()))
    se = float(np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)))
    return dict(metric=name, coarse_mean=float(a.mean()), refined_mean=float(b.mean()),
                absolute_difference=error, standard_error_of_difference=se,
                tolerance=tolerance, passed=error <= 6*se+tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replicates', type=int, default=256)
    args = parser.parse_args()
    if args.replicates < 64:
        parser.error('Use at least 64 realizations per seed.')
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    model = m.Model()
    results = dict(failure=[], weak=[], energy=[], comparisons=[])
    samples = {r:[] for r in (0,1)}
    weak_samples = {(kind,r):[] for kind in ('barbed','inner') for r in (0,1)}
    for seed in (20261028,20261029,20261030):
        for resolution, dt in enumerate((1.e-7,5.e-8)):
            print(f'Failure seed {seed}, dt {dt:g}', flush=True)
            data, summary = m.simulate_rig(model,40.,duration=.002,dt=dt,
                replicates=args.replicates,failure=True,seed=seed+10000*resolution)
            results['failure'].append(dict(seed=seed+10000*resolution,**summary))
            samples[resolution].append(dict(alive=data['alive'][-1],
                penetration=data['penetration'][-1]*1000))
        for kind in ('barbed','inner'):
            for resolution, dt in enumerate((1.e-5,2.5e-6)):
                print(f'Weak {kind}, seed {seed}, dt {dt:g}', flush=True)
                run = m.weak_sites(model,kind,replicates=4096,dt=dt,
                                  seed=seed+10000*resolution)
                results['weak'].append(dict(seed=seed+10000*resolution,**run['summary']))
                weak_samples[kind,resolution].append(run['alive'])
    coarse = np.concatenate([s['alive'] for s in samples[0]])
    fine = np.concatenate([s['alive'] for s in samples[1]])
    for metric, a, b in (
        ('any link loss',~coarse.all(axis=1),~fine.all(axis=1)),
        ('all links lost',~coarse.any(axis=1),~fine.any(axis=1)),
        ('mean live links',coarse.sum(axis=1),fine.sum(axis=1)),
        ('final penetration nm',np.concatenate([s['penetration'] for s in samples[0]]),
         np.concatenate([s['penetration'] for s in samples[1]]))):
        results['comparisons'].append(compare(a,b,metric,tolerance=1. if 'nm' in metric else 0.))
    for kind in ('barbed','inner'):
        results['comparisons'].append(compare(np.concatenate(weak_samples[kind,0]),
            np.concatenate(weak_samples[kind,1]),f'{kind} weak survival'))
    for dt in (1.e-7,5.e-8,2.5e-8):
        print(f'Energy balance dt {dt:g}',flush=True)
        _, summary = m.simulate_rig(model,20.,dt=dt,nquad=32)
        results['energy'].append(dict(dt_s=dt,
            residual_pNnm=summary['work_balance_residual_pNnm']))
    residuals = np.abs([s['residual_pNnm'] for s in results['energy']])
    results['energy_refinement_passed'] = bool(np.all(residuals[1:] < .7*residuals[:-1]))
    results['passed'] = all(c['passed'] for c in results['comparisons']) and results['energy_refinement_passed']
    results['scope'] = ('Supported translational toy rig and fixed-actin weak sites. '
        'Independent streams, three seeds per resolution; 6-SE comparisons, '
        'plus 1 nm tolerance for mean penetration. A pass excludes large observed '
        'resolution effects at this sample size; it does not establish equivalence '
        'or production-network stability.')
    results['wall_seconds'] = time.perf_counter()-start
    results['source_sha256'] = hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()
    (args.output/'convergence.json').write_text(json.dumps(results,indent=2))
    lines = ['# Additional standalone validation','',results['scope'],'',
        '| Metric | Coarse | Refined | Difference | SE of difference | Passed |',
        '|---|---:|---:|---:|---:|---|']
    for c in results['comparisons']:
        lines.append(f'| {c["metric"]} | {c["coarse_mean"]:.5g} | {c["refined_mean"]:.5g} | '
            f'{c["absolute_difference"]:.3g} | {c["standard_error_of_difference"]:.3g} | {c["passed"]} |')
    lines += ['',f'Energy-balance residuals (pN nm), dt / dt/2 / dt/4: {residuals.tolist()}.',
              '',f'Overall checks passed: {results["passed"]}. Wall time: {results["wall_seconds"]:.1f} s.']
    (args.output/'convergence.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(passed=results['passed'],comparisons=results['comparisons'])),flush=True)
    if not results['passed']:
        raise SystemExit('Resolution checks failed; inspect convergence.json.')


if __name__ == '__main__':
    main()
