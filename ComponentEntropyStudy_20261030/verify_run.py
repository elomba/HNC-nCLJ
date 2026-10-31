#!/usr/bin/env python3
"""Fail a sequential allocation before it proceeds past an incomplete stage."""
from pathlib import Path
import sys
p=Path(sys.argv[1]);log=(p/'lammps_run.log').read_text(errors='replace')
assert 'Total wall time:' in log and 'ERROR' not in log,p
assert (p/'final.restart').is_file() and (p/'data.final').is_file(),p
if p.name=='transport':
    import numpy as np
    import json
    m=json.loads((p.parent/'manifest.json').read_text())
    a=np.loadtxt(p/'msd_common.dat')
    x=np.array([s['x'] for s in m['composition']['species']])
    assert np.allclose(a[:,1],a[:,2:]@x,rtol=5e-5,atol=2e-5),p
import json
(p/'verified.json').write_text(json.dumps(dict(normal_completion=True))+'\n')
print('Verified',p,flush=True)
