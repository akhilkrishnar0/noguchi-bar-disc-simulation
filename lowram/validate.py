"""Machine-readable science gate for the low-N exploratory run."""
import json
from pathlib import Path
import numpy as np
p=Path(__file__).parent
f=np.genfromtxt(p/'check2/diagnostics.csv',delimiter=',',names=True)
g=np.genfromtxt(p/'fine/diagnostics.csv',delimiter=',',names=True)
height_percent=100*(f['zrms'][-1]/f['zrms'][0]-1)
relative={k:float(100*(g[k][-1]-f[k][-1])/abs(f[k][-1])) for k in ('zrms','meanR','A2')}
criteria={'height_stationarity_abs_percent_lt_5':abs(height_percent)<5,
          'short_timestep_height_agreement_abs_percent_lt_1':abs(relative['zrms'])<1,
          'short_timestep_radius_agreement_abs_percent_lt_1':abs(relative['meanR'])<1}
result={'status':'FAILED_EQUILIBRIUM' if not all(criteria.values()) else 'SHORT_SCREEN_PASSED_NOT_VALIDATED',
        'early_height_change_percent':float(height_percent),
        'fine_vs_regular_at_75myr_percent':relative,'criteria':criteria,
        'interpretation':'Do not attribute bar or GLSB formation to an unvalidated, strongly thickening galaxy. A2 alone is not proof of a persistent stellar bar.',
        'limitations':['1500 particles leads to shot noise and enhanced collisional relaxation',
                       'approximate circular/Jeans velocities, no exact distribution function',
                       'vertical moment matching is necessary but not sufficient',
                       'no full force/softening/seed convergence or phase-coherent bar identification',
                       'imposed bar injects energy and angular momentum']}
(p/'equilibrium_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
