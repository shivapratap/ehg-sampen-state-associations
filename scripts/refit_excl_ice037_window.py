import sys; sys.path.insert(0,'src')
import pandas as pd, numpy as np
from reml import fit_table, load_windows
df=load_windows('results/conditional/windows_gi_rms_frozen_full.csv')
W='ice037_p_2of2|W60|6094-18094'
print('window present:', (df.window_id==W).sum(), df[df.window_id==W][['rms','gi_total','frozen_total','window_class']].to_dict('records'))
print('rms rank/max:', df.rms.max(), df.rms.nlargest(3).values)
rows=[]
for label,d in [('all',df),('excl_ice037_window',df[df.window_id!=W])]:
    print(label,len(d),d.patient_id.nunique(),d.recording_id.nunique(),d.window_class.value_counts().to_dict())
    for o in ['frozen_total','gi_total']:
        for a in (False,True):
            t=fit_table(d,o,a); t.insert(0,'set',label); rows.append(t)
out=pd.concat(rows)
s=out[out.term.isin(['UC_minus_IN','FM_minus_IN'])][['set','outcome','adjusted_log_RMS','term','estimate','se','ci_low','ci_high','p']]
print(s.to_string(float_format=lambda v:f'{v:.4f}'))
out.to_csv('results/conditional/refit_excl_ice037_window.csv',index=False)
