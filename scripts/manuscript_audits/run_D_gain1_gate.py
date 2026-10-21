from pathlib import Path
import sys,json,time
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/'src')); from features import _window_features
z=np.load(ROOT/'results/manuscript/audits/waveforms/audited_windows_EHG9_fwh.npz',allow_pickle=True); X=z['signals']; B=z['broad_signals']; meta=z['metadata']; ids=[m['window_id'] for m in meta]
cached=pd.read_csv(ROOT/'work/primary/event_unit_analysis_60s_EHG9.csv'); cached=cached[cached.window_id.isin(ids)].set_index('window_id'); feats=[f for f in pd.read_csv(ROOT/'work/primary/mixed_model_results_60s_fwh_EHG9.csv').feature.unique() if f in cached.columns]
rows=[]; start=time.time()
for i,wid in enumerate(ids):
    vals,_=_window_features(X[i],B[i])
    for f in feats: rows.append({'window_id':wid,'feature':f,'computed':vals[f],'cached':cached.loc[wid,f],'abs_diff':abs(vals[f]-cached.loc[wid,f])})
r=pd.DataFrame(rows); r['pass']=np.isclose(r.computed,r.cached,atol=1e-10,rtol=1e-8); out=ROOT/'results/manuscript/audits/D'; r.to_csv(out/'gain1_gate_detail.csv',index=False); s={'status':'PASS' if bool(r['pass'].all()) else 'FAIL','windows':len(ids),'features':len(feats),'rows':len(r),'max_abs_diff':float(r.abs_diff.max()),'n_fail':int((~r['pass']).sum()),'wall_seconds':time.time()-start}; (out/'gain1_gate_summary.json').write_text(json.dumps(s,indent=2)); print(json.dumps(s,indent=2))
