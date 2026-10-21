from pathlib import Path
import sys,time,json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2]; src=ROOT/'src'; sys.path.insert(0,str(src)); from features import _window_features; import config
sys.path.insert(0,str(ROOT/'src')); from gi_sampen.gi_profile import gi_summaries_rounded
z=np.load(ROOT/'results/manuscript/audits/waveforms/audited_windows_EHG9_fwh.npz',allow_pickle=True); X=z['signals']; B=z['broad_signals']; meta=z['metadata']; gains=[.25,.5,1,1.3,2,4,7.1]; out=ROOT/'results/manuscript/audits/D'; out.mkdir(exist_ok=True)
base=[]; start=time.time()
for i in range(60):
 v,_=_window_features(X[i],B[i]); g=gi_summaries_rounded(X[i]); base.append({**v,'gi_total':g['gi_total'],'gi_sd':g['gi_sd'],'window_id':meta[i]['window_id']})
base=pd.DataFrame(base); feats=[f for f in pd.read_csv(ROOT/'work/primary/mixed_model_results_60s_fwh_EHG9.csv').feature.unique() if f in base.columns]+['gi_total','gi_sd']
rows=[]
for i in range(60):
 for gain in gains:
  v,_=_window_features(X[i]*gain,B[i]*gain); g=gi_summaries_rounded(X[i]*gain); vals={**v,'gi_total':g['gi_total'],'gi_sd':g['gi_sd']}
  for f in feats:
   expected=base.loc[i,f]*(gain if ('mean_absolute_value' in f or 'root_mean_square' in f or 'peak_to_peak' in f or 'waveform_length' in f) else 1.0); obs=vals[f]; err=abs(obs-expected) if np.isfinite(obs) and np.isfinite(expected) else np.nan; rows.append({'window_id':meta[i]['window_id'],'feature':f,'gain':gain,'observed':obs,'expected':expected,'abs_error':err,'pass':bool(np.isfinite(err) and err<=1e-10+1e-8*abs(expected))})
res=pd.DataFrame(rows); res.to_csv(out/'gain_audit_window_results.csv',index=False); res.groupby(['feature','gain'],as_index=False).agg(finite_count=('observed',lambda x:int(np.isfinite(x).sum())),max_abs_error=('abs_error','max'),pass_count=('pass','sum')).to_csv(out/'gain_audit_summary.csv',index=False)
laws=pd.DataFrame([{'feature':f,'law':('degree 1' if any(x in f for x in ['mean_absolute_value','root_mean_square','peak_to_peak','waveform_length']) else 'degree 0'),'source':'feature definition / empirical audit'} for f in feats]); laws.to_csv(out/'feature_gain_laws.csv',index=False)
(out/'execution_record.json').write_text(json.dumps({'status':'complete','n_windows':60,'gains':gains,'features':len(feats),'wall_seconds':time.time()-start},indent=2)); print('done',len(rows),time.time()-start)
