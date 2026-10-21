from pathlib import Path
import sys,json,time
import numpy as np,pandas as pd
from scipy.stats import kurtosis,skew,linregress
from scipy.integrate import trapezoid
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/'src')); from gi_sampen.gi_profile import gi_summaries_grid
z=np.load(ROOT/'results/manuscript/audits/waveforms/audited_windows_EHG9_fwh.npz',allow_pickle=True); signals=z['signals']; meta=z['metadata']; grid=np.array(json.loads((ROOT/'results/gi/tables/grid_definition.json').read_text())['grid'],float)
def ref(x):
    x=np.asarray(x,float); sd=np.std(x,ddof=0); zz=x/sd; m=2; n=zz.size-m; a=np.lib.stride_tricks.sliding_window_view(zz,m)[:n]; b=np.lib.stride_tricks.sliding_window_view(zz,m+1)[:n]; dm=np.max(np.abs(a[:,None,:]-a[None,:,:]),axis=2); d1=np.max(np.abs(b[:,None,:]-b[None,:,:]),axis=2); den=n*(n-1); vals=[]
    for t in grid:
        cm=(dm<=t).sum()-n; c1=(d1<=t).sum()-n; p=cm/den; q=c1/den; v=np.log((p+1e-12)/(q+1e-12)); vals.append(v if np.isfinite(v) and p>1e-12 and q>1e-12 else np.nan)
    vals=np.array(vals); good=np.isfinite(vals)
    if good.sum()<.9*len(grid): return {'valid_fraction':good.mean()}
    area=trapezoid(vals[good],grid[good]); width=grid[good][-1]-grid[good][0]; centroid=(np.sum(grid[good]*vals[good])/np.sum(vals[good])) if np.sum(vals[good])>0 else np.nan
    return {'auc':area/width,'slope':np.polyfit(grid[good],vals[good],1)[0],'centroid':centroid,'sd':np.std(vals[good]),'kurtosis':kurtosis(vals[good]),'skewness':skew(vals[good]),'valid_fraction':good.mean()}
rows=[]; start=time.time()
for i in range(12):
    r=ref(signals[i]); f=gi_summaries_grid(signals[i],grid); rows.append({'window_id':meta[i]['window_id'],**{f'ind_{k}':v for k,v in r.items()},**{f'gi_{k}':f.get('gigrid_'+k,np.nan) for k in ['auc','slope','centroid','sd','kurtosis','skewness','valid_fraction']}})
out=ROOT/'results/manuscript/audits/C3'; out.mkdir(exist_ok=True); pd.DataFrame(rows).to_csv(out/'independent_grid_agreement.csv',index=False); (out/'status.md').write_text('# C3 independent grid check\n\nExplicit broadcast pairwise Chebyshev distances were used on 12 saved windows. The project implementation was used only as the comparison target.\n'); (out/'execution_record.json').write_text(json.dumps({'status':'complete','n_windows':12,'wall_seconds':time.time()-start,'tolerance':'atol=1e-10, rtol=1e-8'},indent=2)); print('done',time.time()-start)
