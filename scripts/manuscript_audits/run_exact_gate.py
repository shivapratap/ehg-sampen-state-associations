import sys, json, time
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
src=ROOT/'src'; sys.path.insert(0,str(src))
import models
def run(frame, features):
    out=[]
    for f in features:
        rows,audit=models._fit_one(frame,f); out += rows
    return pd.DataFrame(out)
fwh=ROOT/'work/primary/event_unit_analysis_60s_EHG9.csv'; frozen=ROOT/'work/primary/mixed_model_results_60s_fwh_EHG9.csv'
gi=ROOT/'work/gi_features_60s_EHG9.csv'; gifrozen_path=ROOT/'results/gi/tables/gi_model_results.csv'
start=time.time(); df=pd.read_csv(fwh); fr=pd.read_csv(frozen); feats=[x for x in fr.feature.dropna().unique() if x in df.columns]
got=run(df,feats); joined=got.merge(fr[['feature','contrast','estimate','p_value']],on=['feature','contrast'],suffixes=('_new','_saved'))
joined['estimate_abs_diff']=(joined.estimate_new-joined.estimate_saved).abs(); joined['p_abs_diff']=(joined.p_value_new-joined.p_value_saved).abs(); joined.to_csv(ROOT/'results/manuscript/audits/exact_gate_fwh_detail.csv',index=False)
dg=pd.read_csv(gi); gifrozen=pd.read_csv(gifrozen_path); gf=gifrozen[gifrozen.feature.isin(['gi_total','gi_sd'])]; gg=run(dg,['gi_total','gi_sd']); jg=gg.merge(gf[['feature','contrast','estimate','p_value']],on=['feature','contrast'],suffixes=('_new','_saved')); jg['estimate_abs_diff']=(jg.estimate_new-jg.estimate_saved).abs(); jg['p_abs_diff']=(jg.p_value_new-jg.p_value_saved).abs(); jg.to_csv(ROOT/'results/manuscript/audits/exact_gate_gi_detail.csv',index=False)
summary={'fwh_features':len(feats),'fwh_rows':len(joined),'fwh_max_estimate_abs_diff':float(joined.estimate_abs_diff.max()),'fwh_max_p_abs_diff':float(joined.p_abs_diff.max()),'gi_rows':len(jg),'gi_max_estimate_abs_diff':float(jg.estimate_abs_diff.max()),'gi_max_p_abs_diff':float(jg.p_abs_diff.max()),'pass':bool(joined.estimate_abs_diff.max()<=1e-8 and joined.p_abs_diff.max()<=1e-8 and jg.estimate_abs_diff.max()<=1e-8 and jg.p_abs_diff.max()<=1e-8),'wall_seconds':time.time()-start}
(ROOT/'results/manuscript/audits/exact_gate_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
