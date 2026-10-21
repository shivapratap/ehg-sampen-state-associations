#!/usr/bin/env python3
"""Secondary fetal-movement evidence consolidation (Phase E3)."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "fm"
RESULTS = ROOT / "results" / "fm"
FEATURES = ROOT / "work" / "primary" / "confirmatory_feature_registry_60s_fwh.csv"
LONG = ROOT / "work" / "feature_qc" / "window_features_60s_long.csv"

def _bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]; out=np.empty(len(p)); out[o]=np.clip(q,0,1); return out

def run_paired_uc_vs_fm():
    reg=pd.read_csv(FEATURES); fs=reg.loc[reg.primary_confirmatory.astype(bool),'feature'].tolist()
    if len(fs)!=34: raise ValueError(f'Expected 34 features, found {len(fs)}')
    d=pd.read_csv(LONG,usecols=['patient_id','channel','window_class']+fs); d=d[(d.channel=='EHG9')&d.window_class.isin(['contraction','fetal_movement'])]
    rows=[]; pop=[]
    for f in fs:
        m=d.groupby(['patient_id','window_class'])[f].median().unstack().dropna(subset=['contraction','fetal_movement']); diff=m.contraction-m.fetal_movement
        stat,p=wilcoxon(m.contraction,m.fetal_movement,zero_method='wilcox') if len(m) else (np.nan,np.nan)
        rows.append({'feature':f,'n_paired_women':len(m),'median_uc_minus_fm':diff.median(),'wilcoxon_statistic':stat,'p_value':p})
        pop.append({'feature':f,'n_paired_women':len(m),'paired_women':';'.join(m.index.astype(str))})
    r=pd.DataFrame(rows); r['q_value']=_bh(r.p_value); r['supported']=r.q_value<.05; r.to_csv(RESULTS/'paired_uc_vs_fm_34_features.csv',index=False); pd.DataFrame(pop).to_csv(RESULTS/'paired_uc_vs_fm_population.csv',index=False); return r

def run_fm_vs_in():
    src=ROOT/'work'/'primary'/'mixed_model_results_60s_fwh_EHG9.csv'; d=pd.read_csv(src); d=d[d.contrast.eq('fetal_movement_vs_baseline')].copy(); d.to_csv(RESULTS/'fm_vs_in_34_features.csv',index=False); return d

def run_clean_fm_sensitivity():
    src=ROOT/'results'/'robustness'/'overlap'/'reported_fm_feature_stability.csv'; d=pd.read_csv(src); d.to_csv(RESULTS/'clean_fm_sensitivity_summary.csv',index=False); return d

def summarize_fm_evidence(paired, fm, clean):
    rows=[]
    for _,r in paired.iterrows(): rows.append({'feature':r.feature,'analysis':'paired UC-versus-FM','effect_direction':'positive' if r.median_uc_minus_fm>0 else 'negative','estimate':r.median_uc_minus_fm,'p_value':r.p_value,'q_value':r.q_value,'supported':r.supported,'sensitivity_status':'paired','interpretation':'physiological differentiation' if r.supported else 'not supported'})
    profile={'fwh_entropy_mediansampen','fwh_entropy_sdsampen','fwh_entropy_kurtosissampen','fwh_entropy_skewnesssampen'}
    for _,r in fm[fm.feature.isin(profile)].iterrows(): rows.append({'feature':r.feature,'analysis':'FM versus IN','effect_direction':'positive' if r.estimate>0 else 'negative','estimate':r.estimate,'p_value':r.p_value,'q_value':r.q_value,'supported':r.primary_validated,'sensitivity_status':'clean FM overlap available','interpretation':'SampEn profile FM finding'})
    out=pd.DataFrame(rows); out.to_csv(RESULTS/'fm_evidence_summary.csv',index=False); return out

def main():
    WORK.mkdir(parents=True,exist_ok=True); RESULTS.mkdir(parents=True,exist_ok=True)
    p=run_paired_uc_vs_fm(); f=run_fm_vs_in(); c=run_clean_fm_sensitivity(); summarize_fm_evidence(p,f,c)
    print(f'Paired women: {p.n_paired_women.min()}; supported: {int(p.supported.sum())}/34')

if __name__=='__main__': main()
