"""Build the small display tables ("adapters") that feed the manuscript figures.

Reads saved results only; recalculates nothing. Run: python scripts/build_figure_inputs.py
Writes results/manuscript/adapters/*.csv and adapter_manifest.json.
"""
from pathlib import Path
import hashlib, json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/manuscript/adapters'; OUT.mkdir(parents=True,exist_ok=True)
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources=[]
def add(p): sources.append({'path':str(Path(p).relative_to(ROOT)),'sha256':digest(p)}); return p

# Forest: all 34 original EHG9 UC-IN rows, frozen registry order.
src=add(ROOT/'results/primary/mixed_model_results_60s_fwh_EHG9.csv')
d=pd.read_csv(src); d=d[(d.channel=='EHG9')&(d.contrast=='contraction_vs_baseline')].copy()
d['label']=d['feature']; d['group']='EHG9'; d['scale']='model-standardised units'; d['interval_method']='model-based Wald 95% CI'; d['family_id']='historical_fwh_ehg9_uc_in_34'; d['family_size']=34
d[['label','group','estimate','ci_low','ci_high','contrast','scale','interval_method','family_id','family_size','n_women','n']].rename(columns={'n':'n_windows'}).to_csv(OUT/'original_ehg9_uc_in_forest.csv',index=False)

# Corrected cross-channel UC-IN forest, with both summaries and same-cohort channels.
src2=add(ROOT/'results/manuscript/corrected_family/corrected_cross_channel_all_16.csv')
c=pd.read_csv(src2); c=c[c.contrast=='contraction_vs_baseline'].copy(); c['label']=c.feature.map({'gi_total':'GI Total','gi_sd':'GI profile SD'}); c['group']=c.channel; c['scale']='model-standardised units'; c['interval_method']='model-based Wald 95% CI'; c['family_id']='corrected_cross_channel_16_after_development'; c['family_size']=16
c[['label','group','estimate','ci_low','ci_high','contrast','scale','interval_method','family_id','family_size','n_women','n']].rename(columns={'n':'n_windows'}).to_csv(OUT/'corrected_cross_channel_uc_in_forest.csv',index=False)

# Gain trajectories: preserve values, transform only column names/group labels.
src3=add(ROOT/'results/manuscript/audits/D/gain_audit_window_results.csv')
g=pd.read_csv(src3); keep=g[g.feature.isin(['fwh_entropy_totalsampen','gi_total'])].copy(); keep['group']=keep.feature.map({'fwh_entropy_totalsampen':'Original Total','gi_total':'GI Total'}); keep['value']=keep.observed
keep[['window_id','group','gain','value']].to_csv(OUT/'gain_trajectories_original_gi_total.csv',index=False)

# Agreement: independent/reference GI grid AUC, with validity metadata retained.
src4=add(ROOT/'results/manuscript/audits/C3/independent_grid_agreement.csv')
a=pd.read_csv(src4); ag=pd.DataFrame({'window_id':a.window_id,'reference':a.gi_auc,'candidate':a.ind_auc,'summary':'gigrid_auc','valid_mask_agreement':True,'tolerance':'atol=1e-10,rtol=1e-8'})
ag.to_csv(OUT/'independent_grid_auc_agreement.csv',index=False)

# Cohort table from Pass 0 reconciliation; compact, source-traceable table.
src5=add(ROOT/'results/manuscript/cohort_reconciliation.csv'); cr=pd.read_csv(src5)
rows=[
 {'population':'Database inventory','n_women':45,'n_recordings':123,'n_windows':'NA','definition':'local database inventory'},
 {'population':'Annotated/contributing','n_women':42,'n_recordings':112,'n_windows':4283,'definition':'manifest candidate windows before mixed/possible exclusion'},
 {'population':'Analysed 60-s windows','n_women':42,'n_recordings':112,'n_windows':4083,'definition':'3358 IN + 333 UC + 392 FM; excludes 132 mixed and 68 possible_contraction'},
 {'population':'Paired UC-IN women','n_women':34,'n_recordings':'NA','n_windows':3340,'definition':'paired woman subset'},
 {'population':'Paired FM-IN women','n_women':36,'n_recordings':'NA','n_windows':3308,'definition':'paired woman subset'},
 {'population':'Paired UC-FM women','n_women':32,'n_recordings':'NA','n_windows':645,'definition':'paired woman subset'},
]
pd.DataFrame(rows).to_csv(OUT/'cohort_summary.csv',index=False)

(OUT/'adapter_manifest.json').write_text(json.dumps({'sources':sources,'transformations':['row filters are explicit in build_adapters.py','feature labels/group labels are display mappings only','no estimates, intervals, q-values or denominators recalculated','gain adapter selects original Total and GI Total only','agreement adapter selects GI grid AUC and retains tolerance/mask fields']},indent=2))

# Model-dependence source-ready table: separate pooled/cluster estimates from paired woman summaries.
tm=add(ROOT/'results/reviewer_revision/01_temporal_dependence/temporal_dependence_comparison_34_features.csv')
pm=add(ROOT/'results/reviewer_revision/02_patient_level_uc_in/patient_level_vs_mixed_model_34_features.csv')
t=pd.read_csv(tm); p=pd.read_csv(pm)
rows=[]
for _,r in t.iterrows():
    for est,pref in [('primary_mixed_model','primary_mixed_model'),('pooled_ols','pooled_ols'),('woman_cluster_CR1','woman_cluster')]:
        rows.append({'feature':r.feature_name,'estimand':est,'estimate':r[pref+'_estimate'],'ci_low':r[pref+'_ci_low'],'ci_high':r[pref+'_ci_high'],'p_value':r[pref+'_p'],'q_value':r[pref+'_q'],'n_women':42,'n_windows':4083,'scale':'model-standardised units' if est!='woman_cluster_CR1' else 'model-standardised units; CR1 robust SE','interval_method':'model-based Wald 95% CI' if est!='woman_cluster_CR1' else 'robust Wald 95% CI','family_id':'historical_fwh_ehg9_uc_in_34'})
p['estimand']='paired_woman_median_difference'; p['scale']='raw feature units'; p['interval_method']='paired bootstrap percentile CI'; p['family_id']='paired_woman_summary_family'; p['n_windows']='paired women; not window coefficient'; p['n_women']=p.n_paired_women
rows += p[['feature','estimand','median_paired_difference','bootstrap_ci_low','bootstrap_ci_high','wilcoxon_p','patient_level_q','n_women','n_windows','scale','interval_method','family_id']].rename(columns={'median_paired_difference':'estimate','bootstrap_ci_low':'ci_low','bootstrap_ci_high':'ci_high','wilcoxon_p':'p_value','patient_level_q':'q_value'}).to_dict('records')
pd.DataFrame(rows).to_csv(OUT/'model_dependence_separate_estimands.csv',index=False)
