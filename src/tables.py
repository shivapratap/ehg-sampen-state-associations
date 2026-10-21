from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'; O=R/'tables'; O.mkdir(parents=True,exist_ok=True)
P=ROOT/'results'/'primary'; W=ROOT/'work'; rows=[]
def save(name,df,sources,analysis,notes=''):
    df.to_csv(O/name,index=False)
    for s in sources: rows.append({'table_file':name,'source_result_file':str(s.relative_to(ROOT)),'analysis':analysis,'notes':notes})

def main():
    reg=pd.read_csv(W/'primary/confirmatory_feature_registry_60s_fwh.csv'); primary=pd.read_csv(P/'mixed_model_results_60s_fwh_EHG9.csv')
    pop=pd.DataFrame([{'quantity':'Women','value':42},{'quantity':'Recordings','value':112},{'quantity':'Candidate 60-s windows','value':4283},{'quantity':'Primary analysis units','value':4083},{'quantity':'UC units','value':333},{'quantity':'FM units','value':392},{'quantity':'IN units','value':3358}]); save('table_A_population.csv',pop,[W/'primary/confirmatory_feature_registry_60s_fwh.csv'],'population','Validated canonical counts')
    pri=primary[primary.contrast.eq('contraction_vs_baseline')].copy(); pri['beta']=pri['estimate']; pri['ci95']=pri.apply(lambda x:f"{x.ci_low:.3f} to {x.ci_high:.3f}",axis=1); mainf=pri[pri.primary_validated.astype(bool)][['feature','beta','ci95','p_value','q_value','primary_validated']]; save('table_B_primary_EHG9.csv',mainf,[P/'mixed_model_results_60s_fwh_EHG9.csv'],'primary EHG9','Supported features')
    save('table_S1_primary_EHG9_full_34.csv',pri,[P/'mixed_model_results_60s_fwh_EHG9.csv'],'primary EHG9','Full 34-feature table')
    cc=pd.read_csv(R/'cross_channel/corrected_cross_channel_comparison.csv'); save('table_D_cross_channel.csv',cc,[R/'cross_channel/corrected_cross_channel_comparison.csv'],'cross-channel','Phase E1 corrected results')
    ent=pd.read_csv(R/'entropy/corrected_amplitude_adjusted_model_results.csv'); save('table_E_entropy_amplitude.csv',ent,[R/'entropy/corrected_amplitude_adjusted_model_results.csv',R/'entropy/scale_sensitivity_summary.csv'],'entropy/amplitude','Phase E2 results')
    paired=pd.read_csv(R/'fm/paired_uc_vs_fm_34_features.csv'); save('table_F_fm_secondary.csv',paired,[R/'fm/paired_uc_vs_fm_34_features.csv',R/'fm/fm_vs_in_34_features.csv',R/'fm/clean_fm_sensitivity_summary.csv'],'FM secondary','Paired and clean-FM evidence')
    save('table_S2_cross_channel_detail.csv',pd.read_csv(R/'cross_channel/corrected_model_results.csv'),[R/'cross_channel/corrected_model_results.csv'],'cross-channel','Detailed effects')
    save('table_S3_entropy_scale.csv',pd.read_csv(R/'entropy/scale_sensitivity_summary.csv'),[R/'entropy/scale_sensitivity_summary.csv'],'entropy scale','Scale sensitivity')
    save('table_S4_rms_adjusted.csv',ent,[R/'entropy/corrected_amplitude_adjusted_model_results.csv'],'RMS adjustment','Detailed models')
    save('table_S5_fm_paired.csv',paired,[R/'fm/paired_uc_vs_fm_34_features.csv'],'FM paired','Detailed paired results')
    save('table_sources.csv',pd.DataFrame(rows),[], 'traceability')
if __name__=='__main__': main()
