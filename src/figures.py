from pathlib import Path
import pandas as pd, matplotlib.pyplot as plt, numpy as np
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'; O=R/'figures'; O.mkdir(parents=True,exist_ok=True); src=[]
def record(name,path,analysis,notes=''): src.append({'figure_file':name,'source_result_file':str(path.relative_to(ROOT)),'analysis':analysis,'notes':notes})
def forest():
 d=pd.read_csv(R/'primary/mixed_model_results_60s_fwh_EHG9.csv'); d=d[(d.contrast=='contraction_vs_baseline') & d.primary_validated.astype(bool)].sort_values('estimate'); fig,ax=plt.subplots(figsize=(7,5)); y=np.arange(len(d)); ax.errorbar(d.estimate,y,xerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o'); ax.axvline(0,color='k',lw=.8); ax.set_yticks(y); ax.set_yticklabels(d.feature); ax.set_xlabel('Standardized UC versus IN effect (95% CI)'); fig.tight_layout(); fig.savefig(O/'figure_2_primary_EHG9.pdf'); plt.close(fig); record('figure_2_primary_EHG9.pdf',R/'primary/mixed_model_results_60s_fwh_EHG9.csv','primary EHG9')
def cross():
 d=pd.read_csv(R/'cross_channel/corrected_model_results.csv'); d=d[d.selected_for_primary_inference.astype(bool)]; fig,ax=plt.subplots(figsize=(8,5));
 for f,g in d.groupby('feature'): ax.plot(g.channel,g.estimate,marker='o',label=g.feature_label.iloc[0])
 ax.axhline(0,color='k',lw=.7); ax.set_ylabel('UC versus IN standardized effect'); ax.set_xlabel('Channel'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(O/'figure_4_cross_channel.pdf'); plt.close(fig); record('figure_4_cross_channel.pdf',R/'cross_channel/corrected_model_results.csv','cross-channel')
def entropy():
 d=pd.read_csv(R/'entropy/scale_sensitivity_summary.csv'); d=d[(d.stratum=='overall')&(d.feature=='total_sampen')]; fig,ax=plt.subplots(figsize=(6,4)); ax.plot(d.scale_factor,d.median_relative_change_percent,marker='o'); ax.axhline(0,color='k',lw=.7); ax.set_xscale('log',base=2); ax.set_xlabel('Scale factor'); ax.set_ylabel('Median change in Total SampEn (%)'); fig.tight_layout(); fig.savefig(O/'figure_5_sampen_scale.pdf'); plt.close(fig); record('figure_5_sampen_scale.pdf',R/'entropy/scale_sensitivity_summary.csv','entropy scale')
def main():
 forest(); cross(); entropy(); pd.DataFrame(src).to_csv(O/'figure_sources.csv',index=False)
if __name__=='__main__': main()
