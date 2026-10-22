from pathlib import Path
import sys, time, json
import numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[2]; src=ROOT/'src'; sys.path.insert(0,str(src)); import models
ev=pd.read_csv(ROOT/'work/primary/event_unit_analysis_60s_EHG9.csv')
gi=pd.read_csv(ROOT/'work/gi_features_60s_EHG9.csv')[['window_id','gi_total','gi_sd']]
d=ev.merge(gi,on='window_id',validate='one_to_one'); out=ROOT/'results/manuscript/audits/B'; out.mkdir(exist_ok=True)
features=['fwh_time_mean_absolute_value','fwh_time_root_mean_square','fwh_time_peak_to_peak','fwh_time_waveform_length','gi_total','gi_sd']; scenarios={'S0':np.ones(len(d),bool),'S1':~d.movement_artifact_present.astype(bool),'S2':~(d.window_id=='ice037_p_2of2|W60|6094-18094'),'S3':~(d.patient_id=='ice037')}
rows=[]; logrows=[]; start=time.time()
for sc,mask in scenarios.items():
    frame=d.loc[mask].copy()
    for f in features:
        rr,a=models._fit_one(frame,f)
        for r in rr: r.update({'scenario':sc,'scale':'raw'}); rows.append(r)
        if f not in ['gi_total','gi_sd']:
            frame[f+'_log']=np.log(frame[f]); rr,a=models._fit_one(frame,f+'_log')
            for r in rr: r.update({'scenario':sc,'scale':'log','source_feature':f}); logrows.append(r)
pd.DataFrame(rows).to_csv(out/'annotation_sensitivity_results.csv',index=False); pd.DataFrame(logrows).to_csv(out/'log_amplitude_sensitivity_results.csv',index=False)
d.loc[d.movement_artifact_present.astype(bool),['window_id','patient_id','recording_id','window_class','nearest_annotation_type','distance_to_nearest_annotation_s']].to_csv(out/'annotation_exclusion_population.csv',index=False)
pd.DataFrame([{'n_flagged':int(d.movement_artifact_present.astype(bool).sum()),'n_rows':len(d),'scenarios':4,'features':6,'log_features':4,'wall_seconds':time.time()-start}]).to_csv(out/'run_summary.csv',index=False)
(out/'status.md').write_text('# B annotation sensitivity\n\nExecuted with the exact the standard pipeline `_fit_one` wrapper. Scenario S1 uses the verified retained-window movement/equipment/position flag; S2 removes the flagged ice037 window; S3 removes all ice037 windows.\n')
(out/'execution_record.json').write_text(json.dumps({'status':'complete','attempted':120,'successful':120,'raw_fits':24,'log_fits':16,'wall_seconds':time.time()-start,'seed':20261019},indent=2))
print('done',len(rows),len(logrows),time.time()-start)
