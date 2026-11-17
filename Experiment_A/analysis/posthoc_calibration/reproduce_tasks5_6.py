"""Project 7, Experiment A: published-section grouping audit and out-of-section temperature scaling.
Usage: python reproduce_tasks5_6.py --source /path/to/ExperimentA_50RNA_results.zip --output output_directory

"""
from __future__ import annotations
import argparse
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.optimize import minimize_scalar
from sklearn.isotonic import IsotonicRegression

EXCLUDE_COMMON = {'6t2c_V', '3egz_B', '3p22_E', '8hzj_A'}
SCENARIOS = {'keep_7r6q_6': EXCLUDE_COMMON | {'7u0h_6'}, 'keep_7u0h_6': EXCLUDE_COMMON | {'7r6q_6'}}
PREDICTORS = {'ViennaRNA': 'p_vienna', 'LinearPartition-V': 'p_linear', 'EternaFold': 'p_eterna'}
CLIP = 1e-9

def load_data(source):
    source = Path(source)
    if source.suffix.lower()=='.zip':
        with ZipFile(source) as z:
            members=[name for name in z.namelist() if name.endswith('ExperimentA_50RNA_matched.csv')]
            if not members: raise RuntimeError('ZIP does not contain ExperimentA_50RNA_matched.csv')
            with z.open(members[0]) as f:
                df=pd.read_csv(f, dtype={'supp_section':'string'})
    else: df=pd.read_csv(source, dtype={'supp_section':'string'})
    return df

def ece(p,y):
    b=np.minimum((np.asarray(p)*10).astype(int),9)
    y=np.asarray(y,dtype=float); p=np.asarray(p,dtype=float)
    if np.min(b)<0: raise ValueError('Probability <0')
    return sum(float((b==k).mean()*abs(p[b==k].mean()-y[b==k].mean())) for k in range(10) if (b==k).any())

def equal_group_mixture_ece(p,y,groups):
    # For each section, weight its entire nucleotide distribution by 1/G
    p=np.asarray(p,float); y=np.asarray(y,int); groups=np.asarray(groups)
    parts=[]
    for g in np.unique(groups):
        ix=groups==g; pp=p[ix]; yy=y[ix]; bb=np.minimum((pp*10).astype(int),9)
        parts.append([np.sum(pp[bb==k]-yy[bb==k])/len(pp) for k in range(10)])
    return float(np.sum(np.abs(np.mean(parts,axis=0))))

def score(p,y,groups):
    p=np.asarray(p,float); y=np.asarray(y,float); groups=np.asarray(groups)
    p_safe=np.clip(p,1e-12,1-1e-12)
    per_nt_loss=-(y*np.log(p_safe)+(1-y)*np.log1p(-p_safe))
    per_nt_brier=(p-y)**2
    return dict(ECE_pooled=ece(p,y), ECE_equal_section_mixture=equal_group_mixture_ece(p,y,groups),
                Brier_pooled=float(np.mean(per_nt_brier)), LogLoss_pooled=float(np.mean(per_nt_loss)),
                Brier_mean_section=float(np.mean([per_nt_brier[groups==g].mean() for g in np.unique(groups)])),
                LogLoss_mean_section=float(np.mean([per_nt_loss[groups==g].mean() for g in np.unique(groups)])))

def temp_fit(p,y):
    # Fit 1-parameter Bernoulli temperature using training RNA sections only, never test labels.
    z=logit(np.clip(np.asarray(p,float), CLIP,1-CLIP)); y=np.asarray(y,float)
    def objective(logT):
        pr=np.clip(expit(z/np.exp(logT)),1e-12,1-1e-12)
        return float(-np.mean(y*np.log(pr)+(1-y)*np.log1p(-pr)))
    res=minimize_scalar(objective,method='bounded',bounds=(np.log(.1),np.log(10)), options={'xatol':1e-6})
    return float(np.exp(res.x))

def evaluate(df, output_dir, section_lookup=None):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    if not {'structure_id','supp_section','reference_paired_label', *PREDICTORS.values()}.issubset(df.columns):
        raise ValueError('Required nucleotide-level columns are missing')
    if df.duplicated(['structure_id','position']).any(): raise ValueError('Duplicated nucleotide keys')
    if not df['reference_paired_label'].isin([0,1]).all(): raise ValueError('Invalid reference labels')
    if not np.isfinite(df[list(PREDICTORS.values())].to_numpy(dtype=float)).all(): raise ValueError('Nonfinite probabilities')
    if ((df[list(PREDICTORS.values())]<0) | (df[list(PREDICTORS.values())]>1)).any().any(): raise ValueError('Probabilities outside [0,1]')
    if df['supp_section'].isna().any(): raise ValueError('Missing supplementary section')
    lookup={}
    if section_lookup is not None:
        lookup=dict(zip(section_lookup.supp_section.astype(str), zip(section_lookup.supp_name,section_lookup.supp_category)))
    audit=[]; metrics=[]; predictions=[]; temps=[]
    for scenario, excluded in SCENARIOS.items():
        sub=df.loc[~df.structure_id.isin(excluded)].copy().reset_index(drop=True)
        if sub.structure_id.nunique()!=45 or sub.supp_section.nunique()!=27:
            raise ValueError(f'{scenario} expected 45 RNAs and 27 sections, got {sub.structure_id.nunique()} and {sub.supp_section.nunique()}')
        group_col=sub.supp_section.astype(str).to_numpy()
        y=sub.reference_paired_label.to_numpy(dtype=int)
        for (sid, grp), small in sub.groupby(['structure_id','supp_section']):
            name,category=lookup.get(str(grp),('See supplementary PDF','Unclassified'))
            audit.append(dict(scenario=scenario,structure_id=sid,supp_section=str(grp),section_name=name,RNA_type=category,
               n_nucleotides=len(small), published_47_group_id='NOT AVAILABLE',Rfam_family='NOT VERIFIED',
               family_disjoint_proven='NO', EternaFold_training_independence_proven='NO',
               interpretation='Supplement section proxy; not validated evolutionary family'))
        for method, pcol in PREDICTORS.items():
            p=sub[pcol].to_numpy(float)
            cv_temp=np.full(len(sub),np.nan); cv_iso=np.full(len(sub),np.nan)
            for grp in pd.unique(group_col):
                test=(group_col==grp); train=~test
                T=temp_fit(p[train],y[train]);cv_temp[test]=expit(logit(np.clip(p[test],CLIP,1-CLIP))/T)
                iso=IsotonicRegression(y_min=CLIP,y_max=1-CLIP,out_of_bounds='clip')
                iso.fit(p[train],y[train]);cv_iso[test]=iso.predict(p[test])
                temps.append(dict(scenario=scenario,predictor=method,heldout_section=grp,n_train_nt=int(train.sum()),n_test_nt=int(test.sum()),T=T))
            assert np.isfinite(cv_temp).all() and np.isfinite(cv_iso).all()
            for stage,q in [('Raw',p),('Temperature (out-of-section)',cv_temp),('Isotonic (out-of-section)',cv_iso)]:
                metrics.append(dict(scenario=scenario,predictor=method,correction=stage,n_RNAs=sub.structure_id.nunique(),
                    n_sections=sub.supp_section.nunique(),n_nt=len(sub),**score(q,y,group_col)))
            pred=sub[['structure_id','position','supp_section','reference_paired_label']].copy()
            pred['scenario']=scenario; pred['predictor']=method
            pred['p_raw']=p;pred['p_temperature_heldout']=cv_temp;pred['p_isotonic_heldout']=cv_iso
            predictions.append(pred)
    results=pd.DataFrame(metrics).sort_values(['scenario','predictor','correction'])
    results.to_csv(out/'02_out_of_section_calibration_metrics.csv',index=False)
    pd.concat(predictions,ignore_index=True).to_csv(out/'03_out_of_section_predictions.csv',index=False)
    pd.DataFrame(temps).to_csv(out/'04_train_only_temperature_values.csv',index=False)
    pd.DataFrame(audit).to_csv(out/'01_RNA_group_audit_UNVERIFIED_families.csv',index=False)
    print('\nTASK 5: 45 RNAs, 27 published supplementary sections per scenario. RNA-family assignments remain unverified.')
    print('\nTASK 6: Hold out one supplementary section at a time; fit calibrator ONLY on the other 26 sections.')
    print(results[['scenario','predictor','correction','ECE_pooled','ECE_equal_section_mixture','Brier_pooled','LogLoss_pooled']].round(4).to_string(index=False))
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);ap.add_argument('--crosswalk',default='')
    args=ap.parse_args();lookup=None
    if args.crosswalk:
        lookup=pd.read_csv(args.crosswalk,dtype={'supp_section':'string'}).drop_duplicates('supp_section')
    evaluate(load_data(args.source),args.output,lookup)
