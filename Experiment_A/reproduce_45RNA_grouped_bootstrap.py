"""Reproduce grouped-bootstrap summaries from the two included matched CSVs.
"""
import numpy as np
import pandas as pd
from pathlib import Path
root = Path(__file__).resolve().parent
out = root / "reproduced_tables"
out.mkdir(exist_ok=True)

methods={'ViennaRNA':'p_vienna','LinearPartition-V':'p_linear','EternaFold':'p_eterna'}
scenarios = {
    'Keep 7r6q_6': '45RNA_matched_keep_7r6q_6.csv',
    'Keep 7u0h_6': '45RNA_matched_keep_7u0h_6.csv',
}
metric_names=['ECE_nt_pooled','Brier_nt_pooled','LogLoss_nt_pooled','ECE_equal_section_mixture','ECE_mean_section','Brier_mean_section','LogLoss_mean_section']
n_boot=5000; seed=42
rng=np.random.default_rng(seed)
summary=[]; paired=[]; sectionrows=[]
for scenario, fname in scenarios.items():
    d = pd.read_csv(root / fname, dtype={"supp_section": "string"})
    groups=sorted(d.supp_section.unique(),key=lambda s:[int(x) for x in s.split('.')]); G=len(groups)
    n_records=d.structure_id.nunique(); n_nucs=len(d)
    print('scenario',scenario,'RNA',n_records,'groups',G,'nt',n_nucs)
    assert n_records==45 and G==27
    # Frozen bootstrap by seed, shared within each scenario and predictor. Separate draws between scenarios via RNG.
    draws=rng.integers(0,G,size=(n_boot,G)); C=np.zeros((n_boot,G),dtype=np.int16)
    for k in range(n_boot): C[k]=np.bincount(draws[k],minlength=G)
    # aggregate group/metric in array
    n_group=np.zeros(G,dtype=np.float64)
    group_arrays={}
    for method,col in methods.items():
        all_nbin=np.zeros((G,10)); all_sump=np.zeros((G,10)); all_sumy=np.zeros((G,10))
        brier_sums=np.zeros(G); loss_sums=np.zeros(G)
        for gi,g in enumerate(groups):
            sub=d[d.supp_section==g]
            y=sub.reference_paired_label.to_numpy(float)
            p=sub[col].to_numpy(float)
            assert np.all(np.isfinite(p)) and np.min(p)>=0 and np.max(p)<=1
            bins=np.minimum(np.floor(p*10).astype(int),9)
            n_group[gi]=len(sub)
            all_nbin[gi]=np.bincount(bins,minlength=10)
            all_sump[gi]=np.bincount(bins,weights=p,minlength=10)
            all_sumy[gi]=np.bincount(bins,weights=y,minlength=10)
            brier_sums[gi]=np.sum((p-y)**2)
            p_safe=np.clip(p,1e-12,1-1e-12)
            loss_sums[gi]=np.sum(-y*np.log(p_safe)-(1-y)*np.log(1-p_safe))
            sectionrows.append({'scenario':scenario,'supp_section':g,'predictor':method,'nt':len(sub),'RNAs':sub.structure_id.nunique(),'ECE_section':np.sum(np.abs(all_sump[gi]-all_sumy[gi]))/n_group[gi], 'Brier_section':brier_sums[gi]/n_group[gi], 'LogLoss_section':loss_sums[gi]/n_group[gi]})
        # group aggregates to metrics for one resample
        def calc(c):
            m = np.asarray(c,dtype=float)
            N=m@n_group
            diff=m@(all_sump-all_sumy)
            ece=np.abs(diff).sum(axis=-1)/N
            brier=(m@brier_sums)/N
            loss=(m@loss_sums)/N
            # equal-group mixture: normalize bins by section length, aggregate signed deviations, abs only after pooling bins
            ebins= m@((all_sump-all_sumy)/n_group[:,None])
            equal_ece=np.sum(np.abs(ebins),axis=-1)/m.sum(axis=-1)
            wg_ece=m@(np.sum(np.abs(all_sump-all_sumy),axis=1)/n_group)/m.sum(axis=-1)
            group_brier=m@(brier_sums/n_group)/m.sum(axis=-1)
            group_loss=m@(loss_sums/n_group)/m.sum(axis=-1)
            return np.stack([ece,brier,loss,equal_ece,wg_ece,group_brier,group_loss],axis=-1)
        observed=calc(np.ones(G)); boots=calc(C)
        group_arrays[method]=(observed,boots)
        for j,name in enumerate(metric_names):
            lo,hi=np.quantile(boots[:,j],[.025,.975])
            summary.append({'scenario':scenario,'predictor':method,'metric':name,'estimate':observed[j], 'CI_lower':lo,'CI_upper':hi,'n_RNAs':n_records,'n_nucleotides':n_nucs,'n_supp_sections':G, 'bootstrap_iterations':n_boot,'seed':seed})
    for m in ['ViennaRNA','LinearPartition-V']:
        obs=group_arrays[m][0]-group_arrays['EternaFold'][0]
        boot=group_arrays[m][1]-group_arrays['EternaFold'][1]
        for j,name in enumerate(metric_names):
            lo,hi=np.quantile(boot[:,j],[.025,.975])
            paired.append({'scenario':scenario,'contrast':m+' minus EternaFold','metric':name,'improvement_estimate':obs[j],'CI_lower':lo,'CI_upper':hi,'n_RNAs':n_records,'n_supp_sections':G, 'bootstrap_iterations':n_boot,'seed':seed})
        print('contrast',m,'pooled ECE improvement',obs[0],'95%ci',np.quantile(boot[:,0],[.025,.975]),'equal mix',obs[3],np.quantile(boot[:,3],[.025,.975]),'within group',obs[4],np.quantile(boot[:,4],[.025,.975]))
pd.DataFrame(summary).to_csv(out/'45RNA_group_bootstrap_per_predictor.csv',index=False)
pd.DataFrame(paired).to_csv(out/'45RNA_group_bootstrap_paired_differences.csv',index=False)
pd.DataFrame(sectionrows).to_csv(out/'45RNA_each_supp_section_metrics.csv',index=False)
# Count group wins under within-group ECE for either scenario
sections=pd.DataFrame(sectionrows)
for scen in scenarios:
 x=sections.query('scenario == @scen').pivot(index='supp_section',columns='predictor',values='ECE_section')
 print('scenario',scen,'Eterna group wins vs Vienna:',(x.EternaFold<x.ViennaRNA).sum(),'of',len(x),'ties', (x.EternaFold==x.ViennaRNA).sum())
print('saved outputs:',*[str(p) for p in out.iterdir()],sep='\n')
