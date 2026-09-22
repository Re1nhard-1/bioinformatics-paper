"""Four-level absolute-error and paired-penalty profiles for every held-out donor."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scientific_style import COLORS,configure_style,figure_mm,save_figure_bundle
from run_deconvolution_pilot import sha
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('review',type=Path);a=p.parse_args()
    source=a.review.resolve();record=json.loads((source/'review.json').read_bytes())
    for x in record['artifacts']:assert sha(ROOT/x['path'])==x['sha256']
    donor=pd.read_csv(source/'donor_metrics.csv',dtype={'held_out':str})
    effects=pd.read_csv(source/'donor_penalties.csv',dtype={'held_out':str})
    levels=[20,30,40,50];labels=['20:20:20','30:15:15','40:10:10','50:5:5']
    donors=sorted(donor.held_out.unique(),key=int);markers=['o','s','^','v','D','P']
    font=configure_style();fig=figure_mm(183,165);axes=fig.subplots(2,2)
    fig.subplots_adjust(left=.11,right=.97,bottom=.18,top=.83,hspace=.58,wspace=.30)
    fig.text(.11,.965,'Reference allocation at a fixed budget of 60 cells per type',weight='bold')
    fig.text(.11,.933,'Grey: individual held-out donors; coloured: equal-donor mean')
    handles=[Line2D([],[],color='#8a8a8a',marker=m,mfc='white',lw=.6,ms=3,label=d) for d,m in zip(donors,markers)]
    fig.legend(handles=handles,ncol=6,loc='upper left',bbox_to_anchor=(.11,.9),columnspacing=1.7)
    data=[]
    absolute_bounds=[np.floor(donor.mae_pp.min()*2)/2-.2,np.ceil(donor.mae_pp.max()*2)/2+.2]
    penalty_bounds=[min(-.08,float(effects.mae_penalty_pp.min())-.08),max(.1,float(effects.mae_penalty_pp.max())+.08)]
    for col,(method,title,color,mean_marker) in enumerate([('music_weighted','Weighted MuSiC',COLORS['vermillion'],'s'),('music_nnls','Ordinary NNLS',COLORS['blue'],'o')]):
        values=donor[donor.method.eq(method)].pivot(index='held_out',columns='level',values='mae_pp').loc[donors,levels]
        penalties=values.sub(values[20],axis=0)
        expected=effects[effects.method.eq(method)].pivot(index='held_out',columns='level',values='mae_penalty_pp').loc[donors,[30,40,50]]
        assert np.max(np.abs(penalties[[30,40,50]].to_numpy()-expected.to_numpy()))<1e-12
        for row,(frame,quantity) in enumerate([(values,'absolute_mae'),(penalties,'mae_minus_balanced')]):
            ax=axes[row,col];panel='abcd'[row*2+col]
            for d,marker in zip(donors,markers):
                y=frame.loc[d].to_numpy()
                ax.plot(range(4),y,color='#969696',marker=marker,mfc='white',mec='#777777',lw=.55,ms=3,mew=.5,zorder=1)
                for level,value in zip(levels,y):data.append({'panel':panel,'method':method,'held_out':d,'summary':'donor','quantity':quantity,'level':level,'value_pp':float(value)})
            mean=frame.mean(axis=0).to_numpy()
            ax.plot(range(4),mean,color=color,marker=mean_marker,lw=1.1,ms=4,mfc=color,zorder=3)
            for level,value in zip(levels,mean):data.append({'panel':panel,'method':method,'held_out':'all','summary':'equal_donor_mean','quantity':quantity,'level':level,'value_pp':float(value)})
            ax.set_xticks(range(4),labels);ax.set_xlim(-.15,3.15)
            ax.tick_params(axis='x',pad=7)
            ax.set_title(title,loc='left',pad=10)
            ax.text(-.18,1.065,panel,transform=ax.transAxes,fontsize=8,weight='bold')
            ax.set_ylim(*(absolute_bounds if row==0 else penalty_bounds))
            lo,hi=ax.get_ylim()
            ax.set_yticks([v for v in ax.get_yticks() if lo <= v <= hi])
            ax.set_ylabel(('Absolute MAE' if row==0 else 'MAE minus balanced')+'\n(percentage points)')
            if row==1:
                ax.axhline(0,color=COLORS['neutral'],lw=.6,ls='--',zorder=0)
                ax.set_xlabel('Cells per donor within each type',labelpad=7)
    fig.text(.11,.081,'All six donor profiles are shown; segments connect tested allocations only.')
    fig.text(.11,.048,'Shared reference donors; synthetic mixtures; three draws; no population error bars.')
    caption="""# Fixed-budget reference allocation gradient

**a,b**, Absolute mean absolute error (MAE) of the official weighted MuSiC and internal ordinary NNLS outputs. Both panels share the same vertical scale. **c,d**, Within-donor MAE change relative to the identical 20:20:20 reference; both panels share a common change scale and zero line. Each grey trace is one of six held-out donors, identified by marker shape in the legend. Coloured traces give the equal-donor mean (vermilion squares: weighted; blue circles: ordinary). All 112 displayed values, including overlapping zeros, are retained in the source CSV. Differences use unrounded values.

Each allocation has 60 cells per type from the same three reference donors. Unequal allocations rotate all three dominant donors, averaged within each of the ten reference triples and three saved random draws for each held-out donor. All 60 raw-count synthetic targets per held-out donor are unchanged. The original balanced reference is shared by dominant-donor comparisons and is not independent replicated baseline data. A donor-specific random cell prefix is reused as that donor's allotment changes; the complete reference sets are not nested across increasing dominance.

Lines connect four tested allocations; they are not fitted continuous curves, error bars or thresholds. Six source SLE donors and overlapping reference sets define the scope. Three conditional reference draws do not establish population uncertainty. No population P values or confidence intervals are shown. The intermediate official-MuSiC predictions were added after observing the original endpoints, so this is an exploratory extension. The separate 5:5:5 thinning experiment changes total budget and is not part of this figure. Full reference-arrangement and block-level reversals are reported in accompanying tables.
"""
    out=ROOT/'results/music_gradient_figures'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');stem=out/'fixed_budget_gradient'
    assert len(data)==112
    save_figure_bundle(fig,stem,pd.DataFrame(data),caption,inputs=[source/'review.json',source/'donor_metrics.csv',source/'donor_penalties.csv'],producer=Path(__file__),details={'title':'MuSiC fixed-budget donor allocation gradient','font':font,'plotted_values':112,'no_population_inference':True})
    plt.close(fig);print(stem.relative_to(ROOT).as_posix())
if __name__=='__main__':main()
