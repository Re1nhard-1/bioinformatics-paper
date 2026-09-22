"""Show all external donor profiles, absolute errors and the frozen primary contrast."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scientific_style import configure_style,figure_mm,save_figure_bundle,sha256,COLORS

ROOT=Path(__file__).resolve().parents[1]
LEVELS=['balanced','ratio2','ratio4','ratio10'];LABELS=['1:1:1','2:1:1','4:1:1','10:1:1']
METHODS=['music_weighted','music_nnls'];NAMES={'music_weighted':'Weighted MuSiC','music_nnls':'Ordinary NNLS'}
COL={'music_weighted':COLORS['vermillion'],'music_nnls':COLORS['blue']};MARK={'music_weighted':'s','music_nnls':'o'}
def ticks_in_bounds(ax,axis='y',nbins=4):
    low,high=ax.get_ylim() if axis=='y' else ax.get_xlim();ticks=MaxNLocator(nbins=nbins).tick_values(low,high)
    ticks=ticks[(ticks>=low)&(ticks<=high)]
    (ax.set_yticks if axis=='y' else ax.set_xticks)(ticks)
def panel(ax,label):ax.text(-.17,1.08,label,transform=ax.transAxes,fontsize=8,fontweight='bold',va='top')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('review',type=Path);a=p.parse_args();rv=a.review.resolve();info=json.loads((rv/'review.json').read_bytes())
    assert info['status']=='passed complete external descriptive audit'
    for item in info['artifacts']:assert sha256(ROOT/item['path'])==item['sha256']
    run=json.loads((ROOT/info['run']).read_bytes());inp=json.loads((ROOT/run['input_manifest']).read_bytes());donors=inp['donors']
    donor=pd.read_csv(rv/'donor_metrics.csv');agg=pd.read_csv(rv/'aggregate_metrics.csv');pen=pd.read_csv(rv/'donor_penalties.csv');primary=pd.read_csv(rv/'primary_endpoint_all_donors.csv')
    assert len(donor)==336 and len(agg)==24 and len(primary)==14
    out=ROOT/'results/music_external_figures'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    style=configure_style(7);fig=figure_mm(183,134);axes=fig.subplots(2,3)
    fig.subplots_adjust(left=.085,right=.985,bottom=.13,top=.92,wspace=.24,hspace=.40)
    ymax=float(donor.mae_pp.max())*1.08;plotted=[]
    for row,method in enumerate(METHODS):
        for col,budget in enumerate([60,120,300]):
            ax=axes[row,col];d=donor[(donor.method==method)&(donor.budget==budget)]
            for who in donors:
                vals=d[d.held_out==who].set_index('level').loc[LEVELS,'mae_pp'].to_numpy()
                ax.plot(range(4),vals,color='#B5B5B5',linewidth=.5,alpha=.85,zorder=1)
                plotted.extend({'panel':chr(97+row*3+col),'held_out':who,'method':method,'budget':budget,'level':l,'mae_pp':float(v),'summary':'donor'} for l,v in zip(LEVELS,vals))
            means=agg[(agg.method==method)&(agg.budget==budget)].set_index('level').loc[LEVELS,'mae_pp'].to_numpy()
            assert np.max(abs(means-d.groupby('level').mae_pp.mean().loc[LEVELS].to_numpy()))<1e-10
            ax.plot(range(4),means,color=COL[method],marker=MARK[method],linewidth=1,markersize=3,zorder=3)
            plotted.extend({'panel':chr(97+row*3+col),'held_out':'mean','method':method,'budget':budget,'level':l,'mae_pp':float(v),'summary':'mean'} for l,v in zip(LEVELS,means))
            ax.set(xlim=(-.15,3.15),ylim=(0,ymax),xticks=range(4),xticklabels=LABELS,title=f'{NAMES[method]}\n{budget} reference cells / type')
            ticks_in_bounds(ax);panel(ax,chr(97+row*3+col))
            if col==0:ax.set_ylabel('Mean absolute error (pp)')
            if row==1:ax.set_xlabel('Donor allocation ratio')
    fig.text(.085,.025,'Grey lines: all 14 donors. Coloured line: equal-donor mean. Shared y scale; no population error bars.',fontsize=7)
    caption='''# External absolute-error profiles

All 14 selected normal/na donors are shown at all three reference budgets and four allocation ratios. Panels a–c show weighted MuSiC (vermillion squares); d–f show its ordinary NNLS output using the same reference basis (blue circles). Thin grey lines show individual held-out-donor means over 60 fixed synthetic targets, three reference-cell draws, four predefined reference triples and, at unequal allocations, three dominant-donor choices. Bold coloured lines are equal-weight means of those 14 donor values. All six panels share an absolute-MAE scale in percentage points (pp); allocation ratios are discrete tested conditions, not a fitted continuous relationship. Reference budget is total cells per type across three donors; whole references have six times this many cells. No population confidence intervals or significance tests are shown. The selected high-inventory donors, shared reference draws/folds and synthetic targets limit interpretation; these are not measured-bulk validation data.
'''
    save_figure_bundle(fig,out/'external_absolute_profiles',pd.DataFrame(plotted),caption,[rv/'review.json',rv/'donor_metrics.csv',rv/'aggregate_metrics.csv'],Path(__file__),{'title':'External reference allocation and absolute error','style':style,'donors':14,'manual_visual_review':'pending'})
    plt.close(fig)
    fig=figure_mm(183,111);axes=fig.subplots(1,3,gridspec_kw={'width_ratios':[1,1,1.2]});fig.subplots_adjust(left=.082,right=.975,bottom=.23,top=.87,wspace=.52)
    strong=pen[pen.level=='ratio10'];low=min(0,float(strong.mae_penalty_pp.min()));high=max(0,float(strong.mae_penalty_pp.max()));pad=max((high-low)*.1,.025)
    plotted=[]
    for j,method in enumerate(METHODS):
        ax=axes[j];d=strong[strong.method==method]
        for who in donors:
            vals=d[d.held_out==who].set_index('budget').loc[[60,120,300],'mae_penalty_pp'].to_numpy()
            ax.plot(range(3),vals,color='#B5B5B5',linewidth=.5,marker='.',markersize=2,alpha=.85)
            plotted.extend({'panel':chr(97+j),'held_out':who,'method':method,'budget':b,'value_pp':float(v),'quantity':'strong_imbalance_penalty','summary':'donor'} for b,v in zip([60,120,300],vals))
        means=d.groupby('budget').mae_penalty_pp.mean().loc[[60,120,300]].to_numpy()
        ax.plot(range(3),means,color=COL[method],marker=MARK[method],linewidth=1,markersize=3)
        plotted.extend({'panel':chr(97+j),'held_out':'mean','method':method,'budget':b,'value_pp':float(v),'quantity':'strong_imbalance_penalty','summary':'mean'} for b,v in zip([60,120,300],means))
        ax.axhline(0,color='black',linewidth=.55,linestyle=':');ax.set(xlim=(-.15,2.15),ylim=(low-pad,high+pad),xticks=range(3),xticklabels=['60','120','300'],xlabel='Reference cells / type',title=NAMES[method]);ticks_in_bounds(ax);panel(ax,chr(97+j))
        if j==0:ax.set_ylabel('10:1:1 minus balanced MAE (pp)')
    ax=axes[2];p=primary.set_index('held_out').loc[donors];values=p.mae_penalty_change_pp.to_numpy()
    independent=strong[strong.method=='music_weighted'].pivot(index='held_out',columns='budget',values='mae_penalty_pp').loc[donors]
    assert np.max(abs(values-(independent[300]-independent[60]).to_numpy()))<1e-10
    ax.scatter(values,np.arange(14),c=COL['music_weighted'],marker='s',s=12,zorder=3)
    mean=float(values.mean());ax.axvline(0,color='black',linewidth=.55,linestyle=':');ax.axvline(mean,color=COL['music_weighted'],linewidth=.8,linestyle='--')
    xmin=min(0,float(values.min()));xmax=max(0,float(values.max()));pad=max((xmax-xmin)*.15,.025)
    ax.set(xlim=(xmin-pad,xmax+pad),ylim=(13.7,-.7),yticks=range(14),yticklabels=donors,title='Primary budget contrast',xlabel='Weighted I(300) − I(60) (pp)');ticks_in_bounds(ax,'x',3)
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,labelsize=6);panel(ax,'c')
    plotted.extend({'panel':'c','held_out':who,'method':'music_weighted','budget':0,'value_pp':float(v),'quantity':'I300_minus_I60','summary':'donor'} for who,v in zip(donors,values))
    plotted.append({'panel':'c','held_out':'mean','method':'music_weighted','budget':0,'value_pp':mean,'quantity':'I300_minus_I60','summary':'mean'})
    fig.text(.082,.09,'a, b: all donor profiles and their mean; identical y scales. c: every donor; dashed line is the mean.',fontsize=7)
    fig.text(.082,.045,'I(B) = strongest imbalance penalty at budget B. Negative budget contrast indicates a smaller penalty at B = 300.',fontsize=7)
    caption='''# External strongest imbalance penalty and primary budget contrast

Panels a and b show MAE(10:1:1) minus MAE(balanced), I(B), for weighted MuSiC and ordinary NNLS, in percentage points (pp), at all three frozen reference budgets. Grey profiles show all 14 held-out donors; coloured profiles show their equal-weight mean. The paired balanced baseline shares target, reference triple, block and budget. Individual donor values average all 60 targets, three draws, four triples and three dominant choices. The budgets are drawn at equally spaced categorical positions; no continuous trend is fitted. Panel c shows the predeclared primary contrast, weighted I(300)-I(60), for every donor in the frozen ring order. Its coloured dashed vertical line is the donor mean; dotted black lines mark zero. Negative values mean reduced strongest-imbalance penalty at the larger budget. Absolute errors must be read alongside the companion figure. Shared source donors/reference cells/targets and a high-inventory normal/na inclusion condition prevent these computational replicates from constituting independent population samples. No P values or population uncertainty intervals are supplied; targets are synthetic cell-count sums, not measured bulk samples.
'''
    save_figure_bundle(fig,out/'external_budget_penalty',pd.DataFrame(plotted),caption,[rv/'review.json',rv/'donor_penalties.csv',rv/'primary_endpoint_all_donors.csv'],Path(__file__),{'title':'External imbalance sensitivity and budget contrast','style':style,'donors':14,'primary_mean_pp':mean,'manual_visual_review':'pending'})
    plt.close(fig);print('OUTPUT',out.relative_to(ROOT).as_posix())

if __name__=='__main__':main()
