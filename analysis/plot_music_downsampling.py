"""Plot all held-out donor summaries from the locked thinning policy comparison."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scientific_style import COLORS, configure_style, figure_mm, save_figure_bundle
from run_deconvolution_pilot import sha

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review',type=Path)
    args=parser.parse_args()
    source=args.review.resolve()
    review=json.loads((source/'review.json').read_bytes())
    for item in review['artifacts']:assert sha(ROOT/item['path'])==item['sha256']
    donor=pd.read_csv(source/'donor_metrics.csv',dtype={'held_out':str})
    effects=pd.read_csv(source/'donor_policy_effects.csv',dtype={'held_out':str})
    donors=sorted(donor.held_out.unique(),key=int)
    for method in ['music_weighted','music_nnls']:
        p=donor[donor.method.eq(method)].pivot(index='held_out',columns='policy',values='mae_pp')
        check=effects[effects.method.eq(method)].set_index('held_out').mae_thin_minus_keep_pp
        assert np.max(np.abs((p.thin15-p.keep60)-check))<1e-12
    font=configure_style()
    fig=figure_mm(183,108)
    axes=fig.subplots(1,2)
    fig.subplots_adjust(left=.105,right=.965,bottom=.24,top=.72,wspace=.33)
    fig.text(.105,.955,'Keeping cells versus donor-balanced thinning',weight='bold')
    color=COLORS['vermillion']
    fig.legend(handles=[Line2D([],[],marker='s',linestyle='none',color=color,label='Keep 50:5:5 (60 cells/type)'),
                        Line2D([],[],marker='s',linestyle='none',color=color,markerfacecolor='white',label='Thin to 5:5:5 (15 cells/type)')],
               loc='upper left',bbox_to_anchor=(.105,.91))
    fig.legend(handles=[Line2D([],[],marker='s',linestyle='none',color=color,label='Official MuSiC: weighted'),
                        Line2D([],[],marker='o',linestyle='none',color=COLORS['blue'],label='Official ordinary NNLS')],
               loc='upper left',bbox_to_anchor=(.60,.91))
    plotted=[]
    a,b=axes
    for i,d in enumerate(donors):
        row=donor[donor.held_out.eq(d)&donor.method.eq('music_weighted')].set_index('policy')
        keep=float(row.loc['keep60','mae_pp']); thin=float(row.loc['thin15','mae_pp'])
        a.plot([keep,thin],[i,i],color=color,lw=.8)
        a.plot(keep,i,marker='s',color=color,linestyle='none',ms=4)
        a.plot(thin,i,marker='s',color=color,mfc='white',linestyle='none',ms=4)
        for policy,value in [('keep60',keep),('thin15',thin)]:
            plotted.append({'panel':'a','held_out':d,'method':'music_weighted','quantity':'absolute_mae','policy':policy,'value_pp':value})
        for method,marker,c,offset in [('music_weighted','s',color,-.13),('music_nnls','o',COLORS['blue'],.13)]:
            value=float(effects[effects.held_out.eq(d)&effects.method.eq(method)].mae_thin_minus_keep_pp.iloc[0])
            b.plot(value,i+offset,marker=marker,color=c,linestyle='none',ms=4)
            plotted.append({'panel':'b','held_out':d,'method':method,'quantity':'thin_minus_keep_mae','policy':'paired','value_pp':value})
    for ax in axes:
        ax.set_yticks(range(6),donors)
        ax.set_ylim(5.55,-.55)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y',length=0)
    b.set_yticklabels([])
    a.set_ylabel('Held-out donor')
    a.set_title('Weighted MuSiC absolute error',loc='left',pad=13)
    b.set_title('Error change after thinning',loc='left',pad=13)
    a.text(-.20,1.07,'a',transform=a.transAxes,fontsize=8,weight='bold')
    b.text(-.10,1.07,'b',transform=b.transAxes,fontsize=8,weight='bold')
    a.set_xlabel('Mean absolute error\n(percentage points)')
    b.set_xlabel('Thin minus keep MAE\n(percentage points)')
    a.set_xlim(5.3,8.1); a.set_xticks([5.5,6,6.5,7,7.5,8])
    b.axvline(0,color=COLORS['neutral'],ls='--',lw=.6,zorder=0)
    b.set_xlim(-.40,.75); b.set_xticks([-.3,0,.3,.6])
    fig.text(.105,.092,'Positive values in b indicate worse error after discarding 75% of reference cells.')
    fig.text(.105,.045,'Six SLE source donors; synthetic mixtures; three draws; descriptive paired results.')
    caption='''# Keeping all cells versus donor-balanced thinning

**a**, Official weighted MuSiC MAE for each of six held-out donors, retaining all 50:5:5 cells per type (filled vermilion squares) or using its nested 5:5:5 subset (open squares). Lines connect paired summaries, not uncertainty intervals. **b**, Thin-minus-keep MAE for weighted MuSiC (vermilion squares) and the ordinary NNLS output (blue circles). Positive values indicate larger error after thinning. All values are percentage points; donor order is numeric, not effect-selected. The ordinary output is shown as a secondary comparison; weighted MuSiC remains the prespecified primary endpoint.

The simulated inventory contains 60 cells per type from the same three donors. Thinning retains 15/type and discards75% of reference cells; changes in sample size, expression/library-size estimation and gene support are part of the total policy effect. MuSiC already uses donor-equal basis construction; this is not an experiment correcting unequal donor weights. The available inventory is simulated from a larger deposited cohort, not a measured clinical resource limit.

Each donor summary averages all ten reference triples, three existing random draws, 60 fixed raw-count targets and all three dominant-donor choices. The identical thinned reference is shared across the three parent allocations; these comparisons are dependent. The plot has no population error bars, P values or confidence intervals. The previously studied20:20:20 reference requires additional cells from the two sparse donor/type groups and is not an available preprocessing option; it is not plotted here. This is same-cohort exploratory follow-up, not independent replication. Full arrangement-level reversals and the unchanged original results are retained in the accompanying review tables.
'''
    out=ROOT/'results/music_downsampling_figures'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    stem=out/'donor_balanced_thinning'
    save_figure_bundle(fig,stem,pd.DataFrame(plotted),caption,
        inputs=[source/'review.json',source/'donor_metrics.csv',source/'donor_policy_effects.csv'],
        producer=Path(__file__),details={'title':'MuSiC donor-balanced thinning policy','font':font,'plotted_values':24,'no_population_inference':True})
    plt.close(fig)
    print(stem.relative_to(ROOT).as_posix())

if __name__=='__main__':main()
