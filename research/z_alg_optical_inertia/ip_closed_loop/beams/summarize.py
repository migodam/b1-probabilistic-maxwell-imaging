"""Summarize immutable witness and distance checks without rerunning physics."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PAIR = HERE.parent/'pair'
STUDIES=('main_32','main_96','confirm_32','confirm_96','power_32','power_96')
LABELS=dict(known_pair_coherent='Direct coherent',positive_quadratic='Positive-state coherent',
            indefinite_quadratic='Indefinite coherent',uniform_angle='Uniform angle',
            random_angle='Random angle',fisher_angle='Fisher angle',
            unconditioned_source_eigenbeam='Source-only eigenbeam')
COLORS=dict(known_pair_coherent='#386079',positive_quadratic='#3173a5',
            indefinite_quadratic='#155056',uniform_angle='#b58552',
            random_angle='#a76659',fisher_angle='#806687',
            unconditioned_source_eigenbeam='#557a64')


def read(path):
    return json.loads(path.read_text())


def quick_file(name):
    folder=HERE/name
    if not folder.exists():
        return None
    return folder/('results_final_secant.json' if name.endswith('96') else 'results_final.json')


def pair_file(name):
    ext=dict(main_32='_witnessed',main_96='_final').get(name,'')
    return HERE/f'paired_replay_{name}{ext}.json'


def report_pair(name, out):
    root=PAIR/name
    original=read(root/'proposal.json')
    future=np.load(root/'future_offline.npz')
    row_reports={}
    for i,row in out['rows'].items():
        key=i+'__'
        policies=original['rows'][i]['policies']
        trueidx=int(future[key+'label'])
        truth=future[key+'clean']
        first=future[key+'born']
        noise=future[key+'error']
        policies_rows={}
        for j,policy in enumerate(policies):
            tests={}
            for mode in ('positive','indefinite'):
                check=row['outcomes'][mode]['observed_wrong_class_checks'][policy]
                bounds=row['outcomes'][mode]['certified'][policy]
                spread=bounds.get('lower_distance_over_radius',0.)
                result=sorted(c['state'] for c in check)==['EXCLUDED','WITNESSED']
                tests[mode]=dict(strict_certified_distance_over_radius=spread,
                                 actually_verified_unique_class=result,
                                 state=[c['state'] for c in check],truth_label=trueidx,
                                 true_label_retained=check[trueidx]['state']!='EXCLUDED')
            data_norm=np.linalg.norm(truth[:,j])
            policies_rows[policy]=dict(checks=tests,
                noise_over_clean_port_norm=float(np.linalg.norm(noise)/data_norm),
                born_relative_error=float(np.linalg.norm(truth[:,j]-first[:,j])/data_norm),
                old_and_common_new_witnessed=row['per_action_common_reading_witness'][policy]['both_classes_have_same_legal_old_and_new_reading'])
        records=original['rows'][i]
        row_reports[i]=dict(observation_protocol_id=f'pair/{name}/ipcl_beam_{int(i):03d}__s{out["side"]}__tx180_rx180_16x16__added1',
                            fluence=records['fluence'], wall_seconds=row['full_row_cost']['wall_seconds'],
                            cpu_seconds=row['full_row_cost']['cpu_seconds'],
                            new_noise_radius=records['new_radius'],
                            generator=records.get('physical_source',original['physical_source']),
                            policies=policies_rows)
    return row_reports


def report_separate(name,out):
    origin=read(HERE/name/'proposal.json')
    counts={}
    for j,r in out['rows'].items():
        for a,v in r['policies'].items():
            k=counts.setdefault(a,dict(actual_strict=0,observed_exact_ambiguity=0,guaranteed_separation=0))
            k['actual_strict']+=int(v['actually_exactly_one_witnessed_and_other_excluded'])
            k['observed_exact_ambiguity']+=int(v['ambiguity_lower_witness_statement'] is not None)
            k['guaranteed_separation']+=int(v['verified_prospective_full_box_gap_over_new_radius_lower']>0)
    all_modes={i:r['all_coherent_actions'] for i,r in out['rows'].items()}
    return dict(counts=counts,n_materials=len(out['rows']),
                number_impossible_under_entire_dictionary=sum(r.get('worst_case_separation_impossible',False) for r in all_modes.values()),
                total_online_design_wall_seconds=origin['wall_seconds'],
                wall_seconds=out['wall_seconds'],cpu_seconds=out['cpu_seconds'],
                necessary_fluence_rows=all_modes)


def pair_counts(rows):
    if not rows:
        return {}
    labels=next(iter(rows.values()))['policies'].keys()
    return {name:dict(actual=sum(r['policies'][name]['checks']['indefinite']['actually_verified_unique_class'] or
                                 r['policies'][name]['checks']['positive']['actually_verified_unique_class'] for r in rows.values()),
                      actual_negative_state_readout=sum(r['policies'][name]['checks']['indefinite']['actually_verified_unique_class'] for r in rows.values()),
                      actual_positive_state_readout=sum(r['policies'][name]['checks']['positive']['actually_verified_unique_class'] for r in rows.values()),
                      prospective=sum(any(v['strict_certified_distance_over_radius']>2 for v in r['policies'][name]['checks'].values()) for r in rows.values()),
                      exact_worst_case_counterexamples=sum(r['policies'][name]['old_and_common_new_witnessed'] for r in rows.values()))
            for name in labels}


def plots(power, native):
    if not power or not native:
        return
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42})
    effect=pair_counts(power)
    names=list(effect)
    xx=np.arange(len(names))
    fig,axes=plt.subplots(1,2,figsize=(13,5),width_ratios=(1.25,1),layout='constrained')
    ax=axes[0]
    width=.34
    ax.bar(xx-width/2,[effect[k]['actual'] for k in names],width=width,
           color=[COLORS[k] for k in names],label='Actual decisions (either frozen proof)')
    ax.bar(xx+width/2,[effect[k]['exact_worst_case_counterexamples'] for k in names],width=width,
           color='white',edgecolor=[COLORS[k] for k in names],linewidth=1.2,
           hatch='////',label='Exact common-reading counterexamples')
    ax.set_xticks(xx, labels=[LABELS[n] for n in names],rotation=35,ha='right')
    ax.set_ylim(0,14.6)
    ax.set_yticks(np.arange(0,13,2))
    ax.set_ylabel('Independent materials (out of 12)')
    ax.set_title('32² — equal-fluence remeasurement')
    ax.grid(axis='y',alpha=.22)
    ax.legend(loc='upper right',frameon=False,fontsize=9)
    ax=axes[1]
    xx=np.arange(len(native))
    styles=(('known_pair_coherent','positive','D'),('positive_quadratic','positive','o'),
            ('indefinite_quadratic','indefinite','o'),('uniform_angle','positive','s'),('fisher_angle','positive','v'))
    for shift,(name,cert,marker) in zip((-.20,-.10,0,.10,.20),styles,strict=True):
        kvals=[r['policies'][name]['checks'][cert]['strict_certified_distance_over_radius'] for r in native.values()]
        ax.plot(xx+shift,kvals,marker,ms=7,label=LABELS[name],color=COLORS[name])
    ax.axhline(2,color='#504b42',ls='--',lw=1.25,label='Two noise radii')
    ax.set_xticks(xx,labels=[f'Case {j}' for j in native])
    ax.set_ylim(0,max(2.7,ax.get_ylim()[1]))
    ax.set_ylabel('Certified class distance / new noise radius')
    ax.set_title('96² — prospective whole-class distance')
    ax.legend(loc='lower right',frameon=False,fontsize=9)
    ax.grid(axis='y',alpha=.22)
    for suffix in ('png','pdf'):
        fig.savefig(HERE/f'beam_closed_loop.{suffix}',dpi=175)
    plt.close(fig)



def native_maps():
    path=PAIR/'power_96'
    if not (path/'proposal.json').exists():
        return
    doc=read(path/'proposal.json')
    bank=np.load(path/'public_classes.npz')
    fig,axes=plt.subplots(len(doc['rows']),3,figsize=(8.4,5.7),layout='constrained')
    foreground_max=max(np.max(bank[k+'__centers']) for k in doc['rows'])
    centered_max=max(np.max(np.abs(np.diff(bank[k+'__centers'],axis=0))) for k in doc['rows'])
    first_image=None
    diff_image=None
    for row,(k,v) in enumerate(doc['rows'].items()):
        img=np.zeros((2,96**2))
        img[:,bank[k+'__cells']]=bank[k+'__centers']
        allowed_background=np.zeros(96**2)
        allowed_background[bank[k+'__cells'][bank[k+'__bg']]]=1.
        background=allowed_background.reshape(96,96)
        base=img.reshape(2,96,96)
        for col in range(2):
            first_image=axes[row,col].imshow(base[col],extent=(-1,1,-1,1),origin='lower',
                                         cmap='viridis',vmin=0,vmax=foreground_max)
            axes[row,col].contour(background,levels=[.5],origin='lower',extent=(-1,1,-1,1),
                                  colors=['#e2dbec'],linewidths=.8,linestyles='dotted')
        diff_image=axes[row,2].imshow(base[1]-base[0],extent=(-1,1,-1,1),origin='lower',
                                    cmap='PuOr',vmin=-centered_max,vmax=centered_max)
        axes[row,0].set_ylabel(f'Case {k} • y (m)')
        for c in range(3):
            axes[row,c].set_xticks((-1,0,1))
            axes[row,c].set_yticks((-1,0,1))
            if row==len(doc['rows'])-1:
                axes[row,c].set_xlabel('x (m)')
    for ax,title in zip(axes[0],('Class 0 centre','Class 1 centre','Contrast difference'),strict=True):
        ax.set_title(title)
    fig.colorbar(first_image,ax=axes[:,:2],shrink=.78,label='Real contrast')
    fig.colorbar(diff_image,ax=axes[:,2],shrink=.78,label='Contrast difference')
    for suffix in ('png','pdf'):
        fig.savefig(HERE/f'native_material_pairs.{suffix}',dpi=175)
    plt.close(fig)

def main():
    base={}
    paired={}
    for name in STUDIES:
        a=quick_file(name)
        if a is not None and a.exists():
            base[name]=report_separate(name,read(a))
        f=pair_file(name)
        if f.exists():
            paired[name]=report_pair(name,read(f))
    sparse={}
    for side in (32,96):
        a=read(HERE/f'angle_bound_power_{side}.json')['rows']
        b=read(HERE/f'sparse_ports_power_{side}.json')['rows']
        c=read(HERE/f'sparse_physical_power_{side}.json')['rows']
        sparse[str(side)]={i:dict(ordinary_ports_ambiguous=sum(v['midpoint_is_legal_for_both'] for v in a[i]['angles']),
               lower_support_need_by_physical_witness=next((int(k)+1 for k,r in reversed(c[i]['caps'].items()) if r['all_supports_guaranteed_ambiguous']),1),
               upper_support_by_frozen_forms={k:v['verified_sufficient_port_count'] for k,v in row['cases'].items()})
               for i,row in b.items()}
    extra=read(HERE/'sparse_secant_power_96_receiver.json')['rows']
    for i,row in extra.items():
        sparse['96'][i]['optimized_two_port_real_box_gap_over_noise_lower']=row['best']['whole_box_gap_over_noise_lower']
    content=dict(secant_arm_trials=base,quadratic_pair_rows=paired,
                 quadratic_counts={k:pair_counts(v) for k,v in paired.items()},
                 sparse_support_summary=sparse)
    (HERE/'beam_summary.json').write_text(json.dumps(content,indent=2,allow_nan=False)+'\n')
    with (HERE/'paired_results.csv').open('w',newline='') as out:
        w=csv.writer(out)
        w.writerow(('study','case','policy','fluence','noise_radius','certification_form',
                    'certified_dual_distance_over_radius','class_0','class_1','truth_label',
                    'true_label_retained','common_reading_counterexample','born_relative','noise_over_clean'))
        for name,records in paired.items():
            for i,row in records.items():
                for p,pol in row['policies'].items():
                    for tag,cert in pol['checks'].items():
                        w.writerow((name,i,p,row['fluence'],row['new_noise_radius'],tag,
                                    cert['strict_certified_distance_over_radius'],*cert['state'],
                                    cert['truth_label'],cert['true_label_retained'],
                                    pol['old_and_common_new_witnessed'],pol['born_relative_error'],
                                    pol['noise_over_clean_port_norm']))
    plots(paired.get('power_32'),paired.get('power_96'))
    native_maps()


if __name__=='__main__':main()
