"""Pure offline/material-action diagnostics. No Book/Target/backend import."""
import os,json
from pathlib import Path
import numpy as np
from .evaluation import checkpoint_identity


def check_representative_gates(physical,prepared_path,online_gate_path,receipt):
    physical=Path(physical);prepared_path=Path(prepared_path);receipt=Path(receipt)
    results=json.loads((physical/'RESULTS.json').read_text());closed=json.loads(receipt.read_text());started=json.loads((receipt.parent/'START.json').read_text())
    prepared=json.loads(prepared_path.read_text());snapshot=json.loads((physical/'FROZEN_REQUEST.json').read_text());online=json.loads(Path(online_gate_path).read_text())
    if results.get('status')!='COMPLETE' or results.get('GT_read') is not False or closed.get('status')!='COMPLETE' or started.get('kind')!='physics':raise ValueError('closedCOMPLETE actual representative physics Book required')
    if started.get('job')!=physical.name or receipt.parent.name!=physical.name or started.get('stage')!='II':raise ValueError('receipt must belong to this actual representative job')
    if started['pid']==os.getpid():raise ValueError('offline evaluator must be a separate process')
    if prepared!=snapshot or prepared.get('status')!='FROZEN_READY' or prepared.get('GT_read') is not False or not prepared['request'].get('frozen_estimator') or prepared['request'].get('no_GT') is not True:raise ValueError('actual frozen request mismatch')
    if online.get('status')!='ONLINE_COMPLETE' or online.get('truth_read') is not False or 'completed_time_ns' not in online:raise ValueError('actual observation-only inference gate required')
    if online['pid']==os.getpid():raise ValueError('offline evaluation must differ from inference process')
    source=prepared['source_predictions'];p=Path(source['path']);st=p.stat()
    if (st.st_size,st.st_mtime_ns)!=(source['size'],source['mtime_ns']):raise ValueError('frozen estimator arrays changed')
    # Require direct provenance to an actual inferred array, not an unsealed
    # assertion frozen_estimator=true. Field/generator gates explicitly seal;
    # ordinary risk gate's canonical estimator_id.npz is accepted with its gate.
    identity=checkpoint_identity(p);sealed=[]
    sealed.extend(online.get('field_artifact_identities',[]))
    sealed.extend([r['identity'] for r in online.get('artifacts',[])])
    canonical=Path(online_gate_path).parent/(prepared['request']['estimator_id']+'.npz')
    if identity not in sealed and p.resolve()!=canonical.resolve():raise ValueError('MISSING_ACTUAL_INFERENCE_ARTIFACT_LINK; no GT opened')
    freeze=online.get('freeze',{});inference_freeze=freeze.get('freeze_time_ns')
    if inference_freeze is None or not inference_freeze<online['completed_time_ns']:raise ValueError('inference method freeze must precede inference')
    preparation_time=prepared_path.stat().st_mtime_ns;start_ns=int(started['start_unix']*1e9)
    if st.st_mtime_ns>online['completed_time_ns']:raise ValueError('prediction arrays postdate online closure')
    if not start_ns<=(physical/'RESULTS.json').stat().st_mtime_ns<=receipt.stat().st_mtime_ns:raise ValueError('actual run result must precede its Book closure')
    if not online['completed_time_ns']<=preparation_time<=start_ns:raise ValueError('inference -> requestfreeze -> physicalstart chronology required')
    if (physical/'FROZEN_REQUEST.json').stat().st_mtime_ns<preparation_time or receipt.stat().st_mtime_ns<start_ns:raise ValueError('request/closure chronology failed')
    if set(prepared['request']['parents'])-set(online.get('parents',freeze.get('parents',[]))):raise ValueError('declared request parents absent from inference gate')
    return prepared,results,dict(physical_pid=started['pid'],inference_pid=online['pid'],receipt=str(receipt.resolve()),
                                 inference_gate=str(Path(online_gate_path).resolve()),request_freeze=str(prepared_path.resolve()),chronology='VERIFIED_BEFORE_GT')


def action_properties(chi,config,roi,tolerance=1e-10):
    chi=np.asarray(chi);roi=np.asarray(roi,int)
    if chi.shape!=(576,) or not np.iscomplexobj(chi) or not np.isfinite(chi).all():raise ValueError('finite arbitrary complex material576 required')
    host=complex(*config['host_chi']);outside=np.ones(576,bool);outside[roi]=False;rmax=config['inclusion_increment_real_interval'][1];imax=config['inclusion_increment_imag_interval'][1]
    bounds=(chi[roi].real>=host.real-tolerance)&(chi[roi].real<=host.real+rmax+tolerance)&(chi[roi].imag>=host.imag-tolerance)&(chi[roi].imag<=host.imag+imax+tolerance)
    return dict(passive_imag_nonnegative=bool(np.all(chi.imag>=-tolerance)),ROI_outside_host_fixed=bool(np.all(abs(chi[outside]-host)<=tolerance)),
        original_component_rectangle=bool(np.all(bounds)),rectangle_violating_cells=int((~bounds).sum()),
        decoder_realizability='NOT_PROVEN',P1_legal_world='NOT_CLAIMED',role='arbitrary material risk action; Full response does not establish prior realizability')
