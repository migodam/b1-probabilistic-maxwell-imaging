"""Synthetic accounting/decoder/isolation checks. No physics solves."""
import tempfile,json,time
from pathlib import Path
import numpy as np
from research.p3l_r5_parallel_ab.shared import physics_adapter as m

def run():
    tests={};original=m.SHARED
    target=m.Target(job='decoder_qa')
    x=np.random.default_rng(1106).normal(size=66);v=np.random.default_rng(1107).normal(size=66)
    D=target.decoder_jacobian(x);eps=1e-6
    finite=(target.decoder(x+eps*v)[0]-target.decoder(x-eps*v)[0])/(2*eps)
    err=float(np.linalg.norm(finite-D@v)/np.linalg.norm(D@v))
    assert err<1e-7;tests['pure_decoder_directional_relative_error']=err
    z=np.arange(256).reshape(32,8)+1j*np.arange(256,512).reshape(32,8)
    assert np.array_equal(m.pack(z)[:64],np.r_[z[:,0].real,z[:,0].imag]);tests['source_major_pack']='PASS'
    auth=json.loads((original/'cost_ledger/AUTHORIZATION.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        m.SHARED=Path(tmp);(m.SHARED/'cost_ledger').mkdir();m.write(m.SHARED/'cost_ledger/AUTHORIZATION.json',auth)
        b=m.Ledger('A','test');b.kind='forward';b.reserve(states=1,forward_rhs=8)
        assert b.totals()['states']==1 and b.totals()['forward_rhs']==8
        assert m.Ledger('A','restart').totals()['forward_rhs']==8;tests['restart_counting']='PASS'
        try:b.reserve(states=32768)
        except m.BudgetStop:pass
        else:raise AssertionError('cap not enforced')
        assert b.totals()['states']==1;tests['cap_rejection_no_reset']='PASS'
        try:
            with b.compute('derivative','synthetic_failure'):
                b.reserve_rhs(8,'fake_no_solver')
                raise RuntimeError('synthetic failure')
        except RuntimeError:pass
        assert b.totals()['derivative_rhs']==8 and b.totals()['failures']==1;tests['failures_retained']='PASS'
        b.append(dict(event='COMPUTE_BEGIN',id='interrupted',kind='forward',purpose='test'))
        t=m.Ledger('A','restart2').totals();assert t['unfinished_actions']==1 and t['wall_seconds']>0;tests['interrupted_action_conservative_time']='PASS'
        b.append(dict(event='COMPUTE_END',id='interrupted',cpu_seconds=1.,wall_seconds=86401.,error='synthetic'))
        try:b._check(b.totals())
        except m.BudgetStop:pass
        else:raise AssertionError('wall cap not enforced')
        tests['controlled_wall_cap']='PASS'
        b.path.write_text(b.path.read_text()+'{"partial":')
        try:m.Ledger('A','partial').totals()
        except m.BudgetStop:pass
        else:raise AssertionError('partial tail ignored')
        tests['partial_tail_fails_closed']='PASS'
    m.SHARED=original
    try:target.load_observation('/tmp/oracle/forbidden.npz')
    except ValueError:pass
    else:raise AssertionError('oracle path accepted')
    tests['online_oracle_isolation']='PASS'
    return dict(status='PASS',physics_solves=0,tests=tests)
if __name__=='__main__':m.write(m.OUT/'audit/ACCOUNTING_TESTS.json',run())
