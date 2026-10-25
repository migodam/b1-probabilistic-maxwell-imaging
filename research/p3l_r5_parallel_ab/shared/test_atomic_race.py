"""Synthetic concurrent cap contention; no Maxwell backend import/solve."""
import json,tempfile,multiprocessing as mp
from pathlib import Path
from . import physics_adapter as m

def contender(shared,queue):
    m.SHARED=Path(shared)
    b=m.Ledger('A','synthetic_contender')
    try:b.reserve(states=1,forward_rhs=8);queue.put('reserved')
    except m.BudgetStop:queue.put('denied')

def run():
    auth=json.loads((m.SHARED/'cost_ledger/AUTHORIZATION.json').read_text())
    auth['A'].update(forward_states=1,forward_rhs=8,all_rhs=8)
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/'cost_ledger';p.mkdir();m.write(p/'AUTHORIZATION.json',auth)
        ctx=mp.get_context('spawn');q=ctx.Queue();children=[ctx.Process(target=contender,args=(tmp,q)) for _ in range(4)]
        for c in children:c.start()
        result=[q.get(timeout=20) for _ in children]
        for c in children:c.join(20);assert c.exitcode==0
        assert result.count('reserved')==1 and result.count('denied')==3,result
        lines=[json.loads(r) for r in (p/'LEDGER.jsonl').read_text().splitlines()]
        assert sum(r.get('states',0) for r in lines)==1
    m.write(m.OUT/'audit/ATOMIC_BUDGET_TEST.json',dict(status='PASS',contenders=4,reserved=1,denied=3,physics_RHS=0,synthetic_RHS_reservations=8,production_ledger_untouched=True))
if __name__=='__main__':run()
