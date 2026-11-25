"""Independent dimensional, dense-coordinate, singular and budget chart controls."""
import json
from pathlib import Path
import sys

from flint import acb_mat, ctx
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'src'))
from fixed_trace import adj, block, mat, positive_dimension
from rectangular_chart import AffineEqualities, evaluate_chart, pullback_border, rectangular_equality_charts


def contains_zero(a):
    return all(x.contains(0) for x in a.entries())


def vcat(a,b):
    return acb_mat(a.nrows()+b.nrows(),a.ncols(),a.entries()+b.entries())


def check_cover(cover,lower,upper):
    """Verify the binary partition by rebuilding each parent from stored cuts."""
    actual={x.path:(x.lower,x.upper) for x in cover.leaves}
    cuts={x['path']:x for x in cover.subdivisions}
    current={'root':(tuple(lower),tuple(upper))}
    while current:
        label,(lo,hi)=current.popitem()
        if label not in cuts:
            assert actual.pop(label)==(lo,hi)
            continue
        split=cuts.pop(label)
        i,mid=split['parameter'],split['cut']
        assert lo[i]<mid<hi[i]
        current[label+'L']=(lo,hi[:i]+(mid,)+hi[i+1:])
        current[label+'R']=(lo[:i]+(mid,)+lo[i+1:],hi)
    assert not cuts and not actual


def dense_control():
    A0=np.array([[1.,0.,-.5,.25],[0.,1.,.25,-.5],[0.,0.,1.,-1.]])
    A1=np.array([[0.,0.,.125,.25],[0.,0.,-.375,.125],[0.,0.,.25,.25]])
    B0=np.array([[.5,-.25],[.125,.75],[.5,-.5]])
    B1=np.array([[.125,.25],[-.125,.25],[.25,-.125]])
    # Both exact coordinate changes are dyadic and invertible.
    R=np.array([[1.,.25,0.],[.5,1.,0.],[.25,.5,1.]])
    H=np.array([[1.,.25,-.25,.125],[.5,1.,.25,0.],[.125,0.,1.,.25],[.25,0.,-.125,1.]])
    ident=AffineEqualities(A0,B0,(A1,),(B1,))
    dense=AffineEqualities(R@A0@H,R@B0,(R@A1@H,),(R@B1,))
    lo,hi=(-.125,),(.25,)
    def chart(system):
        return rectangular_equality_charts(system,lo,hi,pivot_size=2,
                                           preferred_pivots=[((0,1),(0,1))])
    i_cover,d_cover=chart(ident),chart(dense)
    for cover in (i_cover,d_cover):
        check_cover(cover,lo,hi)
        assert cover.status=='OUTWARD_COVER' and not cover.subdivisions
        assert cover.leaves[0].remaining_rows==(2,)
    try:
        evaluate_chart(dense,i_cover.leaves[0],[0.])
    except ValueError:
        pass
    else:
        raise AssertionError('A chart must bind the original affine coefficients.')
    Q=np.eye(6)*2
    Q[0,4]=Q[4,0]=.25
    Q[1,5]=Q[5,1]=-.125
    H2=np.eye(6)
    H2[:4,:4]=H
    Qdense=H2.T@Q@H2
    threshold=.25
    points=[-.125,0.,.125,.25]
    for x in points:
        A=mat(A0+x*A1)
        B=mat(B0+x*B1)
        Z=mat(np.eye(2))
        # Independently solve all three equations, with x_2 as a prescribed free row.
        fixed=mat([[.375,.625]])
        independent_columns=(0,1,3)
        numeric_rhs=B*Z-block(A,range(3),[2])*fixed
        direct=block(A,range(3),independent_columns).solve(numeric_rhs)
        xi=vcat(vcat(block(direct,[0,1],range(2)),fixed),block(direct,[2],range(2)))
        xd=mat(H).solve(xi)
        ci=evaluate_chart(ident,i_cover.leaves[0],[x])
        cd=evaluate_chart(dense,d_cover.leaves[0],[x])
        out=[]
        for whole,leaf,q in ((xi,ci,Q),(xd,cd,Qdense)):
            params=vcat(block(whole,leaf.free_columns,range(2)),Z)
            assert contains_zero(leaf.lift*params-vcat(whole,Z))
            assert contains_zero(leaf.consistency*params)
            bad=acb_mat(params)
            bad[0,0]+=1
            assert not contains_zero(leaf.consistency*bad)
            pull=pullback_border(q,leaf,threshold)
            raw=pullback_border(q,leaf)
            assert pull.n==2 and pull.consistency.nrows()==1
            assert contains_zero(pull.full-raw.full+mat(np.diag([0.,0.,threshold,threshold])))
            assert positive_dimension(pull.full)==[4,0]
            out.append(adj(params)*pull.full*params)
        assert contains_zero(out[0]-out[1])
        # Check the entire root enclosures also contain the rational parameter slice.
        for whole,root in ((xi,i_cover.leaves[0]),(xd,d_cover.leaves[0])):
            p=vcat(block(whole,root.free_columns,range(2)),Z)
            assert contains_zero(root.lift*p-vcat(whole,Z))
    return dict(points=points,source_columns=2,retained_nonpivot_equations=1,
                dense_vs_identity_direct_solution=True,inconsistent_free_coordinate_rejected=True,
                source_threshold_pullback_preserved=True,positive_dimension_agrees=True,changed_system_rejected=True,
                identity=i_cover.record(),dense=d_cover.record())


def subdivision_control():
    variable=AffineEqualities([[0.,1.]],[[1.]],([[1.,-1.]],),([[0.]],))
    adaptive=rectangular_equality_charts(variable,[0.],[1.],pivot_size=1)
    check_cover(adaptive,[0.],[1.])
    assert adaptive.status=='OUTWARD_COVER' and len(adaptive.leaves)==2
    assert len(adaptive.subdivisions[0]['failed_attempts'])==2
    assert {x.column_pivot for x in adaptive.leaves}=={(0,),(1,)}
    budgets={}
    for name,restriction in [('depth',dict(max_depth=0)),('leaves',dict(max_leaves=1)),('time',dict(seconds=0.))]:
        row=rectangular_equality_charts(variable,[0.],[1.],pivot_size=1,**restriction)
        check_cover(row,[0.],[1.])
        assert row.status=='UNRESOLVED' and len(row.leaves)==1 and row.leaves[0].lift is None
        budgets[name]=row.record()
    return dict(adaptive=adaptive.record(),budgets=budgets)


def singular_control():
    crossing=AffineEqualities([[0.]],[[1.]],([[1.]],),([[0.]],))
    result=rectangular_equality_charts(crossing,[-.5],[.5],max_depth=3)
    check_cover(result,[-.5],[.5])
    zero_leaves=[x for x in result.leaves if x.lower[0]<=0<=x.upper[0]]
    assert result.status=='UNRESOLVED' and zero_leaves
    assert all(x.status=='UNRESOLVED' and x.lift is None for x in zero_leaves)
    assert any(x.status=='OUTWARD_CHART' for x in result.leaves)
    exact_zero=rectangular_equality_charts(crossing,[0.],[0.])
    assert exact_zero.status=='UNRESOLVED' and exact_zero.leaves[0].lift is None
    # At xi=0 a nonzero source is impossible; a zero source admits nonzero states.
    assert not contains_zero(mat([[0.]])*mat([[1.]])-mat([[1.]]))
    assert contains_zero(mat([[0.]])*mat([[1.]])-mat([[0.]]))
    return dict(pencil='xi x = z',zero_parameter_nonzero_source_impossible=True,
                zero_parameter_zero_source_nonzero_state_feasible=True,
                no_exclusion_from_singularity=True,mixed_box=result.record(),isolated_zero=exact_zero.record())


def main():
    receipt=dict(status='PASS',precision_bits=ctx.prec,dense_identity=dense_control(),
                 subdivision=subdivision_control(),rank_drop=singular_control())
    path=HERE/'chart_check.json'
    path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print('PASS: dense/identity equivalence, discarded-row counterexample, split cover, rank drop, and all budgets.')


if __name__=='__main__':
    ctx.prec=128
    ctx.threads=1
    with threadpool_limits(1):
        main()
