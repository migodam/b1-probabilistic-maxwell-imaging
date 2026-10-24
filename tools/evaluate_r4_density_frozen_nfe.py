"""Fresh evaluation with development-frozen NFE32, original evaluator unchanged."""
import evaluate_r4_density as engine

def variants(frozen):
    rows=[]
    for method in frozen['methods']:
        if method['kind']=='cfm':
            if method.get('sample_nfes')!=[32]:raise ValueError('fresh CFM NFE32 must be frozen before observations')
            rows.append(dict(id=method['id']+'__nfe32',family=method['family']+'__nfe32',seed=method['seed'],nfe=32,method=method))
        else:
            if 'sample_nfes' in method:raise ValueError('NFE belongs only to CFM')
            rows.append(dict(id=method['id'],family=method['family'],seed=method['seed'],nfe=None,method=method))
    return rows

engine.variants=variants
if __name__=='__main__':engine.main()
