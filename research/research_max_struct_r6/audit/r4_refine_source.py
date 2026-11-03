def refine(model,Y,m,maximum=1000):
    before=Counter(model.cost);history=[];scale=1e6;start=time.perf_counter()
    best=[m.copy(),model.physical_score(m,Y)]
    def objective(x):
        value,g=model.data_value_gradient(x,Y)
        history.append([float(value),float(la.norm(g,np.inf))])
        if value<best[1]:best[:]=[x.copy(),float(value)]
        return scale*value,scale*g
    result=minimize(objective,m,method='L-BFGS-B',jac=True,bounds=[(.02,20.)]*model.n,
        options=dict(maxiter=maximum,maxfun=4*maximum,ftol=1e-12,gtol=.01,maxls=30))
    return best[0],dict(stop=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),
        history=history,cost=dict(Counter(model.cost)-before),wall_seconds=time.perf_counter()-start)
