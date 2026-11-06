from pathlib import Path
import json,os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'RESULTS/vector3d';FIG=ROOT/'FIGURES/vector3d'
FIG.mkdir(exist_ok=True,parents=True)
def load(p):return json.loads(p.read_text())
names=['four_cubes','hollow_cube_rods','lossy_coated_cube'];methods=['LSIE_sourcewise','CIE_sourcewise','Dual_sourcewise','CIE_global','A_global','Aplus_global']
fig,axes=plt.subplots(4,7,figsize=(17,9));titles=['truth','LSIE','CIE','dual β','global CIE','A','A+']
for row,(name,part) in enumerate([(names[0],'real'),(names[1],'real'),(names[2],'real'),(names[2],'imag')]):
    target=OUT/name;truth=np.load(target/'truth_for_scoring.npz')['chi'];arrays=[truth]
    for m in methods:
        p=target/(m+'_prediction.npz');arrays.append(np.load(p)['prediction'] if p.exists() else None)
    for col,a in enumerate(arrays):
        ax=axes[row,col]
        if a is None:ax.text(.5,.5,'Unavailable\ncomplex arc',ha='center',va='center');ax.set_axis_off();continue
        q=getattr(a,part).reshape((6,)*3);im=ax.imshow(q[:,:,3].T,origin='lower',vmin=0,vmax=4 if part=='real' else 1.2,cmap='magma',extent=[-.3,.3,-.3,.3]);ax.set_title(titles[col] if row==0 else '',fontsize=10);ax.set_xticks([]);ax.set_yticks([])
        if col==0:ax.set_ylabel(name+' '+part,fontsize=9)
fig.suptitle('Cross-grid 3D full-vector inversions: 48 PR cycles + 300 physical steps (orthonormal current coordinates)\nCentral z slice; scaled .6λ domain, 60 dipoles, 10% noise; material accuracy is not implied by residual',fontsize=12)
fig.tight_layout(rect=[0,0,.94,.94]);cb1=fig.add_axes([.955,.33,.009,.54]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,4),cmap='magma'),cax=cb1,label='Reχ');cb2=fig.add_axes([.955,.06,.009,.20]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,1.2),cmap='magma'),cax=cb2,label='Imχ');fig.savefig(FIG/'material_slices.png',dpi=170);plt.close(fig)

# Offline score sealed trajectories; never used for selecting an iterate.
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for kind in methods+['direct_physical']:
    p=OUT/'four_cubes'/(kind+'_prediction.npz');r=OUT/'four_cubes'/(kind+'_receipt.json')
    if not p.exists():continue
    a=np.load(p);truth=np.load(OUT/'four_cubes/truth_for_scoring.npz')['chi'];rr=load(r);nrmse=np.linalg.norm(a['trials']-truth,axis=1)/np.linalg.norm(truth);curve=rr['curve'];phi=np.array([q['phi'] for q in curve]);inc=np.array([np.argmin(phi[:j+1]) for j in range(len(phi))]);nrmse=nrmse[inc];axes[0].plot([q['physical_LU'] for q in curve],nrmse,label=kind);axes[1].plot([q['cpu_s'] for q in curve],nrmse,label=kind)
axes[0].set_xlabel('Cumulative physical LU (PR work excluded)');axes[1].set_xlabel('Cumulative solver CPU s (all online work)')
for ax in axes:ax.set_ylabel('Material NRMSE, scored only after seal');ax.grid(alpha=.2)
axes[1].legend(fontsize=7);fig.suptitle('Scaled four cubes: public-objective-selected paths; same300-step corrector');fig.tight_layout();fig.savefig(FIG/'four_cube_cost_curves.png',dpi=170);plt.close(fig)
fig=plt.figure(figsize=(14,7))
for i,(name,kind) in enumerate([(n,k) for n in names for k in ['truth','A_global']]):
    ax=fig.add_subplot(2,3,(i//2)+(i%2)*3+1,projection='3d');target=OUT/name;arr=np.load(target/'truth_for_scoring.npz');points=arr['points'];chi=arr['chi'] if kind=='truth' else np.load(target/(kind+'_prediction.npz'))['prediction'];mask=chi.real>.25;ax.scatter(*points[mask].T,c=chi.real[mask],vmin=0,vmax=4,cmap='magma',s=30,marker='s',alpha=.65);ax.set(xlim=(-.3,.3),ylim=(-.3,.3),zlim=(-.3,.3),title=name+' '+kind);ax.set_box_aspect([1,1,1]);ax.tick_params(labelsize=6)
fig.suptitle('Actual3D material voxel centers (display threshold Reχ>.25 only; inverse support unrestricted)');fig.tight_layout();fig.savefig(FIG/'material_voxels.png',dpi=160);plt.close(fig)

rows=load(OUT/'SCORES.json');pr=load(OUT/'PR_SCORES.json') if (OUT/'PR_SCORES.json').exists() else [];r240=load(OUT/'pr240/PR_SCORES.json') if (OUT/'pr240/PR_SCORES.json').exists() else []
complex_scores=[]
for kind in methods+['direct_physical']:
    target=OUT/'lossy_coated_cube';p=target/(kind+'_prediction.npz')
    if not p.exists():continue
    truth=np.load(target/'truth_for_scoring.npz')['chi'];a=np.load(p)['prediction'];complex_scores.append(dict(method=kind,real_chi_NRMSE=float(np.linalg.norm(a.real-truth.real)/np.linalg.norm(truth.real)),imag_chi_NRMSE=float(np.linalg.norm(a.imag-truth.imag)/np.linalg.norm(truth.imag))))
(OUT/'LOSSY_COMPONENT_SCORES.json').write_text(json.dumps(complex_scores,indent=2))
def table(rows):
    text='|Target|Method|Raw NRMSE|Final NRMSE|Residual|CPU s|Physical LU|Stop|\n|---|---|---:|---:|---:|---:|---:|---|\n'
    for r in rows:text+=f"|{r['target']}|{r['method']}|{r['raw_NRMSE']:.4f}|{r['final_NRMSE']:.4f}|{r.get('residual',0):.4f}|{r['cpu_s']:.2f}|{r['counts']['physical_LU']}|{r.get('stop',{}).get('message','see receipt')}|\n"
    return text
lines=['# R8 genuine 3D vector material inversion refinement','', '**Status: actual cross-grid multi-object inversions; no 3D A superiority established.** Material scores use exposed development targets, not blind holdouts. Every prediction and trajectory is written before scorer access. R1–R7 remain untouched.','', '## Original 2019 physical and algorithmic contract','', '[Zhong and Xu 2019, DOI10.3390/jimaging5020027](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/) uses a 3λ cuboid, 60 dipole illuminations and 60 colocated vector receivers on three offset circles of radius3λ, all three electric components (10800 complex observations), fine60³ synthetic/inverse30³ grids and10% additive noise. Four cubes have side.8λ, eps2/3/4/5 and centers (−.7,.7,−.7),(.7,−.7,−.7),(−.7,−.7,.7),(.7,.7,.7)λ. The coated cube outer side2λ has eps1.5+.3i, inner side1λ eps2+.8i. The hollow cube+two rods use eps2.5; exact dimensions are not provided in the article prose, so our declared dimensions are an adaptation. Paper CIE uses SVD DPCS, projected low-frequency FFT APCS, source-normalized objectives, beta6→1, and optional distinct current/material beta.','', '## Implemented 3D contract and adaptations','', 'The scaled set uses .6λ cube,6³ independent isotropic inverse material DOFs,8³ independently assembled synthetic physics,60 finite-distance electric dipoles and60 full-vector receivers with paper circle offsets scaled by.2,10% fixed-seed complex Gaussian noise. All object dimensions are scaled by.2 except the explicitly adapted hollow geometry. Full R1 dyadic near-field coupling, all three current/field components, CM+RR self and complex passive material convention exp(−iωt) are retained. Lossless chi is bounded[0,4]; lossy Rechi[0,4],Imchi[0,1.2]. Exact analytic box volume fractions stabilize declared geometry; volume-averaged chi is an interface approximation. This is smaller electrical size than the paper and is explicitly not a paper reproduction.','', 'Receiver SVD threshold=.1 of largest singular value is fixed from public operator/noise scale. APCS= (I−VV*) IDFT(mask γ); actual3D inverse FFT constructs low-frequency columns with MF1→2. The first controls use an SVD orthonormal range of these projected columns: this preserves the affine family but preconditions the PR coefficient metric and is not the original Fourier gamma coordinates. The final six-arm240-cycle Fourier-coordinate run retains the original projected IDFT columns without whitening. Dense dyadic propagation remains; no FFT-complexity claim. First15 controls solve each reduced quadratic current problem exactly; these are development controls. Final sourcewise/global Polak–Ribiere implementations use the same DPCS/APCS affine family and exact quadratic step length. Global shares one coefficient over all source columns; sourcewise uses per-source coefficients. The dual control uses current beta3/material beta1 in stage2.','', 'A uses self-stripped b=beta−s and the actual full3×3 off-diagonal dyadic H within disjoint2³ voxel blocks: u=beta(j−chi Q)/(1+bchi), R=u+zHu,z=chi/(1+bchi). A material updates include the complete material-dependent derivative. A+ uses real scalar-isotropic common self sI3: conditional R=C+D z_i across all three components and sources. State objective is q|z_i|²+2Re(h z_i); stationary numerator qx(1+Re(b)x)+Re[h(1+conj(b)x)²]. Two coordinate sweeps enumerate legal roots/endpoints/current, compare the complete original quadratic prior and retain the old state if the whole loss increases by more than the declared1e−13 numerical tolerance. This is a state-stationary proposal method, not exact regularized minimization. Complex chi spans a two-real-dimensional region; A+ is NOT_AVAILABLE for loss, while A uses bounded joint Re/Im optimization.','', 'All candidates are selected by original full-Maxwell data loss plus the same smooth/l2 prior; truth never enters selection. Every arm receives max300 full-physics L-BFGS-B correction. Cap exits are not convergence. Direct control has192 initial physical steps plus same300 allowance. Equal iteration counts are not equal compute; LU/RHS, current factors/PR steps, CPU and full trajectories are preserved.','', '## Cross-grid discrepancy and physics qualification','']
for r in load(OUT/'SIMULATION.json')['rows']:lines.append(f"- {r['target']}: independently assembled8³→6³ forward discrepancy {100*r['crossgrid_forward_discrepancy']:.2f}%; additive noise10%.")
lines+=['','Complex off-root physical derivative relative error3.83e−9; LSIE/CIE/A errors≤5.62e−10. This necessary new complex adapter check is separate from frozen identities. Same LU with multiple sources is counted as one factorization with60 RHS, never independent validation. Cross-grid simulations use distinct meshes, assemblies and factorizations; they are distinct discretizations of the same DDA model, not independent continuum solvers.']
if (OUT/'mie_geometry/RECEIPT.json').exists():
    rr=load(OUT/'mie_geometry/RECEIPT.json');lines+=['','One new physically motivated geometry correction addresses R7 boundary bias on its .4λ box (h=.08/.06667/.05λ for5³/6³/8³). These are not the R8 .6λ imaging grids (h=.1/.075λ), so the Mie percentages cannot be assigned to those reconstructions. The correction is: use electrostatic polarizability-preserving Maxwell–Garnett mixing chi_eff=3 f chi/[3+(1−f)chi] with fixed24³ midpoint cell fractions, no fitted parameter, volume rescaling or scan. This is not exact cut-cell interface physics.']
    for r in rr['rows']:lines.append(f"- Sphere radius.145λ eps3: grid{r['nside']}³, h={r['h']:.5f}λ, Mie far-vector error {100*r['mie_far_vector_relative_error']:.3f}%, volume {r['effective_volume']:.7f} versus {r['exact_volume']:.7f}.")
lines+=['','## Actual scaled multi-object results','',table(pr),'','![Actual3D material slices](FIGURES/vector3d/material_slices.png)','', '![Actual3D voxel centers](FIGURES/vector3d/material_voxels.png)','', '![Offline material-cost curves](FIGURES/vector3d/four_cube_cost_curves.png)','', '## Lossy real and imaginary material scores','','|Method|Real χ NRMSE|Imaginary χ NRMSE|\n|---|---:|---:|\n'+''.join(f"|{q['method']}|{q['real_chi_NRMSE']:.4f}|{q['imag_chi_NRMSE']:.4f}|\n" for q in complex_scores),'','## Preserved exact-current adapted controls','',table(rows)]
if (OUT/'hollow_cubic_shell_rods/SCORES.json').exists():lines+=['','## Literal cubic shell and two rods: geometry correction','', 'The earlier `hollow_cube_rods` folder used a .30×.26×.30λ cuboid with a .18×.16×.18λ void. Those arrays are preserved. This additional case uses actual cubes with outer/inner sides .30/.18λ and the same two rods, 60 dipoles, independent 8³ synthetic/6³ inverse grids and 10% noise. These dimensions are an explicit adaptation because the article prose does not specify them.','',table(load(OUT/'hollow_cubic_shell_rods/SCORES.json')),'','![Actual literal cubic-shell slices](FIGURES/vector3d/literal_cubic_shell.png)']
if (OUT/'refinement8/SCORES.json').exists():lines+=['','## Necessary 8³ inverse / 10³ forward refinement','',table(load(OUT/'refinement8/SCORES.json')),'','![Actual refined material slices](FIGURES/vector3d/four_cube_refinement.png)','','Independent 10³→8³ forward mismatch is 10.607%, versus 15.161% for 8³→6³. The existing per-edge-mean discrete prior is retained; this is not an h⁻²-normalized continuum H1 comparison. See [simulation receipt](RESULTS/vector3d/refinement8/SIMULATION.json).']
if r240:lines+=['','## Preserved 240-cycle PR budget with orthonormal coordinates','',table(r240)]
if (OUT/'pr240_fourier_coordinates/PR_SCORES.json').exists():lines+=['','## Final 240-cycle control with original Fourier coordinates','', 'The original projected IDFT columns are retained without whitening. Sourcewise/global PR uses the original coefficient Euclidean metric. All arms share exact last-paid-state caching, 240 current updates and the 300-step physical corrector cap. This removes the coefficient-preconditioning difference; the port remains an adaptation of the original paper.','',table(load(OUT/'pr240_fourier_coordinates/PR_SCORES.json')),'','![Final original-Fourier-coordinate cost curves](FIGURES/vector3d/final_fourier_cost_curves.png)']
if (OUT/'original_size/FROZEN_CONFIG.json').exists():
    lines+=['','## Original 3λ electrical-size attempt','', 'A separate full-vector padded FFT adapter reuses the A18_C convolution form with the R1 dyadic/CM+RR convention, original60-source/60-receiver geometry and original four-cube dimensions. Synthetic16³/inverse12³ are much coarser than the paper60³/30³. At epsilon5, internal wavelength resolution is inadequate; original electrical size alone does not certify physical accuracy. Direct physical optimization is the only arm attempted at this size.']
    if (OUT/'original_size/TRUNCATION.json').exists():lines.append('Measured attempt TRUNCATED after 20 initial iterations and 3 accepted corrector iterations out of the 30-step allowance. No convergence or completed-corrector claim. '+'The sealed partial endpoint has χ NRMSE .24399 and physical data residual .27632. Independent 16³/12³ mismatch is 36.113%, so physical accuracy is not qualified. Total measured cost: 395.87 s wall, 379.59 s CPU, 356040704-byte peak RSS. The 16³ forward used 5292 FFT matvecs/60 primal RHS; the 12³ process used 232717 FFT matvecs, 3100 primal and 1560 adjoint RHS including the interrupted evaluation. No GMRES RHS failed; maximum paid linear residual was 2.00e−6. See [truncation](RESULTS/vector3d/original_size/TRUNCATION.json), [costs](RESULTS/vector3d/original_size/COSTS.json) and [offline score](RESULTS/vector3d/original_size/SCORE.json).')
    elif (OUT/'original_size/SCORE.json').exists():lines.append('Actual scored result: '+json.dumps(load(OUT/'original_size/SCORE.json')))
    elif (OUT/'original_size/FAILED.json').exists():lines.append('Measured attempt FAILED: '+json.dumps(load(OUT/'original_size/FAILED.json')))
    else:lines.append('Status RUNNING; see original_size/PROGRESS.json. No completion or recovery claim.')
lines+=['','## Interpretation, costs and retained failures','', 'Lower data residual does not establish better material recovery. In the exact-current controls, long physical correction worsens hollow material NRMSE from≈.49 to≈.99–1.02 while fitting the10% noise. Lossy correction also worsens material NRMSE; cross-grid model mismatch remains. Results do not establish A/A+ superiority, a better basin, uniqueness or mature full-scale3D CIE reproduction.','', 'Exact-current process:99.743s wall/97.727s CPU,445923328-byte peak RSS. Per-arm receipts contain online cost; setup/SVD/FFT construction is counted separately inside the relevant method. The first PR process failed while serializing an already computed A+ receipt due to a numpy bool; its failure log was retained, bool conversion repaired, and that arm rerun. The first A+ prediction file was overwritten by the same-contract successful rerun; a separate initial-array snapshot is NOT_RETAINED, and bitwise equivalence is not claimed. Failed-attempt CPU was not instrumented separately and remains NOT_MEASURED; do not silently omit this overhead. Plot/report generation is outside solver timing.','', 'Exact material/Y caching reduced the 8³ LSIE physical LU count from 653 to 350 (303 cache hits), and online CPU from 227.39 to 101.28 s, with the same reported residual. This is an engineering reuse benefit; it supplies no new physical information or independent validation.','', 'The first 8³ attempt was interrupted to repair duplicate callback LU work. Its completed uncached LSIE receipt and partial-Dual OS CPU measurement are retained in `RESULTS/vector3d/refinement8_uncached_attempt/`. The observed 344.46 CPU seconds include the completed 227.39-second LSIE arm; these costs must not be double counted. All four final arms use the same state cache and unchanged 48+300 scientific budgets. The sealed 10³ observations were reused with zero new simulation LU. Engineering overhead remains separate from the online comparison.','', 'Reproduce from repository root with project Python and one BLAS thread: `experiments/vector_benchmark.py`, `vector_pr_benchmark.py`, `vector_pr240_benchmark.py`, `vector_mie_geometry.py`, `vector_original_size.py`, then `vector_report.py`, all under `research/research_max_struct_r8/`. Preserve outputs before deliberate reproduction. No GPU, new service, provider change or R1–R7 modification.']
(ROOT/'VECTOR3D_REPORT.md').write_text('\n'.join(lines)+'\n')
print(ROOT/'VECTOR3D_REPORT.md')

# Auditable successful-arm online costs; setup/simulation and failed jobs remain separate.
ledger=[]
for group,folder in [('exact_current',OUT),('PR48',OUT),('PR240',OUT/'pr240'),('refinement8',OUT/'refinement8'),('original_fourier_PR240',OUT/'pr240_fourier_coordinates'),('literal_cubic_shell',OUT/'hollow_cubic_shell_rods')]:
    if not folder.exists():continue
    pattern='*_receipt.json' if group in ['refinement8','literal_cubic_shell'] else '*/*_receipt.json'
    for file in folder.glob(pattern):
        rr=load(file);kind=rr['method']
        if group in ['exact_current','PR48'] and file.parent.name not in names:continue
        if group=='exact_current' and kind not in ['LSIE_TSOM','CIE_TSOM','CIE_dualbeta_TSOM','A_TSOM','direct_physical']:continue
        if group=='PR48' and kind not in methods:continue
        ledger.append(dict(group=group,target=file.parent.name,method=kind,cpu_s=rr['cpu_s'],wall_s=rr['wall_s'],counts=rr['counts']))
(OUT/'ONLINE_COST_LEDGER.json').write_text(json.dumps(dict(successful_arms=ledger,successful_online_cpu_sum_s=sum(r['cpu_s'] for r in ledger),scope='online successful arms only; sequential CPU can be summed; shared setup/simulation/qualification/original3lambda and failed attempts separately documented; no inferred overhead'),indent=2))

# Original final fields are paid solver states, not new forward validations.
field_rows=[]
for file in OUT.rglob('*_prediction.npz'):
 if 'uncached_attempt' in str(file):continue
 with np.load(file) as a:
  if not all(k in a.files for k in ['prediction','J','E']):continue
  chi=a['prediction'];J=a['J'];E=a['E'];r=J-np.repeat(chi,3)[:,None]*E
  field_rows.append(dict(artifact=str(file.relative_to(ROOT)),paid_state_constitutive_relative_residual=float(np.linalg.norm(r)/max(np.linalg.norm(J),1e-30)),new_forward_LU=0,new_forward_RHS=0,scope='consistency of already paid original full-Maxwell endpoint fields; not independent validation'))
(OUT/'PAID_STATE_FIELD_RESIDUALS.json').write_text(json.dumps(field_rows,indent=2))

# Clarify sharedprepare() metadata after algorithm audit; record original text.
metadata_corrections=[]
for file in OUT.rglob('*_receipt.json'):
 if 'uncached_attempt' in str(file):continue
 rr=load(file)
 if rr.get('subspace') and 'exact reduced quadratic current solve' in rr['subspace'].get('implementation','') and 'current_PR_steps' in rr['counts']:
  metadata_corrections.append(dict(artifact=str(file.relative_to(ROOT)),original=rr['subspace']['implementation'],correction='Subspace construction used SVD orthonormal range of projected IFFT; current update wasPR, as actual PRstepcount/noCholeskyrecords show. OriginalpreparedescriptionincorrectlymentionedexactLS.'))
(OUT/'SUBSPACE_METADATA_CORRECTIONS.json').write_text(json.dumps(metadata_corrections,indent=2))

if (OUT/'refinement8/SCORES.json').exists():
 fig,axes=plt.subplots(4,5,figsize=(12,10));cols=['truth','LSIE_sourcewise','Dual_sourcewise','Aplus_global','direct_physical']
 for row,(ns,zz) in enumerate([(6,-.15),(6,.15),(8,-.15),(8,.15)]):
  folder=OUT/'four_cubes' if ns==6 else OUT/'refinement8';truth=np.load(folder/'truth_for_scoring.npz')['chi'];points=np.load(folder/'truth_for_scoring.npz')['points'];zgrid=np.unique(points[:,2]);iz=int(np.argmin(abs(zgrid-zz)))
  for j,kind in enumerate(cols):
   x=truth if kind=='truth' else np.load(folder/(kind+'_prediction.npz'))['prediction'];axes[row,j].imshow(x.real.reshape((ns,)*3)[:,:,iz].T,origin='lower',vmin=0,vmax=4,cmap='magma',extent=[-.3,.3,-.3,.3]);axes[row,j].set_xticks([]);axes[row,j].set_yticks([])
   if row==0:axes[row,j].set_title(kind,fontsize=9)
   if j==0:axes[row,j].set_ylabel(f'{ns}³ inverse,z={zgrid[iz]:.4f}λ',fontsize=9)
 fig.suptitle('Independent-grid refinement: four dielectric cubes (χ=εr−1)\nOriginal edge prior retained; 6³/8³ use orthonormal current-coordinate adaptation');fig.tight_layout(rect=[0,0,.94,.93]);cax=fig.add_axes([.95,.2,.012,.6]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,4),cmap='magma'),cax=cax,label='Reχ');fig.savefig(FIG/'four_cube_refinement.png',dpi=170);plt.close(fig)
if (OUT/'hollow_cubic_shell_rods/SCORES.json').exists():
 folder=OUT/'hollow_cubic_shell_rods';fig,axes=plt.subplots(1,5,figsize=(12,3));truth=np.load(folder/'truth_for_scoring.npz')['chi'];cols=['truth','LSIE_sourcewise','CIE_global','Aplus_global','direct_physical']
 for j,kind in enumerate(cols):
  x=truth if kind=='truth' else np.load(folder/(kind+'_prediction.npz'))['prediction'];axes[j].imshow(x.real.reshape((6,)*3)[:,:,3].T,origin='lower',vmin=0,vmax=4,cmap='magma',extent=[-.3,.3,-.3,.3]);axes[j].set_title(kind,fontsize=9);axes[j].set_xticks([]);axes[j].set_yticks([])
 fig.suptitle('Literal cubic shell+two rods: original Fourier γ, 48 PR +300 physical-step cap');fig.tight_layout(rect=[0,0,.94,.85]);cax=fig.add_axes([.95,.15,.012,.55]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,4),cmap='magma'),cax=cax,label='Reχ');fig.savefig(FIG/'literal_cubic_shell.png',dpi=170);plt.close(fig)

if (OUT/'pr240_fourier_coordinates/PR_SCORES.json').exists():
 folder=OUT/'pr240_fourier_coordinates/four_cubes';truth=np.load(OUT/'four_cubes/truth_for_scoring.npz')['chi'];fig,axes=plt.subplots(1,2,figsize=(12,4.5));ranges=[[],[]]
 for kind in methods:
  a=np.load(folder/(kind+'_prediction.npz'));r=load(folder/(kind+'_receipt.json'));curve=r['curve'];phi=np.array([q['phi'] for q in curve]);inc=np.array([np.argmin(phi[:j+1]) for j in range(len(phi))]);err=np.linalg.norm(a['trials']-truth,axis=1)/np.linalg.norm(truth);err=err[inc]
  for j,key in enumerate(['physical_LU','cpu_s']):
   xx=np.array([q[key] for q in curve]);axes[j].plot(xx,err,label=kind);ranges[j].append((xx[0],xx[-1]))
 for j in [0,1]:
  lo=max(q[0] for q in ranges[j]);hi=min(q[1] for q in ranges[j]);axes[j].axvspan(lo,hi,color='grey',alpha=.1);axes[j].set_ylabel('χ NRMSE, scored offline after seal');axes[j].grid(alpha=.2)
 axes[0].set_xlabel('Actual cumulative physical LU (PR work excluded)');axes[1].set_xlabel('Cumulative online CPU s (all work)');axes[1].legend(fontsize=7);fig.suptitle('Original projected IDFT γ, 240 PR +300 physical steps\nPhysical-Phi incumbents; grey = observed cost range shared by all arms');fig.tight_layout();fig.savefig(FIG/'final_fourier_cost_curves.png',dpi=170);plt.close(fig)
