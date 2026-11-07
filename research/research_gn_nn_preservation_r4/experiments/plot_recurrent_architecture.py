"""Draw the implemented recurrence; no observations, truth or solver calls."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Rectangle

S=Path(__file__).resolve().parents[1]
fig,ax=plt.subplots(figsize=(15,6.4));ax.set(xlim=(0,15),ylim=(0,6.4));ax.axis('off')
def box(x,y,w,h,title,body,color='#eaf4f8'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.10',facecolor=color,edgecolor='#31566a',linewidth=1.3))
    ax.text(x+w/2,y+h-.28,title,ha='center',va='center',fontsize=12,weight='bold')
    ax.text(x+w/2,y+h/2-.12,body,ha='center',va='center',fontsize=10,linespacing=1.6)
def arrow(a,b,**kw):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=15,linewidth=1.5,color='#31566a',**kw))
ax.add_patch(Rectangle((.55,2.05),13.9,3.25,facecolor='#f7fbfd',edgecolor='#79a8bd',linestyle='--'))
ax.text(7.5,6.04,'A17-inspired current–material feedback: the implemented recurrence',ha='center',fontsize=18,weight='bold')
ax.text(7.5,5.62,'Two trained shared-weight steps; fixed Green/receiver multiplications inside the loop; no full Maxwell inverse layer',ha='center',fontsize=11)
box(.85,3.05,2.3,1.55,'Learned proposal',r'$j_a=j^k+N_\theta(m^k,Q^k,\mathrm{gb},\nu)$'+'\nShared U-Net weights')
box(3.65,3.05,2.5,1.55,'Shared real material',r'$q_a=E+Gj_a$'+'\n'+r'$m_a=\mathrm{RowLS}(j_a,q_a)$')
box(6.65,3.05,3.,1.55,'Physics feedback',r'$e=j_a-D_{m_a}q_a,\quad r=Sj_a-Y$'+'\n'+r'$g=S^*r/a+\rho(e-G^*D_{m_a}e)/b_0$')
box(10.15,3.05,3.85,1.55,'Current / material refresh',r'$j^{k+1}=j_a-\tau g$'+'\n'+r'$Q^{k+1}=E+Gj^{k+1}$'+'\n'+r'$m^{k+1}=\mathrm{RowLS}(j^{k+1},Q^{k+1})$')
for a,b in (((3.15,3.83),(3.65,3.83)),((6.15,3.83),(6.65,3.83)),((9.65,3.83),(10.15,3.83))):arrow(a,b)
arrow((12,3.02),(12,2.49));arrow((12,2.49),(1.95,2.49));arrow((1.95,2.49),(1.95,3.02))
ax.text(6.8,2.55,'Carry current and material; refresh auxiliary field',ha='center',va='bottom',fontsize=10)
box(.9,.40,4.7,1.08,'Public strict GN initialization',r'$m^0,\ j^0,\ Q^0,\ \mathrm{gb},\ Y,\nu$'+'\nFull Maxwell LU is paid here',color='#f5eee4')
arrow((2,1.5),(2,2.08))
box(8.05,.40,5.9,1.08,'Outside the recurrent network',r'Full Maxwell $F(m)$ verification + common TRF correction'+'\nPhysical top1 selection; GT scoring only after predictions seal',color='#f5eee4')
arrow((12,2.03),(12,1.51))
ax.text(7.5,.08,'Quadratic feedback shown. Hinge control changes only the data gradient. Neither the two-step network nor deeper replay has a fixed-point guarantee.',ha='center',fontsize=9)
out=S/'figures/recurrent';out.mkdir(parents=True,exist_ok=True)
for ext in ('png','svg'):fig.savefig(out/f'implemented_recurrent_architecture.{ext}',dpi=160,bbox_inches='tight')
plt.close(fig)
print(out/'implemented_recurrent_architecture.png')
