"""Two-page post-scoring atlas of the ten fixed blind reconstructions."""
from collections import defaultdict
import csv
from decimal import Decimal, ROUND_FLOOR
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt


ROOT = Path(__file__).resolve().parents[1]
MAXIMUM = Decimal("4.5")
CONTRAST = "magma"
COLORS = {
    "outside": "#eeeee9",
    "full": "#1a796e",
    "partial": "#e8aa54",
    "single": "#9b85b8",
    "unresolved": "#899aa9",
}
INK = "#172f39"
MUTED = "#465e66"
RULE = "#d8e2df"
PAPER = "#fbfaf6"
PAGE_WIDTH = 7.2
PAGE_HEIGHT = 16.45
MARGIN = .38
GAP = .15
PANEL = (PAGE_WIDTH-2*MARGIN-2*GAP)/3
LEFT = [MARGIN+j*(PANEL+GAP) for j in range(3)]


def setup_font():
    fonts = {entry.name for entry in font_manager.fontManager.ttflist}
    for name in ("Noto Sans CJK SC", "Arial Unicode MS", "Heiti SC", "WenQuanYi Zen Hei"):
        if name in fonts:
            plt.rcParams.update({"font.family": [name, "DejaVu Sans"], "pdf.fonttype": 42})
            return
    raise RuntimeError("Install a Simplified Chinese font, for example Noto Sans CJK SC")


def read_cases():
    grouped = defaultdict(list)
    with (ROOT/"results"/"region_intervals.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            grouped[row["case"]].append(row)
    if len(grouped) != 10:
        raise ValueError("The finalized atlas expects ten independent objects")

    for case, rows in sorted(grouped.items()):
        with np.load(ROOT/"imaging"/"results"/case/"stationary"/"baseline.npz") as result:
            masks = result["rois"].astype(bool)
            image = result["chi"]
        with np.load(ROOT/"ports"/"private"/f"{case}.npz") as scoring:
            truth = scoring["material"].reshape(image.shape)
        if len(masks) != len(rows) or np.any(masks.sum(axis=0) > 1):
            raise ValueError(f"Final ROI partitions do not match {case}")
        rows = sorted(rows,key=lambda r: int(r["roi"]))
        for i,(mask,row) in enumerate(zip(masks,rows,strict=True),1):
            if int(row["roi"]) != i or int(mask.sum()) != int(row["cells"]):
                raise ValueError(f"Frozen ROI mismatch for {case} ROI {i}")
            if not (Decimal(row["outer_lower"]) == 0 and Decimal(row["outer_upper"]) == MAXIMUM):
                raise ValueError(f"Update the shared outer-interval label for {case} ROI {i}")
            if not np.isclose(truth[mask].max(),float(row["true_peak"]),atol=1e-12,rtol=0):
                raise ValueError(f"Post-freeze truth panel mismatch for {case} ROI {i}")
        yield case,truth,image,masks,rows


def evidence(row):
    low = Decimal(row["known_attained_min"]) if row["known_attained_min"] else None
    high = Decimal(row["known_attained_max"]) if row["known_attained_max"] else None
    if low is None and high is None:
        return "unresolved","未证得"
    if low is None or high is None or low == high:
        point = low if low is not None else high
        return "single",f"单点 {point:.3g}"
    if not 0 <= low < high <= MAXIMUM:
        raise ValueError(f"Invalid attained peaks in {row['case']} ROI {row['roi']}")
    if low == 0 and high == MAXIMUM:
        if row["width_exact"] != "True":
            raise ValueError("Full-width endpoints must match the frozen table")
        return "full","4.5"
    span = (high-low).quantize(Decimal("0.01"),rounding=ROUND_FLOOR)
    return "partial","≥"+format(span.normalize(),"f")


def text(figure,x,y,content,**kwargs):
    defaults = {"color":INK,"va":"center","ha":"left"}
    figure.text(x/PAGE_WIDTH,y/PAGE_HEIGHT,content,**(defaults|kwargs))


def swatch(figure,x,y,color):
    figure.add_artist(Rectangle((x/PAGE_WIDTH,(y-.066)/PAGE_HEIGHT),.132/PAGE_WIDTH,.132/PAGE_HEIGHT,
                    transform=figure.transFigure,facecolor=color,edgecolor=RULE,linewidth=.6))


def panel(figure,x,y):
    ax = figure.add_axes([x/PAGE_WIDTH,y/PAGE_HEIGHT,PANEL/PAGE_WIDTH,PANEL/PAGE_HEIGHT])
    ax.set(xticks=[],yticks=[])
    for spine in ax.spines.values():
        spine.set(color=RULE,linewidth=.8)
    return ax


def annotate_regions(ax,masks):
    size = masks.shape[-1]
    used = []
    for i,mask in enumerate(masks,1):
        distance = distance_transform_edt(np.pad(mask,1))[1:-1,1:-1]
        candidates = np.argwhere(distance==distance.max())
        cy,cx = np.argwhere(mask).mean(axis=0)
        ay,ax0 = min(candidates,key=lambda yx: (yx[0]-cy)**2+(yx[1]-cx)**2)
        anchor = np.array([ax0,ay])/size
        x,y = np.clip(anchor,.10,.89)
        trial = [np.array([x,y])]
        for radius in (.12,.21,.31,.46):
            for phi in np.linspace(0,2*np.pi,16,endpoint=False):
                proposal = np.array([x,y])+radius*np.array([np.cos(phi),np.sin(phi)])
                if np.all((.10<=proposal)&(proposal<=.89)):
                    trial.append(proposal)
        # Small masks remain visible next to an externally placed number badge.
        if mask.sum() < .027*size**2:
            outside = [p for p in trial if np.linalg.norm(p-anchor)>=.12
                       and not masks[:,int(p[1]*size),int(p[0]*size)].any()]
            if outside:
                trial = outside + trial
        position = next((p for p in trial if all(np.linalg.norm(p-u)>.16 for u in used)),trial[0])
        used.append(position)
        point,tip = size*position,size*anchor
        leader = dict(arrowstyle="-",color=INK,lw=.9,shrinkA=2,shrinkB=0) if np.linalg.norm(position-anchor)>.027 else None
        ax.annotate(str(i),xy=tip,xytext=point,ha="center",va="center",
                    color=INK,fontsize=12.5,fontfamily="DejaVu Sans",fontweight="bold",
                    bbox=dict(boxstyle="circle,pad=.24",fc="white",ec=INK,lw=.7),
                    arrowprops=leader,zorder=4)


def draw_case(figure, case_data, top):
    case, truth, image, masks, rows = case_data
    text(figure,MARGIN,top,f"{case}  /  {image.shape[0]}×{image.shape[1]}",
         fontfamily="DejaVu Sans",fontsize=12.5,fontweight="bold")
    y = top-.20-PANEL
    ax = [panel(figure,x,y) for x in LEFT]
    for a,field in zip(ax[:2],(truth,image),strict=True):
        a.imshow(field,origin="lower",cmap=CONTRAST,vmin=0,vmax=float(MAXIMUM),interpolation="nearest")
    state = np.zeros(image.shape,np.uint8)
    classes = list(COLORS)
    for index,(row,mask) in enumerate(zip(rows,masks,strict=True)):
        kind,label = evidence(row)
        state[mask] = classes.index(kind)
        # White contours on the blind image refer to the same fixed masks.
        ax[1].contour(mask.astype(float),levels=[.5],colors=["#d2f5f0"],linewidths=.54)
        # The number under each row is the certified span between the extrema.
        x = LEFT[index]
        swatch(figure,x,y-.21,COLORS[kind])
        text(figure,x+.19,y-.21,f"ROI {index+1} · {label}",fontsize=12.5)
    ax[2].imshow(state,origin="lower",cmap=ListedColormap(list(COLORS.values())),
                 vmin=0,vmax=len(classes)-1,interpolation="nearest")
    annotate_regions(ax[2],masks)


def draw_sheet(cases,page):
    figure = plt.figure(figsize=(PAGE_WIDTH,PAGE_HEIGHT),facecolor=PAPER)
    top = PAGE_HEIGHT-.36
    text(figure,MARGIN,top,"十件对象的盲区域图集",fontsize=24)
    text(figure,PAGE_WIDTH-MARGIN,top,f"{page} / 2",ha="right",fontsize=12,color=MUTED)
    text(figure,MARGIN,top-.40,"每个编号区域的保证外包均为 [0, 4.5]，宽 4.5。",fontsize=14)
    legend = [
        (.00,"full","全宽见证"), (1.38,"partial","部分宽度"), (2.76,"single","单点见证"),
        (4.14,"unresolved","无见证"), (5.38,"outside","区外"),
    ]
    for dx,kind,label in legend:
        swatch(figure,MARGIN+dx,top-.78,COLORS[kind])
        text(figure,MARGIN+dx+.19,top-.78,label,fontsize=12.5)
    bar = figure.add_axes([MARGIN/PAGE_WIDTH,(top-1.12)/PAGE_HEIGHT,
                          (2*PANEL+GAP)/PAGE_WIDTH,.135/PAGE_HEIGHT])
    cb = figure.colorbar(ScalarMappable(norm=Normalize(0,float(MAXIMUM)),cmap=CONTRAST),
                         cax=bar,orientation="horizontal",ticks=[0,1.5,3,4.5])
    cb.ax.tick_params(length=2,labelsize=11,pad=2,colors=MUTED)
    cb.outline.set_linewidth(0)
    text(figure,LEFT[2],top-1.16,"前两列共用对比度色标",fontsize=12,color=MUTED)
    for x,label in zip(LEFT,("真值","盲重建","区域证据"),strict=True):
        text(figure,x+PANEL/2,top-1.60,label,ha="center",fontsize=15)
    next_top = top-1.95
    for row in cases:
        draw_case(figure,row,next_top)
        next_top -= PANEL+.72
    text(figure,MARGIN,.18,"数字给出两端见证的宽度下限。整段可行需要单独证明。",fontsize=11.5,color=MUTED)
    return figure


def main():
    setup_font()
    cases = list(read_cases())
    out = ROOT/"figures"
    with PdfPages(out/"roi_atlas.pdf",metadata={"Title":"十件对象盲区域图集"}) as pdf:
        for page in (1,2):
            fig = draw_sheet(cases[(page-1)*5:page*5],page)
            pdf.savefig(fig)
            fig.savefig(out/f"roi_atlas_{page}.png",dpi=170,facecolor=PAPER)
            plt.close(fig)


if __name__ == "__main__":
    main()
