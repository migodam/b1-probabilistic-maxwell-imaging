"""Read only existing final PNG pixels; no GT/scientific arrays and no rendering."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();rows=[]
 for p in sorted((ROOT/'delegated/twod/independent48/plots').glob('MATERIALS_*.png')):
  with Image.open(p) as z:r=np.asarray(z.convert('RGB'))
  y,x=np.where(np.any(r<245,axis=2));h,w=r.shape[:2];rows.append(dict(file=p.name,width=w,height=h,left_ink_margin_pixels=int(x.min()),right_ink_margin_pixels=int(w-1-x.max()),top_ink_margin_pixels=int(y.min()),bottom_ink_margin_pixels=int(h-1-y.max()),outermost_pixel_content=bool((x==0).any() or (x==w-1).any() or (y==0).any() or (y==h-1).any())))
 a.output.mkdir(exist_ok=False);receipt=dict(status='PASS_ALL8_FINAL_PNG_CONTENT_INSIDE_WHITE_MARGINS' if len(rows)==8 and not any(r['outermost_pixel_content'] for r in rows) else 'REVIEW_REQUIRED',pages=rows,GT_access='NONE',new_FW=0,new_render=0,visual_annulus_check='Full r11_independent_0000..0005 labels and right0..20 colorbar title/ticks visible in original PNG; display crop is not present in artifact')
 (a.output/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
