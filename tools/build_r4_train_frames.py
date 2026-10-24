"""Freeze an actual TRAIN-only prior-relative Maxwell frame and random control."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,save,write
from p3r4.conditional_coordinates import native_to_gaussianized_jacobian,fixed_random_rotation

def main():
    dest=OUT/'inputs/train_frames_v1'
    if dest.exists():raise FileExistsError('TRAIN frame already frozen')
    book=Book('I','build_train_frames_v1',analysis_only=True,max_wall=120)
    status='COMPLETE';error=None
    try:
        dest.mkdir(parents=True);H=np.zeros((66,66));identities=[]
        for k in range(256):
            p=OUT/'results/contexts/train'/f'R2_train_{k:05d}.npz'
            stat=p.stat();identities.append(dict(path=str(p.resolve()),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
            with np.load(p,allow_pickle=False) as z:
                J=native_to_gaussianized_jacobian(z['x'][None],z['J'][None])[0]
            H+=J.T@J/256
        values,frame=np.linalg.eigh((H+H.T)/2)
        signs=np.sign(frame[np.argmax(np.abs(frame),axis=0),np.arange(66)]);frame*=np.where(signs==0,1,signs)
        save(dest/'maxwell_frame.npz',frame=frame,mean_likelihood_spectrum=values)
        save(dest/'random_frame.npz',frame=fixed_random_rotation(20261025113))
        write(dest/'PROVENANCE.json',dict(status='TRAIN_ONLY_FRAME_FROZEN',n_train=256,dimension=66,weak_first=16,
            criterion='eigenvectors of mean whitened Gaussianized Jx.T Jx; prior-relative frame, not material-risk eigenvectors',
            original_acquisition_cost='paid StageI context reconstruction, not free deployment',truth_read=False,
            orthogonality_error=float(np.linalg.norm(frame.T@frame-np.eye(66))),sources=identities,
            random_seed=20261025113,new_complex_RHS=0,new_LU=0))
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:book.close(status,error)
if __name__=='__main__':main()
