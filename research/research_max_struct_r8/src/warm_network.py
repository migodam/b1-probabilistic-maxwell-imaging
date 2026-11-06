"""Small two-scale CNN; material-only prediction, no solver or test labels."""
import torch
from torch import nn
import torch.nn.functional as F

class TinyWarmStart(nn.Module):
    def __init__(self,channels=10,width=16):
        super().__init__()
        self.enc=nn.Sequential(nn.Conv2d(channels,width,3,padding=1),nn.SiLU(),nn.Conv2d(width,width,3,padding=1),nn.SiLU())
        self.low=nn.Sequential(nn.Conv2d(width,width*2,3,padding=1),nn.SiLU(),nn.Conv2d(width*2,width*2,3,padding=1),nn.SiLU())
        self.dec=nn.Sequential(nn.Conv2d(width*3,width,3,padding=1),nn.SiLU(),nn.Conv2d(width,1,3,padding=1))
        nn.init.constant_(self.dec[-1].bias,-2.)

    def forward(self,x):
        hi=self.enc(x);low=self.low(F.avg_pool2d(hi,2))
        z=self.dec(torch.cat((hi,F.interpolate(low,size=hi.shape[-2:],mode='bilinear',align_corners=False)),dim=1))
        return .02+19.98*torch.sigmoid(z[:,0])

def load_predictor(path,device='cpu'):
    checkpoint=torch.load(path,map_location=device,weights_only=False)
    net=TinyWarmStart(**checkpoint['architecture']).to(device);net.load_state_dict(checkpoint['state_dict']);net.eval()
    mean=torch.tensor(checkpoint['mean'],device=device)[None,:,None,None]
    std=torch.tensor(checkpoint['std'],device=device)[None,:,None,None]
    def predict(features):
        x=torch.as_tensor(features,dtype=torch.float32,device=device)
        if x.ndim==3:x=x[None]
        if checkpoint['feature_mode']=='born':
            x=x.clone();x[:,2:6]=0
        with torch.no_grad():out=net((x-mean)/std)
        return out.cpu().numpy()
    return predict,checkpoint
