"""Observation-to-material point estimators. No posterior calibration is implied."""
import numpy as np
import torch
from torch import nn


def observation_features(observation):
    """Source-major Re32 then Im32, matching the frozen measurement contract."""
    a = np.asarray(observation)
    if a.ndim != 3 or a.shape[1:] != (32, 8) or not np.iscomplexobj(a):
        raise ValueError('complex observation [N,32,8] required')
    return np.concatenate((a.real.transpose(0, 2, 1), a.imag.transpose(0, 2, 1)), axis=2).reshape(len(a), 512).astype('float32')


class MaterialRegressor(nn.Module):
    def __init__(self, kind, output_dim=512):
        super().__init__()
        self.kind = kind
        if kind in ('mlp_small', 'mlp_medium'):
            widths = [512, 256, 128] if kind == 'mlp_small' else [512, 512, 256, 256]
            layers = []
            for a, b in zip(widths[:-1], widths[1:]):
                layers += [nn.Linear(a, b), nn.GELU(), nn.Dropout(0.05)]
            self.encoder = nn.Sequential(*layers)
            self.head = nn.Linear(widths[-1], output_dim)
        elif kind in ('cnn', 'cnn_context'):
            self.encoder = nn.Sequential(nn.Conv2d(2, 16, 3, padding=1), nn.GELU(),
                                         nn.Conv2d(16, 32, 3, stride=(2, 1), padding=1), nn.GELU(),
                                         nn.Conv2d(32, 32, 3, stride=2, padding=1), nn.GELU(),
                                         nn.AdaptiveAvgPool2d((4, 2)), nn.Flatten())
            extra = 0
            if kind == 'cnn_context':
                # Absolute receiver/source position and signed phase are preserved
                # in a parallel global measurement branch, using observations only.
                self.context = nn.Sequential(nn.Linear(512, 128), nn.GELU(), nn.Dropout(0.05))
                extra = 128
            self.head = nn.Sequential(nn.Linear(256 + extra, 256), nn.GELU(), nn.Linear(256, output_dim))
        else:
            raise ValueError(kind)
        last = self.head if isinstance(self.head, nn.Linear) else self.head[-1]
        # Initialization is the train material mean, a meaningful update-zero baseline.
        nn.init.zeros_(last.weight)
        nn.init.zeros_(last.bias)

    def forward(self, x):
        if self.kind.startswith('cnn'):
            im = x.reshape(-1, 8, 2, 32).permute(0, 2, 3, 1)
            h = self.encoder(im)
            if self.kind == 'cnn_context':
                h = torch.cat((h, self.context(x)), dim=1)
        else:
            h = self.encoder(x)
        return self.head(h)


class RiskPredictor:
    """Reusable locked inference; accepts observations and never material truth."""
    def __init__(self, payload, device='cpu'):
        self.payload = payload
        self.device = torch.device(device)
        self.kind = payload['kind']
        self.roi_indices = np.asarray(payload['roi_indices'], dtype=int)
        self.model = None
        if self.kind not in ('mean', 'ridge'):
            self.model = MaterialRegressor(self.kind, 2 * len(self.roi_indices)).to(self.device)
            self.model.load_state_dict(payload['state_dict'])
            self.model.eval()

    def predict(self, observation, batch_size=256):
        p = self.payload
        x = (observation_features(observation) - p['x_mean']) / p['x_scale']
        output = []
        with torch.inference_mode():
            for start in range(0, len(x), batch_size):
                b = x[start:start + batch_size]
                if self.kind == 'mean':
                    y = np.broadcast_to(p['y_mean'], (len(b), len(p['y_mean'])))
                elif self.kind == 'ridge':
                    y = b @ p['coef'] + p['y_mean']
                else:
                    y = self.model(torch.from_numpy(b).to(self.device)).cpu().numpy() * p['y_scale'] + p['y_mean']
                c = len(self.roi_indices)
                output.append(y[:, :c] + 1j * y[:, c:])
        return np.concatenate(output)

    def predict_full(self, observation, batch_size=256):
        roi = self.predict(observation, batch_size)
        full = np.full((len(roi), 576), complex(*self.payload['host_chi']), dtype=np.complex64)
        full[:, self.roi_indices] = roi
        return full


def load_predictor(path, device='cpu'):
    # Own checkpoints contain normalization arrays as well as tensor weights.
    return RiskPredictor(torch.load(path, map_location='cpu', weights_only=False), device)
