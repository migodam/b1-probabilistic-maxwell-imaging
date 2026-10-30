"""Explicit process-local mature-SVD adapter for the SciPy TRF implementation.

At most one large SVD uses cuSOLVER. Other ready CPU SVD work does not queue
behind it. This adapts the linear algebra backend, not TRF's trust-region rule.
"""
import threading
import time
import numpy as np
import torch


class HybridSVD:
    def __init__(self, minimum_columns=256, driver="gesvd"):
        from scipy.optimize._lsq import trf
        self.module = trf
        self.original = trf.svd
        self.minimum_columns = minimum_columns
        self.driver = driver
        self.gpu_slot = threading.BoundedSemaphore(1)
        self.local = threading.local()
        self.rows = []

    def __call__(self, matrix, full_matrices=True, **kwargs):
        m, n = matrix.shape
        t = time.perf_counter()
        use_gpu = n >= self.minimum_columns and not full_matrices and self.gpu_slot.acquire(blocking=False)
        if not use_gpu:
            value = self.original(matrix, full_matrices=full_matrices, **kwargs)
            self.rows.append({"backend": "CPU_gesdd", "m": m, "n": n, "wall_s": time.perf_counter() - t})
            return value
        try:
            if getattr(self.local, "shape", None) != (m, n):
                self.local.shape = (m, n)
                self.local.input = torch.empty((m, n), dtype=torch.float64, pin_memory=True)
                self.local.device = torch.empty((m, n), dtype=torch.float64, device="cuda")
                self.local.output = torch.empty(m * n + n + n * n, dtype=torch.float64, pin_memory=True)
                self.local.done = torch.cuda.Event()
            self.local.input.copy_(torch.from_numpy(np.asarray(matrix)))
            self.local.device.copy_(self.local.input, non_blocking=True)
            u, values, vh = torch.linalg.svd(self.local.device, full_matrices=False, driver=self.driver)
            packed = torch.cat((u.reshape(-1), values, vh.reshape(-1)))
            self.local.output.copy_(packed, non_blocking=True)
            self.local.done.record()
            self.local.done.synchronize()
            host = self.local.output.numpy().copy()
            out = host[:m*n].reshape(m, n), host[m*n:m*n+n], host[m*n+n:].reshape(n, n)
            self.rows.append({"backend": "CUDA_" + self.driver, "m": m, "n": n,
                              "wall_s": time.perf_counter() - t,
                              "transfer_bytes": int((m*n*2 + n + n*n) * 8)})
            return out
        finally:
            self.gpu_slot.release()

    def install(self):
        self.module.svd = self

    def restore(self):
        self.module.svd = self.original

    def summary(self):
        return {"policy": "one nonblocking GPU slot; otherwise CPU; columns>=256; same FP64 economy SVD",
                "minimum_columns": self.minimum_columns, "driver": self.driver,
                "counts": {kind: sum(row["backend"] == kind for row in self.rows)
                           for kind in sorted({row["backend"] for row in self.rows})},
                "calls": self.rows,
                "scope": "SciPy private TRF svd alias in this process only; no frozen source edits; no global system setting"}
