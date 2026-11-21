"""Wait for the finite inversion queue, seal the cohort, score, and export."""
import time,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
while len(list((root/'predictions').glob('*/SEALED.json')))<48:time.sleep(20)
for script in ('score.py','pack_results.py'):subprocess.run([sys.executable,str(root/script)],check=True)
