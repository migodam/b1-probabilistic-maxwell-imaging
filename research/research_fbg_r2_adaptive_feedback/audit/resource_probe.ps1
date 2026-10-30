Get-CimInstance Win32_Process | Where-Object { $_.Name -eq 'python.exe' } | Select-Object ProcessId,CommandLine | ConvertTo-Json -Compress
nvidia-smi --query-gpu=name,memory.used,utilization.gpu --format=csv
New-Item -ItemType Directory -Force D:\AI\FBG_R2_20261030 | Out-Null
