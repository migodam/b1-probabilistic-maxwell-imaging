Set-Location D:\AI\FBG_R2_20261030
$env:CUDA_PIPELINE_ROOT='D:\AI\ZALG_CUDA_PIPELINE_20261030'
$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'
[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]21845
Expand-Archive -LiteralPath 'D:\AI\FBG_R2_20261030\development_addons.zip' -DestinationPath 'D:\AI\FBG_R2_20261030' -Force
# PID 21344 is this study's already-running finite Stage4/Stage3 PowerShell.
# Its child PID/command provenance was recorded before this handoff. No kill.
Wait-Process -Id 21344 -Timeout 7200
@{phase='schur_probe_started'; utc=[DateTime]::UtcNow.ToString('o')} | ConvertTo-Json -Compress | Set-Content addon_state.json
& D:\python\python.exe -u -m experiments.schur_feedback_probe schur_probe_remote.json
@{phase='cie_probe_started'; utc=[DateTime]::UtcNow.ToString('o'); previous_exit=$LASTEXITCODE} | ConvertTo-Json -Compress | Set-Content addon_state.json
& D:\python\python.exe -u -m experiments.cie_mechanism_probe cie_probe_remote.json
@{phase='probes_returned'; utc=[DateTime]::UtcNow.ToString('o'); previous_exit=$LASTEXITCODE} | ConvertTo-Json -Compress | Set-Content addon_state.json
