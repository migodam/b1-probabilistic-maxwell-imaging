Set-Location D:\AI\FBG_R2_20261030
$env:CUDA_PIPELINE_ROOT='D:\AI\ZALG_CUDA_PIPELINE_20261030'
$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'
[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]21845
Expand-Archive -LiteralPath 'D:\AI\FBG_R2_20261030\stage4_data.zip' -DestinationPath 'D:\AI\FBG_R2_20261030' -Force
$waitingFbgId = Get-CimInstance Win32_Process | Where-Object { $_.Name -eq 'python.exe' -and $_.CommandLine -like '*src.adaptive_solver stage2_remote.json*' } | Select-Object -First 1 -ExpandProperty ProcessId
if ($waitingFbgId) { Wait-Process -Id $waitingFbgId -Timeout 7200 }
@{phase='stage4_started'; utc=[DateTime]::UtcNow.ToString('o')} | ConvertTo-Json -Compress | Set-Content queue_state.json
& D:\python\python.exe -u -m src.adaptive_solver stage4_remote.json
@{phase='stage4_returned_stage3_started'; utc=[DateTime]::UtcNow.ToString('o'); previous_exit=$LASTEXITCODE} | ConvertTo-Json -Compress | Set-Content queue_state.json
& D:\python\python.exe -u -m src.adaptive_solver stage3_remote.json
@{phase='all_returned'; utc=[DateTime]::UtcNow.ToString('o'); previous_exit=$LASTEXITCODE} | ConvertTo-Json -Compress | Set-Content queue_state.json
