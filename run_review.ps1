param([ValidateSet('figures','full','public')][string]$Mode='figures')
$env:PYTHONUTF8='1'
$env:PYTHONDONTWRITEBYTECODE='1'
& python (Join-Path $PSScriptRoot 'reproduce.py') --mode $Mode
exit $LASTEXITCODE
