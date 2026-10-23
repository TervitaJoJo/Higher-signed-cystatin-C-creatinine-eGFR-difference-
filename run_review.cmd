@echo off
setlocal
set PYTHONDONTWRITEBYTECODE=1
set PYTHONUTF8=1
where python >nul 2>nul
if errorlevel 1 (
  echo Python is required. See README.md.
  exit /b 1
)
python "%~dp0reproduce.py" %*
exit /b %errorlevel%
