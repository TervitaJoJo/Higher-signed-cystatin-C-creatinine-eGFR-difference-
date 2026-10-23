@echo off
call "%~dp0run_review.cmd" --mode full
exit /b %errorlevel%
