@echo off
REM Navigate to the directory where the script resides
cd /d "%~dp0"
REM Run the Python script
python generate_preference.py
REM Pause to keep the console open
pause
