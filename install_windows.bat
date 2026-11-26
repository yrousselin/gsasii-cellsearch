@echo off
REM Installs the unit-cell search plugin into GSAS-II (Windows)
set "DEST=%USERPROFILE%\.GSASII\exports"
if not exist "%DEST%" mkdir "%DEST%"
copy /Y "%~dp0G2export_CellSearch.py" "%DEST%\" >nul
copy /Y "%~dp0g2_cellsearch.py" "%DEST%\" >nul
copy /Y "%~dp0csd_cellsearch_worker.py" "%DEST%\" >nul
echo.
echo Files copied to: %DEST%
echo Restart GSAS-II:
echo   Export - Entire project as - (user) Search databases by unit cell
echo.
pause
