@echo off
rem Runs inside the freshly installed environment. PREFIX is set by constructor.
setlocal

rem The installed prefix is a bare environment: it has no conda and therefore no activate script.
rem Everything here, and the launcher below, works by putting its directories first on PATH.
set "PATH=%PREFIX%;%PREFIX%\Library\bin;%PREFIX%\Scripts;%PATH%"

"%PREFIX%\python.exe" -m pip install --no-deps --no-index --find-links "%PREFIX%" poliscreen
if errorlevel 1 exit /b 1

rem Vina is not a conda package and its 1.2.x series is not on any channel, so the official
rem binary is fetched with its SHA256 verified, exactly as scripts\get_vina.ps1 does.
powershell -NoProfile -ExecutionPolicy Bypass -File "%PREFIX%\scripts\get_vina.ps1" -Dest "%PREFIX%\Scripts"
if errorlevel 1 echo WARNING: Vina was not installed. Run scripts\get_vina.ps1 later.

rem The launcher travels as extra_files, so there is nothing to write here.

endlocal
exit /b 0
