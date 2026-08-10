@echo off
rem Runs inside the freshly installed environment. PREFIX is set by constructor.
setlocal

rem The installed prefix is a bare environment: it has no conda and therefore no activate script.
rem Everything here, and the launcher below, works by putting its directories first on PATH.
set "PATH=%PREFIX%;%PREFIX%\Library\bin;%PREFIX%\Scripts;%PATH%"

"%PREFIX%\python.exe" -m pip install --no-deps --no-index --find-links "%PREFIX%" poliscreen
if errorlevel 1 exit /b 1

rem The design engine, in this same environment. --no-deps because everything it needs is already
rem in the specs; what is missing is ADMET-AI, and admelab reports that itself. Not fatal: the
rem screening runs without it, only the analogue design would go.
"%PREFIX%\python.exe" -m pip install --no-deps --no-index --find-links "%PREFIX%" admelab
if errorlevel 1 echo WARNING: analogue design was not installed.

rem Vina travels inside the installer, with its SHA256 verified in the release workflow. This is
rem only the way back for an installer built without it: downloading it here, on the user's
rem machine, is what used to leave an installation with no docking engine when that one request
rem did not go through, and the warning saying so scrolls past unread.
if not exist "%PREFIX%\Scripts\vina.exe" (
    powershell -NoProfile -ExecutionPolicy Bypass -File "%PREFIX%\scripts\get_vina.ps1" -Dest "%PREFIX%\Scripts"
    if errorlevel 1 echo WARNING: Vina was not installed. Run scripts\get_vina.ps1 later.
)

rem The launcher travels as extra_files, so there is nothing to write here.

endlocal
exit /b 0
