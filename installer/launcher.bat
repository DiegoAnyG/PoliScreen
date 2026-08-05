@echo off
rem Entry point of the installed environment. Without arguments it opens the interface; with them it
rem forwards to the CLI, so `PoliScreen.bat info` checks the installation without activating anything.
setlocal
set "PATH=%~dp0;%~dp0Library\bin;%~dp0Scripts;%PATH%"

if "%~1"=="" (
    "%~dp0Scripts\poliscreen.exe" ui
) else (
    "%~dp0Scripts\poliscreen.exe" %*
)
endlocal
