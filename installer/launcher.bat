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

rem Without this the window closes on its own the moment anything fails, taking the reason with
rem it: it looked like the launcher had merely gone to the background, while the page still being
rem served came from somewhere else entirely.
if errorlevel 1 (
    echo.
    echo PoliScreen stopped. The message above says why; this window stays open so it can be read.
    pause
)
endlocal
