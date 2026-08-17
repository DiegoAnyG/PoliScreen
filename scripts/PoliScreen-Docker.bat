@echo off
rem One double-click for the Docker route, because "download an exe and run it" was the whole
rem reason the Windows build existed. This gets the same experience without the native build's
rem hydrogen-bond divergence: no git, no terminal, no 10-minute conda solve. The image is already
rem built and published, so the first run is a download and every run after that is instant.
setlocal
rem :latest is only published by a tagged release; every other build publishes :edge. Before the
rem first tag exists, asking for :latest fails with a message about the tag rather than about the
rem release that has not happened yet, so fall back rather than leave a double-click dead.
set "IMAGE=ghcr.io/diegoanyg/poliscreen:latest"
set "FALLBACK=ghcr.io/diegoanyg/poliscreen:edge"
set "PROJECTS=%USERPROFILE%\PoliScreen"

rem A UTF-8 code page and a real ESC character, so the block letters and the colour
rem gradient render instead of arriving as literal escape sequences. The prompt trick is
rem the only way to get ESC into a variable in batch. Delayed expansion is deliberately
rem NOT enabled: it would eat any exclamation mark in the messages further down.
chcp 65001 >nul
for /F %%a in ('"prompt $E$S & for %%b in (1) do rem"') do set "ESC=%%a"
set "C1=%ESC%[38;2;0;240;255m"
set "C2=%ESC%[38;2;30;190;255m"
set "C3=%ESC%[38;2;60;140;255m"
set "C4=%ESC%[38;2;100;90;255m"
set "C5=%ESC%[38;2;150;40;255m"
set "C6=%ESC%[38;2;150;40;255m"
set "DIM=%ESC%[38;2;120;160;170m"
set "RESET=%ESC%[0m"

cls
echo.
echo %C1% ██████╗  ██████╗ ██╗     ██╗███████╗ ██████╗██████╗ ███████╗███████╗███╗   ██╗%RESET%
echo %C2% ██╔══██╗██╔═══██╗██║     ██║██╔════╝██╔════╝██╔══██╗██╔════╝██╔════╝████╗  ██║%RESET%
echo %C3% ██████╔╝██║   ██║██║     ██║███████╗██║     ██████╔╝█████╗  █████╗  ██╔██╗ ██║%RESET%
echo %C4% ██╔═══╝ ██║   ██║██║     ██║╚════██║██║     ██╔══██╗██╔══╝  ██╔══╝  ██║╚██╗██║%RESET%
echo %C5% ██║     ╚██████╔╝███████╗██║███████║╚██████╗██║  ██║███████╗███████╗██║ ╚████║%RESET%
echo %C6% ╚═╝      ╚═════╝ ╚══════╝╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝╚══════╝╚══════╝╚═╝  ╚═══╝%RESET%
echo.
echo %DIM%      Reproducible virtual screening   v1.0.1   container setup%RESET%
echo.

rem Docker Desktop is never started from here on purpose: it is the user's machine and starting a
rem background service behind their back is not this script's business. Say what is wrong instead.
where docker >nul 2>&1
if errorlevel 1 (
    echo   Docker is not installed.
    echo.
    echo   Install Docker Desktop, then run this again:
    echo     https://www.docker.com/products/docker-desktop
    echo.
    echo   It is the only thing you have to install. Everything else travels
    echo   inside the image: Python, Vina, RDKit, PLIP, fpocket, the ADMET engine.
    echo.
    pause
    exit /b 1
)

docker version >nul 2>&1
if errorlevel 1 (
    echo   Docker is installed but not running.
    echo.
    echo   Start Docker Desktop and wait for the whale icon to stop animating,
    echo   then run this again.
    echo.
    pause
    exit /b 1
)
echo   Docker              OK

if not exist "%PROJECTS%" mkdir "%PROJECTS%"
echo   Projects folder     %PROJECTS%
echo.

docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
    echo   First run: downloading the image, about 3 GB. This happens once;
    echo   every run after this one starts in seconds.
    echo.
    docker pull %IMAGE%
    if errorlevel 1 (
        echo   No released image yet, trying the development one.
        docker pull %FALLBACK% || (
            echo.
            echo   Could not download the image.
            echo   If the message above says "denied" or "unauthorized", the package is
            echo   private: make it public once, at
            echo     https://github.com/users/DiegoAnyG/packages/container/poliscreen/settings
            echo.
            pause
            exit /b 1
        )
        set "IMAGE=%FALLBACK%"
    )
) else (
    rem Already downloaded. Checking for a newer one costs a few kilobytes -- the manifest, not
    rem the layers -- and skipping it is how somebody keeps running a version fixed months ago
    rem without ever being told. Offline, the copy on disk is used and nothing is said, because
    rem being offline is not an error here.
    echo   Checking for updates...
    docker pull %IMAGE% >nul 2>&1
    if errorlevel 1 (
        echo   Image               local copy (no network^)
    ) else (
        echo   Image               up to date
    )
)

echo.
echo   ---------------------------------------------------------------
echo   Ready. The interface opens at  http://localhost:8501
echo   Close this window to stop PoliScreen.
echo   ---------------------------------------------------------------
echo.

rem Opened on a delay because the page only answers once streamlit is up; opening it immediately
rem showed a connection error and looked like a failure.
start "" /b cmd /c "timeout /t 8 >nul & start "" http://localhost:8501"

rem --init so Ctrl+C reaches streamlit instead of being ignored by PID 1, --rm so nothing is left
rem behind, and the port bound to 127.0.0.1 because the interface has no authentication.
rem Dark interface: set POLISCREEN_THEME=dark before running this, or add the line here. The
rem menu in Streamlit 1.59 no longer carries the switcher when the app defines its own theme, so
rem the choice is made when the container starts.
set "THEME="
if defined POLISCREEN_THEME set "THEME=-e STREAMLIT_THEME_BASE=%POLISCREEN_THEME%"

docker run --rm -it --init -p 127.0.0.1:8501:8501 -v "%PROJECTS%:/data" %THEME% %IMAGE%

endlocal
