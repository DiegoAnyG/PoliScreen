@echo off
rem One double-click for the Docker route, because "download an exe and run it" was the whole
rem reason the Windows build existed. This gets the same experience without the native build's
rem hydrogen-bond divergence: no git, no terminal, no 10-minute conda solve. The image is already
rem built and published, so the first run is a download and every run after that is instant.
setlocal
set "IMAGE=ghcr.io/diegoanyg/poliscreen:latest"
set "PROJECTS=%USERPROFILE%\PoliScreen"

echo PoliScreen (container)
echo.

rem Docker Desktop is never started from here on purpose: it is the user's machine and starting a
rem background service behind their back is not this script's business. Say what is wrong instead.
docker version >nul 2>&1
if errorlevel 1 (
    echo Docker is not responding.
    echo.
    echo   - If Docker Desktop is installed, start it and wait for the whale icon to stop animating.
    echo   - If it is not, install it from https://www.docker.com/products/docker-desktop
    echo.
    pause
    exit /b 1
)

if not exist "%PROJECTS%" mkdir "%PROJECTS%"

docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
    echo First run: downloading the image. This happens once.
    docker pull %IMAGE% || (echo Could not download the image. & pause & exit /b 1)
)

echo Projects are saved in %PROJECTS%
echo Opening http://localhost:8501 - close this window to stop PoliScreen.
echo.

rem Opened on a delay because the page only answers once streamlit is up; opening it immediately
rem showed a connection error and looked like a failure.
start "" /b cmd /c "timeout /t 8 >nul & start "" http://localhost:8501"

rem --init so Ctrl+C reaches streamlit instead of being ignored by PID 1, --rm so nothing is left
rem behind, and the port bound to 127.0.0.1 because the interface has no authentication.
docker run --rm -it --init -p 127.0.0.1:8501:8501 -v "%PROJECTS%:/data" %IMAGE%

endlocal
