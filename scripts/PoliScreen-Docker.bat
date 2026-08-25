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

rem The banner is base64 so that this file stays seven-bit ASCII. Block-drawing characters in a
rem batch file need the console code page changed, and changing it part-way through shifts the
rem parser's byte offset -- lines split mid-command and every fragment comes back as an
rem unrecognised command, which is exactly how the first release failed to start. cmd never sees
rem the characters here: it passes base64 to PowerShell, which decodes and prints them itself.
cls
echo.
powershell -NoProfile -Command "[Console]::OutputEncoding=[Text.Encoding]::UTF8; [Text.Encoding]::UTF8.GetString([Convert]::FromBase64String('G1szODsyOzA7MjQwOzI1NW0g4paI4paI4paI4paI4paI4paI4pWXICDilojilojilojilojilojilojilZcg4paI4paI4pWXICAgICDilojilojilZfilojilojilojilojilojilojilojilZcg4paI4paI4paI4paI4paI4paI4pWX4paI4paI4paI4paI4paI4paI4pWXIOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKVl+KWiOKWiOKWiOKWiOKWiOKWiOKWiOKVl+KWiOKWiOKWiOKVlyAgIOKWiOKWiOKVlxtbMG0KG1szODsyOzMwOzE5MDsyNTVtIOKWiOKWiOKVlOKVkOKVkOKWiOKWiOKVl+KWiOKWiOKVlOKVkOKVkOKVkOKWiOKWiOKVl+KWiOKWiOKVkSAgICAg4paI4paI4pWR4paI4paI4pWU4pWQ4pWQ4pWQ4pWQ4pWd4paI4paI4pWU4pWQ4pWQ4pWQ4pWQ4pWd4paI4paI4pWU4pWQ4pWQ4paI4paI4pWX4paI4paI4pWU4pWQ4pWQ4pWQ4pWQ4pWd4paI4paI4pWU4pWQ4pWQ4pWQ4pWQ4pWd4paI4paI4paI4paI4pWXICDilojilojilZEbWzBtChtbMzg7Mjs3MDsxNDA7MjU1bSDilojilojilojilojilojilojilZTilZ3ilojilojilZEgICDilojilojilZHilojilojilZEgICAgIOKWiOKWiOKVkeKWiOKWiOKWiOKWiOKWiOKWiOKWiOKVl+KWiOKWiOKVkSAgICAg4paI4paI4paI4paI4paI4paI4pWU4pWd4paI4paI4paI4paI4paI4pWXICDilojilojilojilojilojilZcgIOKWiOKWiOKVlOKWiOKWiOKVlyDilojilojilZEbWzBtChtbMzg7MjsxMjA7OTA7MjU1bSDilojilojilZTilZDilZDilZDilZ0g4paI4paI4pWRICAg4paI4paI4pWR4paI4paI4pWRICAgICDilojilojilZHilZrilZDilZDilZDilZDilojilojilZHilojilojilZEgICAgIOKWiOKWiOKVlOKVkOKVkOKWiOKWiOKVl+KWiOKWiOKVlOKVkOKVkOKVnSAg4paI4paI4pWU4pWQ4pWQ4pWdICDilojilojilZHilZrilojilojilZfilojilojilZEbWzBtChtbMzg7MjsxNzA7NDA7MjU1bSDilojilojilZEgICAgIOKVmuKWiOKWiOKWiOKWiOKWiOKWiOKVlOKVneKWiOKWiOKWiOKWiOKWiOKWiOKWiOKVl+KWiOKWiOKVkeKWiOKWiOKWiOKWiOKWiOKWiOKWiOKVkeKVmuKWiOKWiOKWiOKWiOKWiOKWiOKVl+KWiOKWiOKVkSAg4paI4paI4pWR4paI4paI4paI4paI4paI4paI4paI4pWX4paI4paI4paI4paI4paI4paI4paI4pWX4paI4paI4pWRIOKVmuKWiOKWiOKWiOKWiOKVkRtbMG0KG1szODsyOzE3MDs0MDsyNTVtIOKVmuKVkOKVnSAgICAgIOKVmuKVkOKVkOKVkOKVkOKVkOKVnSDilZrilZDilZDilZDilZDilZDilZDilZ3ilZrilZDilZ3ilZrilZDilZDilZDilZDilZDilZDilZ0g4pWa4pWQ4pWQ4pWQ4pWQ4pWQ4pWd4pWa4pWQ4pWdICDilZrilZDilZ3ilZrilZDilZDilZDilZDilZDilZDilZ3ilZrilZDilZDilZDilZDilZDilZDilZ3ilZrilZDilZ0gIOKVmuKVkOKVkOKVkOKVnRtbMG0='))"
echo.
echo      Reproducible virtual screening   v1.1.0   container setup
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
        rem A pull that replaced the image leaves the previous one untagged. Only this
        rem repository's dangling images are removed -- a blanket prune would take other
        rem projects' with it, on a machine that is not ours to tidy.
        for /f %%i in ('docker images ghcr.io/diegoanyg/poliscreen --filter "dangling=true" -q') do (
            docker rmi %%i >nul 2>&1
        )
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

rem A container from a previous window keeps port 8501, so the browser opened to the build that
rem was already there and nothing said the new image had not started. Naming it makes that
rem visible; stopping it is asked for rather than done, because it may be someone's running
rem screening.
docker ps -q --filter "name=poliscreen" > "%TEMP%\poliscreen_running.txt"
for /f %%i in (%TEMP%\poliscreen_running.txt) do set "RUNNING=%%i"
del "%TEMP%\poliscreen_running.txt" >nul 2>&1
if defined RUNNING (
    echo.
    echo   PoliScreen is already running in another window, on this same port.
    echo   That older container is what your browser will show.
    echo.
    choice /C YN /M "Stop it and start this one"
    if errorlevel 2 (
        echo   Left running. Open http://localhost:8501 to use it.
        pause
        exit /b 0
    )
    docker stop poliscreen >nul 2>&1
)

docker run --rm -it --init --name poliscreen -p 127.0.0.1:8501:8501 -v "%PROJECTS%:/data" %THEME% %IMAGE%

endlocal
