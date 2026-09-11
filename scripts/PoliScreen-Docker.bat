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

rem Create desktop shortcut with icon
if not exist "%~dp0PoliScreen.ico" (
    powershell -NoProfile -Command "try { Invoke-WebRequest -Uri 'https://raw.githubusercontent.com/DiegoAnyG/PoliScreen/main/scripts/PoliScreen.ico' -OutFile '%~dp0PoliScreen.ico' -UseBasicParsing } catch {}" >nul 2>&1
)
powershell -NoProfile -Command "$ws = New-Object -ComObject WScript.Shell; $desk = $ws.SpecialFolders.Item('Desktop'); $lnk = Join-Path $desk 'PoliScreen.lnk'; if (-not (Test-Path $lnk)) { $s = $ws.CreateShortcut($lnk); $s.TargetPath = '%~f0'; if (Test-Path '%~dp0PoliScreen.ico') { $s.IconLocation = '%~dp0PoliScreen.ico' }; $s.WorkingDirectory = $env:USERPROFILE; $s.Save() }" >nul 2>&1
echo.

set "CONFIG_FILE=%PROJECTS%\config.env"
set "TOOLS_DIR=%PROJECTS%\tools"
if not exist "%TOOLS_DIR%" mkdir "%TOOLS_DIR%"

set "DO_SETUP=0"
if /i "%~1"=="--setup" set "DO_SETUP=1"
if /i "%~1"=="-s" set "DO_SETUP=1"
if /i "%~1"=="--configure" set "DO_SETUP=1"
if not exist "%CONFIG_FILE%" set "DO_SETUP=1"

if "%DO_SETUP%"=="1" goto :run_questionnaire

echo   Config file         %CONFIG_FILE%
choice /C SC /T 3 /D C /M "  Press S to reconfigure components, or C to continue (3s)" >nul 2>&1
if not errorlevel 2 goto :run_questionnaire
goto :load_config

:run_questionnaire
set "GPU_NAME="
for /f "usebackq delims=" %%g in (`powershell -NoProfile -Command "try { $v = (Get-CimInstance Win32_VideoController | Where-Object { $_.Name -match 'NVIDIA' } | Select-Object -First 1).Name; if ($v) { Write-Output $v } else { Write-Output 'NONE' } } catch { Write-Output 'NONE' }"`) do (
    set "GPU_NAME=%%g"
)
set "HAS_GPU=0"
if defined GPU_NAME (
    if not "%GPU_NAME%"=="NONE" (
        set "HAS_GPU=1"
        set "GPU_STATUS=[GPU Detected: %GPU_NAME%]"
    ) else (
        set "GPU_STATUS=[WARNING: No NVIDIA GPU Detected]"
    )
) else (
    set "GPU_STATUS=[WARNING: No NVIDIA GPU Detected]"
)

echo.
echo   ===================================================================
echo    What do you want to install?
echo   ===================================================================
echo.
echo   Base (Default):
echo     [x] 1. PoliScreen
echo            Core virtual screening: Vina 1.2, RDKit, OpenBabel, PLIP, fpocket.
echo.
echo   Complements:
echo     [ ] 2. ADME-AI (admelab)
echo            Deep-learning pharmacokinetic and toxicity endpoint predictions.
echo.
echo     [ ] 3. ADCP (AutoDock CrankPep)
echo            Conformational sampling and flexible peptide docking (1-20 residues).
echo.
echo     [ ] 4. CAVER Suite
echo            Geometric detection and analysis of active-site access tunnels.
echo.
echo     [ ] 5. GNINA (CNN Scoring)
echo            Convolutional neural network scoring and docking.
echo            Status: %GPU_STATUS%
echo.
echo   -------------------------------------------------------------------
echo   Select complements to install (e.g. 2,4,5 or 'all') [Default: Base only]:
set "SELECTION="
set /p "SELECTION=> "

set "ENABLE_ADMET=0"
set "ENABLE_ADCP=0"
set "ENABLE_CAVER=0"
set "ENABLE_GNINA=0"

if /i "%SELECTION%"=="all" (
    set "ENABLE_ADMET=1"
    set "ENABLE_ADCP=1"
    set "ENABLE_CAVER=1"
    set "ENABLE_GNINA=1"
) else (
    echo %SELECTION% | find "2" >nul && set "ENABLE_ADMET=1"
    echo %SELECTION% | find "3" >nul && set "ENABLE_ADCP=1"
    echo %SELECTION% | find "4" >nul && set "ENABLE_CAVER=1"
    echo %SELECTION% | find "5" >nul && set "ENABLE_GNINA=1"
)

(
    echo # PoliScreen Component Configuration
    echo ENABLE_ADMET=%ENABLE_ADMET%
    echo ENABLE_ADCP=%ENABLE_ADCP%
    echo ENABLE_CAVER=%ENABLE_CAVER%
    echo ENABLE_GNINA=%ENABLE_GNINA%
) > "%CONFIG_FILE%"

echo   Configuration saved to %CONFIG_FILE%
echo.

:load_config
set "ENABLE_ADMET=0"
set "ENABLE_ADCP=0"
set "ENABLE_CAVER=0"
set "ENABLE_GNINA=0"
if exist "%CONFIG_FILE%" (
    for /f "usebackq tokens=1,* delims==" %%A in ("%CONFIG_FILE%") do (
        if "%%A"=="ENABLE_ADMET" set "ENABLE_ADMET=%%B"
        if "%%A"=="ENABLE_ADCP" set "ENABLE_ADCP=%%B"
        if "%%A"=="ENABLE_CAVER" set "ENABLE_CAVER=%%B"
        if "%%A"=="ENABLE_GNINA" set "ENABLE_GNINA=%%B"
    )
)

echo   Components          Base [OK]  ADME-AI [%ENABLE_ADMET%]  ADCP [%ENABLE_ADCP%]  CAVER [%ENABLE_CAVER%]  GNINA [%ENABLE_GNINA%]
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
    echo   Checking for updates...
    docker pull %IMAGE%
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

set "EXTRA_FLAGS="
if "%ENABLE_ADMET%"=="0" (
    set "EXTRA_FLAGS=%EXTRA_FLAGS% -e POLISCREEN_WITH_ADMET=0 -e POLISCREEN_ADME_PYTHON="
) else (
    set "EXTRA_FLAGS=%EXTRA_FLAGS% -e POLISCREEN_WITH_ADMET=1"
)

if "%ENABLE_CAVER%"=="0" (
    set "EXTRA_FLAGS=%EXTRA_FLAGS% -e POLISCREEN_WITH_CAVER=0 -e POLISCREEN_CAVER="
) else (
    set "EXTRA_FLAGS=%EXTRA_FLAGS% -e POLISCREEN_WITH_CAVER=1"
)

if "%ENABLE_GNINA%"=="1" (
    powershell -NoProfile -Command "if (Get-CimInstance Win32_VideoController | Where-Object { $_.Name -match 'NVIDIA' }) { exit 0 } else { exit 1 }" >nul 2>&1
    if not errorlevel 1 (
        set "EXTRA_FLAGS=%EXTRA_FLAGS% --gpus all"
    )
)

docker run --rm -it --init --name poliscreen -p 127.0.0.1:8501:8501 -v "%PROJECTS%:/data" -v "%TOOLS_DIR%:/root/poliscreen_tools" %THEME% %EXTRA_FLAGS% %IMAGE%

endlocal
