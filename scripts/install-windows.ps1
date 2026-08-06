<#
.SYNOPSIS
    Prepares Windows and starts PoliScreen in Docker.

.DESCRIPTION
    Run it as many times as needed: it checks the state of the machine on every run and does
    the next thing that is missing, so after a reboot you simply run it again and it carries
    on. There is no saved state and nothing is scheduled to run by itself.

    Two things it cannot do, and will tell you about instead:

      * Turn on virtualization (VT-x / AMD-SVM). That lives in the BIOS/UEFI firmware and no
        program running under Windows can change it.
      * Skip the reboot Windows requires after enabling WSL. It stops and asks you to reboot.

.PARAMETER NoAdmet
    Leave the ADMET engine out. It is built in by default because without it the analogue
    builder and the ADMET report are disabled, which is the omission users notice first. It
    adds about 1.5 GB, using the CPU build of PyTorch.

.PARAMETER WithAdcp
    Also build peptide docking (ADCP) into the image. Adds about 900 MB and accepts the Scripps
    academic licence on your behalf, on your own machine.

.PARAMETER WithGnina
    Also build the neural-network rescoring (gnina) into the image. Adds about 4.5 GB and needs
    an NVIDIA GPU to actually run.

.PARAMETER Path
    Where to clone PoliScreen when the script is run on its own. Ignored when it is run from
    inside a checkout.

.EXAMPLE
    .\install-windows.ps1
.EXAMPLE
    .\install-windows.ps1 -WithAdcp
#>
[CmdletBinding()]
param(
    [switch]$NoAdmet,
    [switch]$WithAdcp,
    [switch]$WithGnina,
    [string]$Path = "$env:USERPROFILE\PoliScreen"
)

$ErrorActionPreference = 'Stop'
$RepoUrl = 'https://github.com/DiegoAnyG/PoliScreen.git'
$Port = 8501

function Write-Step { param($m) Write-Host "`n==> $m" -ForegroundColor Cyan }
function Write-Ok   { param($m) Write-Host "    OK   $m" -ForegroundColor Green }
function Write-Info { param($m) Write-Host "    ..   $m" -ForegroundColor Gray }
function Write-Warn { param($m) Write-Host "    !!   $m" -ForegroundColor Yellow }

function Stop-WithReason {
    param([string]$Title, [string[]]$Lines)
    Write-Host ""
    Write-Host "  $Title" -ForegroundColor Yellow
    Write-Host ""
    foreach ($l in $Lines) { Write-Host "    $l" }
    Write-Host ""
    exit 1
}

function Test-Admin {
    $id = [Security.Principal.WindowsIdentity]::GetCurrent()
    $p = New-Object Security.Principal.WindowsPrincipal($id)
    return $p.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Test-PortOpen {
    param([int]$Number)
    $c = New-Object Net.Sockets.TcpClient
    try {
        $c.Connect('127.0.0.1', $Number)
        return $c.Connected
    } catch {
        return $false
    } finally {
        $c.Close()
    }
}

Write-Host ""
Write-Host "  PoliScreen installer for Windows" -ForegroundColor White
Write-Host "  Run this again after any reboot; it continues where it stopped."

# ---------------------------------------------------------------- 1. privileges

Write-Step "Checking privileges"
if (-not (Test-Admin)) {
    Stop-WithReason "This script needs to run as Administrator." @(
        "Enabling Windows features and installing Docker Desktop both require it.",
        "",
        "Open the Start menu, type 'PowerShell', right-click it and choose",
        "'Run as administrator'. Then run this script again:",
        "",
        "    $PSCommandPath"
    )
}
Write-Ok "running as Administrator"

# ---------------------------------------------------------------- 2. Windows build

Write-Step "Checking the Windows version"
$build = [int](Get-CimInstance Win32_OperatingSystem).BuildNumber
if ($build -lt 19041) {
    Stop-WithReason "Windows is too old for WSL 2 (build $build, 19041 or newer needed)." @(
        "Update Windows from Settings > Windows Update, then run this script again."
    )
}
Write-Ok "build $build"

# ---------------------------------------------------------------- 3. virtualization

Write-Step "Checking hardware virtualization"
# When a hypervisor is already running, the firmware flag reads False because Hyper-V has
# claimed it, so a running hypervisor is the stronger signal and is checked first.
$hyperv = (Get-CimInstance Win32_ComputerSystem).HypervisorPresent
if ($hyperv) {
    Write-Ok "a hypervisor is already running, so virtualization is enabled"
} else {
    $virt = (Get-CimInstance Win32_Processor | Select-Object -First 1).VirtualizationFirmwareEnabled
    if (-not $virt) {
        Stop-WithReason "Virtualization is disabled in the BIOS/UEFI, and no program can enable it." @(
            "Docker cannot run without it. You have to change it in the firmware:",
            "",
            "  1. Reboot and press the setup key while the manufacturer logo shows.",
            "     Usually F2, F10, Del or Esc; some laptops use F1.",
            "  2. Look for 'Intel Virtualization Technology', 'Intel VT-x', 'AMD-V'",
            "     or 'SVM Mode'. It is often under Advanced, CPU Configuration or Security.",
            "  3. Set it to Enabled, save and exit.",
            "",
            "To confirm afterwards: Task Manager > Performance > CPU, where it should",
            "read 'Virtualization: Enabled'. Then run this script again.",
            "",
            "If you cannot change the firmware (a locked corporate machine, for instance),",
            "Docker is not an option on it. Install with conda instead, see docs/INSTALL.md."
        )
    }
    Write-Ok "enabled in firmware"
}

# ---------------------------------------------------------------- 4. Windows features

Write-Step "Checking the Windows features WSL 2 needs"
$features = @('Microsoft-Windows-Subsystem-Linux', 'VirtualMachinePlatform')
$enabledNow = @()

foreach ($f in $features) {
    $state = (Get-WindowsOptionalFeature -Online -FeatureName $f).State
    if ($state -eq 'Enabled') {
        Write-Ok "$f"
    } else {
        Write-Info "enabling $f"
        Enable-WindowsOptionalFeature -Online -FeatureName $f -All -NoRestart | Out-Null
        $enabledNow += $f
    }
}

if ($enabledNow.Count -gt 0) {
    Stop-WithReason "Windows needs a reboot before those features work." @(
        "Just enabled: $($enabledNow -join ', ')",
        "",
        "Reboot, then run this script again. It will pick up from here."
    )
}

# ---------------------------------------------------------------- 5. WSL

Write-Step "Checking WSL"
$wsl = Get-Command wsl.exe -ErrorAction SilentlyContinue
if (-not $wsl) {
    Stop-WithReason "wsl.exe is missing even though the features are enabled." @(
        "This usually means the reboot is still pending. Reboot and run this script again."
    )
}

# The kernel ships separately from the feature and is what Docker asks to update.
Write-Info "updating the WSL kernel (no effect if it is current)"
& wsl.exe --update 2>&1 | ForEach-Object { Write-Info $_ }
& wsl.exe --set-default-version 2 2>&1 | Out-Null
Write-Ok "WSL 2 is the default version"

# ---------------------------------------------------------------- 6. Docker Desktop

Write-Step "Checking Docker Desktop"
$docker = Get-Command docker.exe -ErrorAction SilentlyContinue
if (-not $docker) {
    $winget = Get-Command winget.exe -ErrorAction SilentlyContinue
    if (-not $winget) {
        Stop-WithReason "Docker Desktop is not installed and winget is unavailable to install it." @(
            "Download it by hand from https://www.docker.com/products/docker-desktop/",
            "install it, then run this script again."
        )
    }
    Write-Info "installing Docker Desktop, this takes a few minutes"
    & winget.exe install -e --id Docker.DockerDesktop `
        --accept-package-agreements --accept-source-agreements
    Stop-WithReason "Docker Desktop was installed and Windows needs a reboot." @(
        "Reboot, let Docker Desktop finish its first start, then run this script again."
    )
}
Write-Ok "docker.exe found"

Write-Step "Waiting for the Docker engine"
& docker.exe info 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    # Docker Desktop installs per-machine or per-user depending on how it was installed, and
    # winget gives the per-user layout, so both have to be looked at.
    $exe = @(
        "$env:ProgramFiles\Docker\Docker\Docker Desktop.exe",
        "$env:LOCALAPPDATA\Programs\DockerDesktop\Docker Desktop.exe"
    ) | Where-Object { Test-Path $_ } | Select-Object -First 1

    if ($exe) {
        Write-Info "starting Docker Desktop"
        Start-Process $exe | Out-Null
    } else {
        Write-Warn "Docker Desktop is installed but its launcher was not found; start it by hand"
    }
    $waited = 0
    while ($waited -lt 180) {
        Start-Sleep -Seconds 5
        $waited += 5
        & docker.exe info 2>&1 | Out-Null
        if ($LASTEXITCODE -eq 0) { break }
        Write-Info "still starting... ${waited}s"
    }
}

& docker.exe info 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Stop-WithReason "The Docker engine did not come up." @(
        "Open Docker Desktop and read what it reports; the error here is only a symptom.",
        "If it says 'Virtualization support not detected', see step 3 above.",
        "Once its whale icon stops animating, run this script again."
    )
}
Write-Ok "engine responding"

# ---------------------------------------------------------------- 7. the sources

Write-Step "Locating the PoliScreen sources"
# When run from inside a checkout, that checkout is used; otherwise one is cloned.
$here = Split-Path -Parent $PSCommandPath
$repo = Split-Path -Parent $here

if (Test-Path (Join-Path $repo 'docker\docker-compose.yml')) {
    Write-Ok "using this checkout: $repo"
} elseif (Test-Path (Join-Path $Path 'docker\docker-compose.yml')) {
    $repo = $Path
    Write-Ok "using the existing checkout: $repo"
} else {
    $git = Get-Command git.exe -ErrorAction SilentlyContinue
    if (-not $git) {
        Stop-WithReason "git is needed to download PoliScreen and is not installed." @(
            "Install it with:  winget install -e --id Git.Git",
            "then run this script again."
        )
    }
    Write-Info "cloning into $Path"
    & git.exe clone $RepoUrl $Path
    $repo = $Path
    Write-Ok "cloned"
}

# ---------------------------------------------------------------- 8. build and run

Write-Step "Building and starting the container"
$compose = @('-f', (Join-Path $repo 'docker\docker-compose.yml'))

# The optional engines are Linux-only, which is why they can go in the container even though
# they can never go in the native Windows installer.
$env:POLISCREEN_WITH_ADMET = if ($NoAdmet) { '0' } else { '1' }
$env:POLISCREEN_WITH_ADCP = if ($WithAdcp) { '1' } else { '0' }
$env:POLISCREEN_WITH_GNINA = if ($WithGnina) { '1' } else { '0' }

if ($NoAdmet) {
    Write-Warn "leaving the ADMET engine out: no analogue builder and no ADMET report"
} else {
    Write-Info "including the ADMET engine (about 1.5 GB, PyTorch CPU build)"
}
if ($WithAdcp)  { Write-Info "including ADCP (about 900 MB, accepts the Scripps academic licence)" }
if ($WithGnina) {
    Write-Info "including gnina (about 4.5 GB)"
    $gpu = Join-Path $repo 'docker\docker-compose.gpu.yml'
    if (Test-Path $gpu) { $compose += @('-f', $gpu) }
    Write-Warn "gnina needs an NVIDIA GPU exposed to Docker; without one it installs but will not run"
}

Push-Location $repo
try {
    & docker.exe compose @compose up --build -d
    if ($LASTEXITCODE -ne 0) {
        Stop-WithReason "The build failed." @(
            "The output above says why. The usual causes are in docs/INSTALL.md,",
            "under 'Common problems'."
        )
    }
} finally {
    Pop-Location
}

Write-Step "Waiting for the interface"
$waited = 0
while ($waited -lt 120) {
    if (Test-PortOpen -Number $Port) { break }
    Start-Sleep -Seconds 3
    $waited += 3
}

if (-not (Test-PortOpen -Number $Port)) {
    Stop-WithReason "The container is up but nothing is listening on port $Port." @(
        "Look at the log with:",
        "    docker compose -f `"$repo\docker\docker-compose.yml`" logs"
    )
}

Write-Ok "listening on port $Port"
Start-Process "http://localhost:$Port"

Write-Host ""
Write-Host "  PoliScreen is running at http://localhost:$Port" -ForegroundColor Green
Write-Host "  Projects are saved in $repo\proyectos"
Write-Host ""
Write-Host "  Stop it with:"
Write-Host "      docker compose -f `"$repo\docker\docker-compose.yml`" down"
Write-Host ""
