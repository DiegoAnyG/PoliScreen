# Installs the official AutoDock Vina binary on Windows.
#
# Not via conda: bioconda only ships autodock-vina 1.1.2 (2011), whose scoring function
# differs from the 1.2.x series. PoliScreen was validated with 1.2.5, so the version is
# pinned explicitly and its SHA256 verified.
#
#   powershell -ExecutionPolicy Bypass -File scripts\get_vina.ps1
#   powershell -ExecutionPolicy Bypass -File scripts\get_vina.ps1 -Version 1.2.7
param(
    [string]$Version = "1.2.5",
    [string]$Dest = ""
)
$ErrorActionPreference = "Stop"

if ([string]::IsNullOrEmpty($Dest)) {
    if ($env:CONDA_PREFIX) { $Dest = Join-Path $env:CONDA_PREFIX "Scripts" }
    else { $Dest = Join-Path $env:LOCALAPPDATA "PoliScreen\bin" }
}

# Known sums of the verified versions: if the binary changed, installation aborts
# en vez de producir resultados distintos en silencio.
$Sums = @{ "1.2.5" = "b5a23b226ee5dc1b028df04082d13a7b0e5b9b12b988bbcf9498cb351e20b38f" }

$Url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v$Version/vina_${Version}_win.exe"
$Out = Join-Path $Dest "vina.exe"

New-Item -ItemType Directory -Force -Path $Dest | Out-Null
Write-Host "Downloading AutoDock Vina $Version -> $Out"
Invoke-WebRequest -Uri $Url -OutFile $Out

if ($Sums.ContainsKey($Version)) {
    $hash = (Get-FileHash -Algorithm SHA256 -Path $Out).Hash.ToLower()
    if ($hash -ne $Sums[$Version]) {
        Remove-Item $Out -Force
        throw "La suma SHA256 no coincide. Descarga corrupta o binario alterado; se aborta."
    }
    Write-Host "Integrity verified (SHA256)."
} else {
    Write-Warning "No hay suma conocida para $Version; no se verifico la integridad."
}

& $Out --version | Select-Object -First 1

if (-not (Get-Command obrms -ErrorAction SilentlyContinue)) {
    Write-Warning "No encuentro 'obrms' (OpenBabel). La estabilidad geometrica de la confianza quedara vacia."
}
Write-Host "Done. Check the installation with: poliscreen info"
