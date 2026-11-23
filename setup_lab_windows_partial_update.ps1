<#
.SYNOPSIS
    Lightweight, reproducible dependency setup for local .sld extraction on a Windows lab PC.
.DESCRIPTION
    Configures everything required to run collect_4d_sld_timelapses.py locally on Windows
    without Administrator privileges:
      1. Enforces TLS 1.2 and safe strict-mode state.
      2. Creates an isolated user-space directory structure.
      3. Detects Java from PATH/JAVA_HOME or an existing FIJI (Fiji.app) installation,
         falling back to a lightweight portable Eclipse Temurin JRE 17 (.zip) download.
      4. Downloads and unpacks Bio-Formats CLI tools (bftools 8.5.0), automatically
         upgrading older installed versions (such as 7.3.0) in place.
      5. Provisions a pinned Python 3.12.7 environment via portable uv (skipping heavy
         HPC-only Cellpose/PyTorch packages by default unless -IncludeHpcPackages is set).
      6. Generates activate_mito4d.ps1 and verifies both Python and Bio-Formats + Java.
.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\setup_lab_windows.ps1
    powershell -ExecutionPolicy Bypass -File .\setup_lab_windows.ps1 -FijiDir "C:\Users\BestmanLab\Fiji.app"
    powershell -ExecutionPolicy Bypass -File .\setup_lab_windows.ps1 -DryRun
#>

[CmdletBinding()]
param(
    [string]$InstallDir = "$env:LOCALAPPDATA\mito_pipeline_tools",
    [string]$FijiDir = "",
    [string]$EnvName = "mito4d_env",
    [string]$PythonVersion = "3.12.7",
    [string]$BfToolsVersion = "8.5.0",
    [switch]$IncludeHpcPackages,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$global:LASTEXITCODE = 0

function Write-Step {
    param(
        [int]$Number,
        [int]$Total,
        [string]$Title,
        [string]$Reason
    )
    Write-Host "`n[$Number/$Total] $Title" -ForegroundColor Cyan
    Write-Host "      Reason: $Reason" -ForegroundColor DarkGray
}

function Find-JavaHome {
    param([string]$ExplicitFijiDir, [string]$LocalJreRoot)

    # 1. Check existing environment variables first.
    foreach ($VarName in @("BF_JAVA_HOME", "JAVA_HOME")) {
        $Val = [Environment]::GetEnvironmentVariable($VarName)
        if ($Val -and (Test-Path (Join-Path $Val "bin\java.exe"))) {
            return (Resolve-Path $Val).Path
        }
    }

    # 2. Check if java.exe is already on the system PATH (excluding dead Windows system32/javapath stubs).
    $JavaCmd = Get-Command java.exe -ErrorAction SilentlyContinue
    if ($JavaCmd -and $JavaCmd.Source -notmatch "Common Files\\Oracle\\Java\\javapath|Windows\\System32") {
        & $JavaCmd.Source -version *> $null
        if ($LASTEXITCODE -eq 0) {
            return (Split-Path (Split-Path $JavaCmd.Source -Parent) -Parent)
        }
    }

    # 3. Search common FIJI (Fiji.app) locations for FIJI's bundled OpenJDK/JRE.
    $FijiCandidates = @()
    if ($ExplicitFijiDir) { $FijiCandidates += $ExplicitFijiDir }
    $FijiCandidates += @(
        "C:\Fiji.app",
        "D:\Fiji.app",
        (Join-Path $env:USERPROFILE "Fiji.app"),
        (Join-Path $env:USERPROFILE "Desktop\Fiji.app"),
        (Join-Path $env:USERPROFILE "Downloads\Fiji.app"),
        (Join-Path $env:LOCALAPPDATA "Fiji.app"),
        "C:\Program Files\Fiji.app"
    )
    foreach ($Root in $FijiCandidates) {
        $FijiJavaDir = Join-Path $Root "java"
        if (Test-Path $FijiJavaDir) {
            $FoundJava = Get-ChildItem -Path $FijiJavaDir -Filter "java.exe" -Recurse -ErrorAction SilentlyContinue |
                Where-Object { $_.Directory.Name -ieq "bin" } |
                Select-Object -First 1
            if ($FoundJava) {
                return $FoundJava.Directory.Parent.FullName
            }
        }
    }

    # 4. Check if a portable JRE was already downloaded into our tool folder in a prior run.
    if (Test-Path $LocalJreRoot) {
        $LocalJava = Get-ChildItem -Path $LocalJreRoot -Filter "java.exe" -Recurse -ErrorAction SilentlyContinue |
            Where-Object { $_.Directory.Name -ieq "bin" } |
            Select-Object -First 1
        if ($LocalJava) {
            return $LocalJava.Directory.Parent.FullName
        }
    }

    return $null
}

$TotalSteps = 6
$EnvPrefix = Join-Path $InstallDir "envs\$EnvName"
$CacheRoot = Join-Path $InstallDir ".cache"
$TmpDir = Join-Path $CacheRoot "tmp"
$UvDir = Join-Path $InstallDir "uv"
$JreRoot = Join-Path $InstallDir "jre"
$BfToolsDir = Join-Path $InstallDir "bftools"
$ActivateScript = Join-Path $InstallDir "activate_mito4d.ps1"

Write-Host "======================================================================" -ForegroundColor White
Write-Host "Windows Lab PC Setup - Local .sld Extraction (4D Mitochondrial Pipeline)" -ForegroundColor White
Write-Host "Install Dir         : $InstallDir"
Write-Host "Python Environment  : $EnvPrefix (Python $PythonVersion)"
Write-Host "Bio-Formats Version : $BfToolsVersion"
Write-Host "HPC Analysis Pkgs   : $IncludeHpcPackages"
Write-Host "Dry Run             : $DryRun"
Write-Host "======================================================================" -ForegroundColor White

if ($DryRun) {
    $PreviewJava = Find-JavaHome -ExplicitFijiDir $FijiDir -LocalJreRoot $JreRoot
    Write-Host "`n[DRY-RUN] No files will be downloaded or modified." -ForegroundColor Yellow
    if ($PreviewJava) {
        Write-Host "[DRY-RUN] Detected existing Java runtime at: $PreviewJava" -ForegroundColor Green
    } else {
        Write-Host "[DRY-RUN] No existing/FIJI Java found; would download portable Temurin JRE 17 (~42 MB) to: $JreRoot"
    }
    Write-Host "[DRY-RUN] Would install/upgrade Bio-Formats $BfToolsVersion CLI (showinf.bat / bfconvert.bat) to: $BfToolsDir"
    Write-Host "[DRY-RUN] Would install portable uv and create Python $PythonVersion venv at: $EnvPrefix"
    Write-Host "[DRY-RUN] Would generate session activation script: $ActivateScript"
    exit 0
}

# ------------------------------------------------------------------------------
# Step 1: Enforce TLS 1.2 and initialize safe PowerShell session state.
# ------------------------------------------------------------------------------
Write-Step -Number 1 -Total $TotalSteps -Title "Configure TLS 1.2 Security Protocol" `
    -Reason "Windows PowerShell 5.1 defaults to older TLS versions that modern servers (Adoptium, OME, GitHub) reject."
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

# ------------------------------------------------------------------------------
# Step 2: Create isolated user-space directories.
# ------------------------------------------------------------------------------
Write-Step -Number 2 -Total $TotalSteps -Title "Create User-Space Workspace and Cache Directories" `
    -Reason "Keeping tools and caches inside a single user-writable folder avoids Administrator prompts and prevents cluttering system paths."
foreach ($Dir in @($InstallDir, (Join-Path $InstallDir "envs"), (Join-Path $CacheRoot "uv"), (Join-Path $CacheRoot "pip"), $TmpDir)) {
    if (-not (Test-Path $Dir)) {
        New-Item -ItemType Directory -Path $Dir -Force | Out-Null
    }
}
$env:UV_CACHE_DIR = Join-Path $CacheRoot "uv"
$env:PIP_CACHE_DIR = Join-Path $CacheRoot "pip"

# ------------------------------------------------------------------------------
# Step 3: Locate Java (System or FIJI) or download lightweight portable JRE 17.
# ------------------------------------------------------------------------------
Write-Step -Number 3 -Total $TotalSteps -Title "Locate Existing Java (FIJI/System) or Download Portable JRE 17" `
    -Reason "Bio-Formats $BfToolsVersion requires Java 8+, so we reuse FIJI's bundled Java if present or fetch a minimal ~42 MB portable JRE."
$JavaHome = Find-JavaHome -ExplicitFijiDir $FijiDir -LocalJreRoot $JreRoot
if ($JavaHome) {
    Write-Host "      Reusing existing Java runtime at: $JavaHome" -ForegroundColor Green
} else {
    Write-Host "      No system or FIJI Java found; downloading lightweight Eclipse Temurin JRE 17 (.zip)..."
    $JreZip = Join-Path $TmpDir "temurin_jre17_win64.zip"
    $JreUrl = "https://api.adoptium.net/v3/binary/latest/17/ga/windows/x64/jre/hotspot/normal/eclipse"
    Invoke-WebRequest -Uri $JreUrl -OutFile $JreZip -UseBasicParsing
    if (Test-Path $JreRoot) { Remove-Item $JreRoot -Recurse -Force }
    Expand-Archive -Path $JreZip -DestinationPath $JreRoot -Force
    Remove-Item $JreZip -Force

    $JavaHome = Find-JavaHome -ExplicitFijiDir $FijiDir -LocalJreRoot $JreRoot
    if (-not $JavaHome) {
        throw "Downloaded JRE archive to $JreRoot, but bin\java.exe could not be located."
    }
    Write-Host "      Installed portable JRE 17 at: $JavaHome" -ForegroundColor Green
}
$env:JAVA_HOME = $JavaHome
$env:BF_JAVA_HOME = $JavaHome
$env:PATH = "$(Join-Path $JavaHome 'bin');$env:PATH"

# ------------------------------------------------------------------------------
# Step 4: Download or upgrade Bio-Formats CLI tools (bftools 8.5.0).
# ------------------------------------------------------------------------------
Write-Step -Number 4 -Total $TotalSteps -Title "Install or Upgrade Bio-Formats CLI Tools (bftools $BfToolsVersion)" `
    -Reason "bftools $BfToolsVersion adds lenient metadata parsing and SlideBook 7 fixes over 7.3.0 while keeping identical showinf/bfconvert CLI syntax."
$ShowInfBat = Join-Path $BfToolsDir "showinf.bat"
$BfConvertBat = Join-Path $BfToolsDir "bfconvert.bat"
$NeedsBfInstall = $true
if ((Test-Path $ShowInfBat) -and (Test-Path $BfConvertBat)) {
    $CurrentBfVer = (& $ShowInfBat -version 2>&1 | Out-String)
    if ($CurrentBfVer -match [regex]::Escape($BfToolsVersion)) {
        $NeedsBfInstall = $false
        Write-Host "      bftools $BfToolsVersion is already installed at: $BfToolsDir"
    } else {
        Write-Host "      Found older bftools installation; upgrading in place to $BfToolsVersion..." -ForegroundColor Yellow
        Remove-Item $BfToolsDir -Recurse -Force
    }
}
if ($NeedsBfInstall) {
    $BfZip = Join-Path $TmpDir "bftools_$BfToolsVersion.zip"
    Invoke-WebRequest -Uri "https://downloads.openmicroscopy.org/bio-formats/$BfToolsVersion/artifacts/bftools.zip" -OutFile $BfZip -UseBasicParsing
    Expand-Archive -Path $BfZip -DestinationPath $InstallDir -Force
    Remove-Item $BfZip -Force
}
$env:PATH = "$BfToolsDir;$env:PATH"

# ------------------------------------------------------------------------------
# Step 5: Provision reproducible Python 3.12.7 environment via portable uv.
# ------------------------------------------------------------------------------
Write-Step -Number 5 -Total $TotalSteps -Title "Provision Python $PythonVersion Environment via Portable uv" `
    -Reason "collect_4d_sld_timelapses.py runs on Python's standard library, so creating a lightweight Python $PythonVersion venv avoids unnecessary multi-gigabyte HPC package downloads."
$UvExe = $null
$ExistingUv = Get-Command uv.exe -ErrorAction SilentlyContinue
if ($ExistingUv) {
    $UvExe = $ExistingUv.Source
} else {
    $PortableUv = Join-Path $UvDir "uv.exe"
    if (-not (Test-Path $PortableUv)) {
        $UvZip = Join-Path $TmpDir "uv-x86_64-pc-windows-msvc.zip"
        Invoke-WebRequest -Uri "https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-pc-windows-msvc.zip" -OutFile $UvZip -UseBasicParsing
        if (-not (Test-Path $UvDir)) { New-Item -ItemType Directory -Path $UvDir -Force | Out-Null }
        Expand-Archive -Path $UvZip -DestinationPath $UvDir -Force
        Remove-Item $UvZip -Force
    }
    $UvExe = $PortableUv
}

$VenvPython = Join-Path $EnvPrefix "Scripts\python.exe"
if (-not (Test-Path $VenvPython)) {
    & $UvExe venv $EnvPrefix --python $PythonVersion
    if ($LASTEXITCODE -ne 0) {
        throw "uv failed to create the Python $PythonVersion environment at $EnvPrefix."
    }
} else {
    Write-Host "      Python environment already exists at: $EnvPrefix"
}

if ($IncludeHpcPackages) {
    Write-Host "      Installing optional HPC Stage 3-5 analysis packages (-IncludeHpcPackages)..."
    $HpcPackages = @(
        "numpy>=1.26.0", "scipy>=1.12.0", "pandas>=2.2.0", "pyyaml>=6.0.1",
        "tifffile>=2024.2.12", "imagecodecs>=2024.1.1", "pillow>=10.2.0",
        "opencv-python-headless>=4.9.0", "scikit-image>=0.22.0",
        "scikit-fmm>=2023.4.2", "networkx>=3.2.1", "cellpose>=3.0.0"
    )
    & $UvExe pip install --python $VenvPython @HpcPackages
    if ($LASTEXITCODE -ne 0) { throw "Optional HPC package installation failed." }
}

# ------------------------------------------------------------------------------
# Step 6: Write activation script and verify Python + Java + Bio-Formats.
# ------------------------------------------------------------------------------
Write-Step -Number 6 -Total $TotalSteps -Title "Generate Activation Script & Verify Local Extraction Stack" `
    -Reason "A single activation script sets JAVA_HOME, BF_JAVA_HOME, and PATH per session without modifying the Windows Registry, and testing showinf.bat confirms readiness."
$ActivateContent = @"
# Source this file in PowerShell before running collect_4d_sld_timelapses.py:
#   . "$ActivateScript"
& "$EnvPrefix\Scripts\Activate.ps1"
`$env:JAVA_HOME = "$JavaHome"
`$env:BF_JAVA_HOME = "$JavaHome"
`$env:PATH = "$(Join-Path $JavaHome 'bin');$BfToolsDir;" + `$env:PATH
Write-Host "[OK] Activated $EnvName | Java: $JavaHome | bftools $BfToolsVersion" -ForegroundColor Green
"@
$ActivateContent | Out-File -FilePath $ActivateScript -Encoding utf8

$PyCheck = & $VenvPython -c "import sys, xml.etree.ElementTree, concurrent.futures, hashlib; print('PYTHON_OK:' + sys.version.split()[0])"
if ($LASTEXITCODE -ne 0 -or ($PyCheck -join " ") -notmatch "PYTHON_OK:") {
    throw "Python verification failed: $PyCheck"
}

$BfCheck = & $ShowInfBat -version 2>&1 | Out-String
if ($LASTEXITCODE -ne 0 -or $BfCheck -notmatch [regex]::Escape($BfToolsVersion)) {
    throw "Bio-Formats $BfToolsVersion + Java verification failed:`n$BfCheck"
}

Write-Host "`n======================================================================" -ForegroundColor Green
Write-Host "[SUCCESS] Local .sld extraction environment verified!" -ForegroundColor Green
Write-Host "  Python      : $($PyCheck -join ' ')"
Write-Host "  Java Home   : $JavaHome"
Write-Host "  Bio-Formats : bftools $BfToolsVersion ($BfToolsDir)"
Write-Host "`nTo activate in PowerShell and run local extraction:"
Write-Host "  . `"$ActivateScript`"" -ForegroundColor Yellow
Write-Host "  python collect_4d_sld_timelapses.py `"D:\RawSLD`" `"D:\Staged_OMETIFF`" --mode scan-only --workers 1" -ForegroundColor Yellow
Write-Host "======================================================================" -ForegroundColor Green
