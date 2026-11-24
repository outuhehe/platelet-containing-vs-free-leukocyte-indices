<#
.SYNOPSIS
  Single entry point for the reproducible MIMIC-IV / eICU analysis pipeline.

.DESCRIPTION
  Runs the numbered R scripts in R\ in order (01..15, excluding 00_common.R).
  Each stage reads the previous stage's output; the run stops at the first failure.
  R routinely writes progress/warnings to stderr, so success is judged solely by
  the Rscript exit code (Start-Process isolates stderr and does not treat it as a
  PowerShell error). Database credentials come from config/db.env (or pgpass) via
  00_common.R; this orchestrator never connects to a database itself.

.PARAMETER SkipBootstrap
  Sets IPW_RUN_BOOTSTRAP=false so the IPW stage skips its 2000-replicate bootstrap.
  Use for a fast logic check; the published run should omit this switch.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File .\run_all.ps1
.EXAMPLE
  powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -SkipBootstrap
#>
[CmdletBinding()]
param(
  [switch]$SkipBootstrap
)

# Continue rather than Stop: native programs (Rscript) writing to stderr is normal.
$ErrorActionPreference = 'Continue'
$Root = $PSScriptRoot
$Rdir = Join-Path $Root 'R'
$LogDir = Join-Path $Root 'outputs\logs'
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

# Locate Rscript.exe: PATH first, then the standard per-version install folder.
$rsCmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
$Rscript = if ($rsCmd) { $rsCmd.Source } else { $null }
if (-not $Rscript) {
  $candidate = Get-ChildItem 'C:\Program Files\R' -Directory -ErrorAction SilentlyContinue |
    Sort-Object Name -Descending |
    ForEach-Object { Join-Path $_.FullName 'bin\Rscript.exe' } |
    Where-Object { Test-Path $_ } | Select-Object -First 1
  if (-not $candidate) { throw 'Rscript.exe not found on PATH or under C:\Program Files\R.' }
  $Rscript = $candidate
}
Write-Host "Using Rscript: $Rscript"

if ($SkipBootstrap) { $env:IPW_RUN_BOOTSTRAP = 'false' }
else { Remove-Item Env:IPW_RUN_BOOTSTRAP -ErrorAction SilentlyContinue }

# Numbered analysis scripts only; 00_common.R is sourced, not run.
$scripts = Get-ChildItem -Path $Rdir -Filter '*.R' |
  Where-Object { $_.Name -match '^(0[1-9]|1[0-9])_' } |
  Sort-Object Name

if (-not (Test-Path (Join-Path $Root 'config\db.env'))) {
  Write-Warning 'config\db.env not found; extraction stages will fall back to pgpass.conf.'
}

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$total = $scripts.Count
$i = 0
foreach ($s in $scripts) {
  $i++
  $base = [System.IO.Path]::GetFileNameWithoutExtension($s.Name)
  $outLog = Join-Path $LogDir "$base.txt"
  $errLog = Join-Path $LogDir "${base}_err.txt"
  Write-Host ("[{0}/{1}] {2} ... " -f $i, $total, $s.Name) -NoNewline
  $t0 = Get-Date
  # Start-Process captures both streams without promoting stderr to an error;
  # -PassThru lets us read the real exit code.
  $proc = Start-Process -FilePath $Rscript -ArgumentList @('"' + $s.FullName + '"') `
    -RedirectStandardOutput $outLog -RedirectStandardError $errLog `
    -WindowStyle Hidden -Wait -PassThru
  $code = $proc.ExitCode
  $secs = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)
  if ($code -ne 0) {
    Write-Host "FAILED (exit $code, ${secs}s)" -ForegroundColor Red
    Write-Host "See $errLog and $outLog" -ForegroundColor Yellow
    Get-Content $errLog -Tail 20
    exit $code
  }
  Write-Host "ok (${secs}s)" -ForegroundColor Green
}
$sw.Stop()
Write-Host ("`nPipeline complete: {0} stages in {1:N1} minutes." -f $total, $sw.Elapsed.TotalMinutes) -ForegroundColor Cyan
Write-Host "Tables: $(Join-Path $Root 'outputs\tables')"
Write-Host "Figures: $(Join-Path $Root 'outputs\figures')"
