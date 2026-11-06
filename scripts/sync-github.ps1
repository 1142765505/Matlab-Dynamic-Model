param([switch]$SkipQuietPeriod)

$ErrorActionPreference = 'Stop'
$repoPath = Split-Path -Parent $PSScriptRoot
$env:GIT_TERMINAL_PROMPT = '0'
$env:GCM_INTERACTIVE = 'Never'

function Invoke-Git {
    param([string[]]$GitArgs)
    $previousPreference = $ErrorActionPreference
    try {
        # Windows PowerShell treats native stderr as ErrorRecord objects,
        # even when Git succeeds (for example, fetch progress).
        $ErrorActionPreference = 'Continue'
        $result = & git -C $repoPath @GitArgs 2>&1
        $gitExitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $previousPreference
    }
    if ($gitExitCode -ne 0) { throw ($result -join "`n") }
    return $result
}

$lock = $null
try {
    $lockPath = Join-Path $repoPath '.git/github-sync.lock'
    try {
        $lock = [System.IO.File]::Open($lockPath, 'OpenOrCreate', 'ReadWrite', 'None')
    } catch [System.IO.IOException] {
        Write-Output 'Another sync is running; skipped.'
        exit 0
    }
    $remote = Invoke-Git -GitArgs @('remote', 'get-url', 'origin')
    if ($remote -ne 'https://github.com/1142765505/Matlab-Dynamic-Model.git') {
        throw 'Unexpected origin URL; sync stopped.'
    }
    $branch = Invoke-Git -GitArgs @('branch', '--show-current')
    if ($branch -ne 'main') { throw 'Switch back to main before syncing.' }
    foreach ($marker in @('MERGE_HEAD', 'CHERRY_PICK_HEAD', 'rebase-merge', 'rebase-apply')) {
        if (Test-Path -LiteralPath (Join-Path $repoPath ".git/$marker")) {
            throw 'An unfinished Git operation exists; sync stopped.'
        }
    }
    $files = @(Get-ChildItem -LiteralPath $repoPath -Recurse -Force -File |
        Where-Object { $_.FullName -notlike "$repoPath\.git\*" })
    if (-not $SkipQuietPeriod -and ($files | Where-Object { $_.LastWriteTimeUtc -gt [DateTime]::UtcNow.AddSeconds(-30) })) {
        Write-Output 'Files were recently saved; deferred until the next check.'
        exit 0
    }
    $largeFiles = @($files | Where-Object { $_.Length -ge 100MB })
    if ($largeFiles.Count) { throw ('Files exceed the GitHub size limit: ' + ($largeFiles.Name -join ', ')) }
    Invoke-Git -GitArgs @('add', '-A') | Out-Null
    & git -C $repoPath diff --cached --quiet
    $diffExit = $LASTEXITCODE
    if ($diffExit -eq 1) {
        $message = 'Auto-sync local changes ' + (Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')
        Invoke-Git -GitArgs @('commit', '-m', $message) | Write-Output
    } elseif ($diffExit -ne 0) {
        throw 'Unable to inspect staged changes.'
    }
    try {
        Invoke-Git -GitArgs @('fetch', 'origin', 'main') | Out-Null
    } catch {
        Write-Output 'Git fetch failed; trying the authenticated GitHub API.'
        & (Join-Path $PSScriptRoot 'push-github-api.ps1')
        return
    }
    & git -C $repoPath merge-base --is-ancestor origin/main HEAD
    if ($LASTEXITCODE -ne 0) {
        throw 'Remote main has new commits. Reconcile local and remote changes before retrying; no force push was attempted.'
    }
    try {
        Invoke-Git -GitArgs @('-c', 'credential.interactive=false', 'push', 'origin', 'main') | Write-Output
    } catch {
        Write-Output 'Git push failed; trying the authenticated GitHub API.'
        & (Join-Path $PSScriptRoot 'push-github-api.ps1')
        return
    }
    $localHead = Invoke-Git -GitArgs @('rev-parse', 'HEAD')
    $remoteHead = Invoke-Git -GitArgs @('ls-remote', 'origin', 'refs/heads/main')
    if (($remoteHead -split '\s+')[0] -ne $localHead) { throw 'Remote verification failed.' }
    Write-Output "Verified GitHub main: $localHead"
} finally {
    if ($null -ne $lock) { $lock.Dispose() }
}
