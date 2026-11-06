# Fallback for networks where Git smart HTTP is unavailable but api.github.com works.
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Net.Http
$repoPath = Split-Path -Parent $PSScriptRoot
$apiRoot = 'https://api.github.com/repos/1142765505/Matlab-Dynamic-Model'

function Read-GitObject {
    param([string]$Arguments, [switch]$Binary)
    $start = New-Object System.Diagnostics.ProcessStartInfo
    $start.FileName = 'git.exe'
    $start.Arguments = $Arguments
    $start.WorkingDirectory = $repoPath
    $start.UseShellExecute = $false
    $start.CreateNoWindow = $true
    $start.RedirectStandardOutput = $true
    $start.RedirectStandardError = $true
    $start.StandardOutputEncoding = [Text.Encoding]::UTF8
    $process = [Diagnostics.Process]::Start($start)
    if ($Binary) {
        $stream = New-Object IO.MemoryStream
        $process.StandardOutput.BaseStream.CopyTo($stream)
        $data = $stream.ToArray()
        $stream.Dispose()
    } else { $data = $process.StandardOutput.ReadToEnd() }
    $errorText = $process.StandardError.ReadToEnd()
    $process.WaitForExit()
    if ($process.ExitCode -ne 0) { throw $errorText }
    $process.Dispose()
    if ($Binary) { return ,$data }
    return $data
}

function Invoke-GitHubApi {
    param([string]$Method, [string]$Path, $Body = $null)
    $request = New-Object System.Net.Http.HttpRequestMessage
    $request.Method = New-Object System.Net.Http.HttpMethod($Method)
    $request.RequestUri = [uri]($apiRoot + $Path)
    if ($null -ne $Body) {
        $json = ConvertTo-Json -InputObject $Body -Depth 20 -Compress
        $request.Content = New-Object System.Net.Http.StringContent($json, [Text.Encoding]::UTF8, 'application/json')
    }
    try {
        $response = $client.SendAsync($request).GetAwaiter().GetResult()
        $responseText = $response.Content.ReadAsStringAsync().GetAwaiter().GetResult()
        if (-not $response.IsSuccessStatusCode) {
            throw ('GitHub API HTTP ' + [int]$response.StatusCode + ': ' + $responseText)
        }
        return ($responseText | ConvertFrom-Json)
    } finally {
        if ($null -ne $response) { $response.Dispose() }
        $request.Dispose()
    }
}

$client = New-Object System.Net.Http.HttpClient
$client.Timeout = [TimeSpan]::FromMinutes(5)
try {
    $credentialOutput = @('protocol=https', 'host=github.com', '') | git -c credential.interactive=false credential fill
    if ($LASTEXITCODE -ne 0) { throw 'GitHub login is required.' }
    $token = ($credentialOutput | Where-Object { $_ -like 'password=*' }) -replace '^password=', ''
    if (-not $token) { throw 'GitHub login is required.' }
    $client.DefaultRequestHeaders.Authorization = New-Object System.Net.Http.Headers.AuthenticationHeaderValue('Bearer', $token)
    $client.DefaultRequestHeaders.UserAgent.ParseAdd('MATLAB-GitHub-Sync/1.0')
    $token = $null
    $credentialOutput = $null
    $remote = Invoke-GitHubApi -Method GET -Path '/git/ref/heads/main'
    $remoteSha = $remote.object.sha
    $localSha = (Read-GitObject 'rev-parse HEAD').Trim()
    & git -C $repoPath merge-base --is-ancestor $remoteSha $localSha
    if ($LASTEXITCODE -ne 0) { throw 'Remote main has new commits; reconcile the histories before syncing.' }
    $commits = @((Read-GitObject "rev-list --reverse ${remoteSha}..${localSha}").Trim() -split "`n" | Where-Object { $_ })
    $uploaded = @{}
    foreach ($entry in ((Read-GitObject "ls-tree -r $remoteSha") -split "`n")) {
        if ($entry -match '^\d+ blob ([0-9a-f]{40})') { $uploaded[$Matches[1]] = $true }
    }
    foreach ($commitSha in $commits) {
        $entries = @()
        foreach ($line in ((Read-GitObject "ls-tree -r -z $commitSha") -split "`0")) {
            if (-not $line) { continue }
            if ($line -notmatch '^(\d+) blob ([0-9a-f]{40})\t(.+)$') { throw 'Unsupported Git tree entry.' }
            $mode, $blobSha, $path = $Matches[1], $Matches[2], $Matches[3]
            if (-not $uploaded.ContainsKey($blobSha)) {
                $bytes = Read-GitObject "cat-file blob $blobSha" -Binary
                $blob = Invoke-GitHubApi -Method POST -Path '/git/blobs' -Body @{encoding='base64';content=[Convert]::ToBase64String($bytes)}
                if ($blob.sha -ne $blobSha) { throw 'Uploaded blob hash mismatch.' }
                $uploaded[$blobSha] = $true
                Write-Output "Uploaded: $path"
            }
            $entries += @{path=$path; mode=$mode; type='blob'; sha=$blobSha}
        }
        $tree = Invoke-GitHubApi -Method POST -Path '/git/trees' -Body @{tree=$entries}
        $expectedTree = (Read-GitObject "show -s --format=%T $commitSha").Trim()
        if ($tree.sha -ne $expectedTree) { throw 'Uploaded tree hash mismatch.' }
        $rawCommit = Read-GitObject "cat-file commit $commitSha"
        $separator = $rawCommit.IndexOf("`n`n")
        $message = $rawCommit.Substring($separator + 2)
        $parents = @([regex]::Matches($rawCommit.Substring(0,$separator), '(?m)^parent ([0-9a-f]{40})$') | ForEach-Object { $_.Groups[1].Value })
        $fields = (Read-GitObject "show -s --format=%an%x00%ae%x00%aI%x00%cn%x00%ce%x00%cI $commitSha").TrimEnd("`r","`n") -split "`0"
        $body = @{message=$message; tree=$tree.sha; parents=$parents; author=@{name=$fields[0];email=$fields[1];date=$fields[2]}; committer=@{name=$fields[3];email=$fields[4];date=$fields[5]}}
        $commit = Invoke-GitHubApi -Method POST -Path '/git/commits' -Body $body
        if ($commit.sha -ne $commitSha) { throw 'Uploaded commit hash mismatch; remote branch was not changed.' }
    }
    if ($commits.Count) {
        $latest = Invoke-GitHubApi -Method GET -Path '/git/ref/heads/main'
        if ($latest.object.sha -ne $remoteSha) { throw 'Remote changed during upload; retry after reconciling.' }
        Invoke-GitHubApi -Method PATCH -Path '/git/refs/heads/main' -Body @{sha=$localSha;force=$false} | Out-Null
    }
    $verified = Invoke-GitHubApi -Method GET -Path '/git/ref/heads/main'
    if ($verified.object.sha -ne $localSha) { throw 'Remote verification failed.' }
    & git -C $repoPath update-ref refs/remotes/origin/main $localSha
    if ($LASTEXITCODE -ne 0) { throw 'Unable to update local remote-tracking ref.' }
    Write-Output "Verified GitHub main via API: $localSha"
} finally {
    $client.Dispose()
    $token = $null
    $credentialOutput = $null
}
