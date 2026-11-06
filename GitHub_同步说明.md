# GitHub 同步

本文件夹对应仓库：https://github.com/1142765505/Matlab-Dynamic-Model

Codex 自动任务每 5 分钟检查一次本文件夹。有修改时自动提交到 `main` 并上传；没有修改时不创建提交。需要电脑开机、Codex 应用保持运行且网络可用。离线期间的修改会在后续检查时重试上传。

每次检查会跳过最近 30 秒仍有文件写入的目录，避免在 MATLAB 连续输出时立即提交。仿真结果（包括 MAT、CSV 和图片）一并上传；MATLAB 自动恢复文件、临时文件及 `.env` 本地配置被忽略。请勿将密码、令牌等私密文件放入同步目录。

手动同步：在本目录的 PowerShell 中运行：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\sync-github.ps1
```

首次使用须在本机完成 GitHub 登录。若登录失效，可运行 `git credential-manager github login --username 1142765505 --device` 重新授权。

当普通 Git 连接失败时，脚本会使用 Git Credential Manager 已保存的凭据，经 GitHub API 上传同一批 Git 对象，逐一校验文件、目录树和提交的哈希值。凭据不会写入文件或日志。

如远端出现本机没有的提交，脚本会停止上传并保留本地文件及提交。需处理远端修改后再同步，脚本不会强制覆盖远端历史。

在 Codex 的自动任务页面可暂停或删除“MATLAB 文件夹自动同步 GitHub”。
