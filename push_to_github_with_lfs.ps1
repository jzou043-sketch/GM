$ErrorActionPreference = "Stop"
Set-Location "C:\Users\15280\OneDrive\Desktop\GM_github_upload"

git init
git lfs install
git remote remove origin 2>$null
git remote add origin https://github.com/jzou043-sketch/GM.git

git add .gitattributes README.md .
git commit -m "Upload GM project files"
git branch -M main
git push -u origin main
