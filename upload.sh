#!/bin/bash
# Carica la cartella pipeline/ su GitHub
# uso:  ./upload.sh <nome_repo> ["messaggio commit"]
# (il repo deve gia esistere su GitHub, anche vuoto)

set -e
cd "$(dirname "$0")"

GITHUB_REPO_URL="$1"
COMMIT_MESSAGE="${2:-Update pipeline $(date +%Y-%m-%d)}"

if [ -z "$GITHUB_REPO_URL" ]; then
    echo "uso: ./upload.sh <nome_repo | https://github.com/manciounipd/<nome_repo>.git> [\"messaggio\"]"
    exit 1
fi

# solo il nome del repo -> costruisco l'indirizzo completo
if [[ "$GITHUB_REPO_URL" != *:* ]]; then
    GITHUB_REPO_URL="https://github.com/manciounipd/${GITHUB_REPO_URL%.git}.git"
fi

# init repo se non esiste
if [ ! -d ".git" ]; then
    git init
fi

git config user.email "enrico.mancin@unipd.it"
git config user.name "manciounipd"

# remote origin (aggiunge o aggiorna)
if git remote | grep -q "^origin$"; then
    git remote set-url origin "$GITHUB_REPO_URL"
else
    git remote add origin "$GITHUB_REPO_URL"
fi

git add .
if git diff --cached --quiet; then
    echo "niente di nuovo da committare"
else
    git commit -m "$COMMIT_MESSAGE"
fi
git branch -M main

# se il repo remoto ha gia un main, unisco le storie prima del push
if git ls-remote --exit-code --heads origin main >/dev/null 2>&1; then
    git pull --no-rebase --no-edit origin main --allow-unrelated-histories
fi

git push -u origin main
echo "fatto -> $GITHUB_REPO_URL"
