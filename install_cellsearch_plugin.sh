#!/bin/bash
# Installe la recherche par maille dans GSAS-II (menu Export) -- macOS/Linux
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.GSASII/exports"
mkdir -p "$DEST"
cp "$HERE/G2export_CellSearch.py" "$HERE/g2_cellsearch.py" \
   "$HERE/csd_cellsearch_worker.py" "$DEST/"
echo "Fichiers copies dans $DEST"

PY="$HOME/g2main/bin/python"
[ -x "$PY" ] || PY="$(command -v python3)"
echo
echo "Installation de pymysql (acces rapide a la COD)..."
"$PY" -m pip install --quiet pymysql || echo "pymysql non installe : la COD passera par le web (plus lent)"
echo
echo "Recherche du Python du CCDC (peut prendre une minute)..."
"$PY" - <<EOF
import sys
sys.path.insert(0, "$DEST")
import g2_cellsearch as cs
cmd, tried = cs.find_csd_python()
if cmd:
    print("OK : la recherche CSD utilisera", cmd[0])
else:
    print("Python du CCDC non trouve automatiquement.")
    for t in tried: print("  essaye :", t[0], "->", t[1][-120:])
    print("La recherche COD et dossier local fonctionne quand meme.")
    print("Pour la CSD, envoyez le resultat de :")
    print('  find ~/CCDC -maxdepth 7 \\\\( -name "python3*" -o -name "run_csd_python_api*" \\\\) | head -20')
EOF
echo
echo "Redemarrez GSAS-II : Export -> Entire project as -> (user) Search databases by unit cell"
