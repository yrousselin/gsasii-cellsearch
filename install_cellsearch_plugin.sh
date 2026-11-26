#!/bin/bash
# Installs the unit-cell search into GSAS-II (Export menu) -- macOS/Linux
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.GSASII/exports"
mkdir -p "$DEST"
cp "$HERE/G2export_CellSearch.py" "$HERE/g2_cellsearch.py" \
   "$HERE/csd_cellsearch_worker.py" "$DEST/"
echo "Files copied to $DEST"

PY="$HOME/g2main/bin/python"
[ -x "$PY" ] || PY="$(command -v python3)"
echo
echo "Installing pymysql (fast COD access)..."
"$PY" -m pip install --quiet pymysql || echo "pymysql not installed: the COD will be queried through the web (slower)"
echo
echo "Looking for the CCDC Python (may take a minute)..."
"$PY" - <<EOF
import sys
sys.path.insert(0, "$DEST")
import g2_cellsearch as cs
cmd, tried = cs.find_csd_python()
if cmd:
    print("OK: the CSD search will use", cmd[0])
else:
    print("CCDC Python not found automatically.")
    for t in tried: print("  tried:", t[0], "->", t[1][-120:])
    print("The COD and local-folder searches still work.")
    print("For the CSD, give its path with --csd-python, or report the output of:")
    print('  find ~/CCDC -maxdepth 7 \\\\( -name "python3*" -o -name "run_csd_python_api*" \\\\) | head -20')
EOF
echo
echo "Restart GSAS-II: Export -> Entire project as -> (user) Search databases by unit cell"
