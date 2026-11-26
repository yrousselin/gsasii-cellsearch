#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""g2_cellsearch.py -- search crystal-structure databases by unit cell.

Starting from a unit cell (an indexing result or a phase of a GSAS-II
project, or values typed on the command line), look for known structures
with the same lattice in:

  * COD   -- Crystallography Open Database (online, free)
  * CSD   -- Cambridge Structural Database (local CCDC installation,
             through the CSD Python API run in the CCDC's own Python)
  * local -- any folder of CIF files (e.g. CIFs exported by hand from
             ICSD Web)

Only numpy is needed. Typical use on the Mac:

  ~/g2main/bin/python g2_cellsearch.py project.gpx
  ~/g2main/bin/python g2_cellsearch.py project.gpx --list
  ~/g2main/bin/python g2_cellsearch.py --cell 8.86 10.45 13.42 101.3 99.9 92.7
  ~/g2main/bin/python g2_cellsearch.py project.gpx --local ~/ICSD_exports --icsd

How cells are compared
----------------------
Every cell (query and database entry) is converted to its primitive
lattice using the centring, then described by a setting-independent
"fingerprint": the cell volume and the lengths of the 7 shortest primitive
lattice vectors.  Two lattices match when all 7 lengths agree within --tol
(relative) and volumes within 3*--tol.  This is insensitive to the choice
of axes/setting and does not suffer from the instabilities of the Niggli
cell.  Supercells/subcells are NOT searched.
"""
from __future__ import annotations

import argparse
import csv
import io
import itertools
import json
import os
import pickle
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser
import hashlib

import numpy as np

__version__ = "1.5"

COD_URL = "https://www.crystallography.net/cod/result"
COD_CIF = "https://www.crystallography.net/cod/{}.cif"
ICSD_CELL_URL = "https://icsd-fiz-karlsruhe-de.inc.bib.cnrs.fr/search/cell.xhtml"
CONFIG = os.path.join(os.path.expanduser("~"), ".GSASII", "cellsearch.json")

# GSAS-II Bravais symbols, same order as GSASIIpwdGUI.bravaisSymb
BRAVAIS = ['Fm3m', 'Im3m', 'Pm3m', 'R3-H', 'P6/mmm', 'I4/mmm', 'P4/mmm',
           'Fmmm', 'Immm', 'Ammm', 'Bmmm', 'Cmmm', 'Pmmm', 'I2/m', 'A2/m',
           'C2/m', 'P2/m', 'P1', 'C1']

# ---------------------------------------------------------------------------
# lattice geometry
# ---------------------------------------------------------------------------
# rows = primitive basis vectors expressed on the conventional axes
_CENT = {
    'P': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    'I': [[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]],
    'F': [[0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
    'A': [[1, 0, 0], [0, .5, .5], [0, -.5, .5]],
    'B': [[.5, 0, .5], [0, 1, 0], [-.5, 0, .5]],
    'C': [[.5, .5, 0], [-.5, .5, 0], [0, 0, 1]],
    'R': [[2 / 3., 1 / 3., 1 / 3.], [-1 / 3., 1 / 3., 1 / 3.],
          [-1 / 3., -2 / 3., 1 / 3.]],          # obverse, hexagonal axes
}
_MULT = {'P': 1, 'I': 2, 'F': 4, 'A': 2, 'B': 2, 'C': 2, 'R': 3}
_GRIDS = {}


def _grid(n1, n2, n3):
    """Primitive integer triples (gcd = 1), one of each +/- pair."""
    key = (n1, n2, n3)
    if key not in _GRIDS:
        g = np.array([c for c in itertools.product(
            range(-n1, n1 + 1), range(-n2, n2 + 1), range(-n3, n3 + 1))
            if c > (0, 0, 0)], dtype=int)
        gcd = np.gcd.reduce(np.abs(g), axis=1)
        _GRIDS[key] = g[gcd == 1].astype(float)
    return _GRIDS[key]


def shortest_lengths(prim, n=7):
    """Lengths of the n shortest primitive lattice vectors (exact)."""
    lens = np.sort(np.sqrt(((_grid(2, 2, 2) @ prim) ** 2).sum(axis=1)))[:n]
    # every vector shorter than lens[-1] has |coef_i| <= lens[-1]*|a_i*|
    rec = np.linalg.inv(prim).T
    bound = np.floor(lens[-1] * np.sqrt((rec ** 2).sum(axis=1)) + 1e-9)
    bound = [int(min(max(b, 1), 40)) for b in bound]
    if max(bound) > 2:
        lens = np.sort(np.sqrt(((_grid(*bound) @ prim) ** 2).sum(axis=1)))[:n]
    return lens


NFP = 7


def cell_vectors(cell):
    """3x3 matrix whose rows are a, b, c (cartesian, Angstrom)."""
    a, b, c, al, be, ga = [float(x) for x in cell[:6]]
    al, be, ga = np.radians([al, be, ga])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 0:
        raise ValueError("impossible cell %s" % (cell,))
    return np.array([[a, 0, 0],
                     [b * np.cos(ga), b * np.sin(ga), 0],
                     [cx, cy, np.sqrt(cz2)]])


def vectors_cell(M):
    a, b, c = [np.linalg.norm(v) for v in M]
    ang = lambda u, v: np.degrees(np.arccos(np.clip(
        np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
    return [a, b, c, ang(M[1], M[2]), ang(M[0], M[2]), ang(M[0], M[1])]


def is_hex_axes(cell, tol=0.5):
    a, b, c, al, be, ga = cell[:6]
    return (abs(al - 90) < tol and abs(be - 90) < tol and abs(ga - 120) < tol
            and abs(a - b) < 0.005 * a)


def normalise_centring(cent, cell):
    """Return one of P I F A B C R (R only kept for hexagonal axes)."""
    c = (cent or 'P').strip().upper()[:1]
    if c == 'H':
        c = 'R'
    if c not in _CENT:
        c = 'P'
    if c == 'R' and not is_hex_axes(cell):
        c = 'P'                      # rhombohedral axes: already primitive
    return c


def reduce_basis(M):
    """Minkowski reduction of a 3x3 row basis: no vector can be shortened
    by adding integer combinations of the two others."""
    M = np.array(M, dtype=float)
    combos = [(s, t) for s in (-1, 0, 1) for t in (-1, 0, 1) if (s, t) != (0, 0)]
    for _ in range(1000):
        M = M[np.argsort([np.dot(v, v) for v in M])]
        changed = False
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                n = np.rint(np.dot(M[i], M[j]) / np.dot(M[j], M[j]))
                if n:
                    new = M[i] - n * M[j]
                    if np.dot(new, new) < np.dot(M[i], M[i]) * (1 - 1e-12):
                        M[i] = new
                        changed = True
        if not changed:
            for i in range(3):
                j, k = [x for x in range(3) if x != i]
                for s, t in combos:
                    new = M[i] + s * M[j] + t * M[k]
                    if np.dot(new, new) < np.dot(M[i], M[i]) * (1 - 1e-12):
                        M[i] = new
                        changed = True
        if not changed:
            break
    return M


class Lattice(object):
    """Primitive reduced lattice + fingerprint of a (cell, centring)."""

    def __init__(self, cell, centring='P'):
        self.cell = [float(x) for x in cell[:6]]
        self.centring = normalise_centring(centring, self.cell)
        conv = cell_vectors(self.cell)
        self.prim = reduce_basis(np.dot(np.array(_CENT[self.centring]), conv))
        self.volume = abs(np.linalg.det(self.prim))
        self.conv_volume = abs(np.linalg.det(conv))
        self.fp = shortest_lengths(self.prim, NFP)
        self.reduced_cell = vectors_cell(self.prim)

    def compare(self, other):
        """(max relative length deviation, relative volume deviation)."""
        dl = float(np.max(np.abs(other.fp - self.fp) / self.fp))
        dv = float(abs(other.volume - self.volume) / self.volume)
        return dl, dv

    def matches(self, other, tol):
        dl, dv = self.compare(other)
        return dl <= tol and dv <= 3 * tol


def fmt_cell(cell):
    return "%.4f %.4f %.4f  %.3f %.3f %.3f" % tuple(cell[:6])


# ---------------------------------------------------------------------------
# reading the cell from a GSAS-II project
# ---------------------------------------------------------------------------
class _Stub(object):
    def __init__(self, *a, **k):
        pass

    def __setstate__(self, state):
        self.__dict__['_state'] = state


class _Unpickler(pickle.Unpickler):
    """Read a .gpx without GSAS-II (unknown classes become stubs)."""

    def find_class(self, module, name):
        root = module.split('.')[0]
        if root in ('numpy', 'builtins', 'collections', 'copyreg',
                    '__builtin__', 'copy_reg', '_codecs'):
            if root == '__builtin__':
                module = 'builtins'
            if root == 'copy_reg':
                module = 'copyreg'
            try:
                return super().find_class(module, name)
            except Exception:
                pass
        return type(name, (_Stub,), {})


def read_gpx(path):
    """Return {top-level tree name: {child name: data, '': own data}}."""
    out = {}
    with open(path, 'rb') as fp:
        while True:
            try:
                item = _Unpickler(fp, encoding='latin-1').load()
            except EOFError:
                break
            if not item:
                continue
            d = {}
            for i, pair in enumerate(item):
                d['' if i == 0 else pair[0]] = pair[1]
            out[item[0][0]] = d
    return out


def cells_from_gpx(path):
    """List of candidate query cells: dicts(label, cell, centring, source)."""
    cands = []
    tree = read_gpx(path)
    for name, d in tree.items():
        if not name.startswith('PWDR'):
            continue
        ucl = d.get('Unit Cells List')
        if not ucl or len(ucl) < 3:
            continue
        controls, cells = ucl[0], ucl[2]
        used = [c for c in cells if len(c) > 11 and c[-2]]
        rows = used or list(cells[:1])
        for c in rows:
            try:
                sym = BRAVAIS[int(c[2])]
                cands.append(dict(
                    label="indexing %s, M20=%.1f [%s]%s" % (
                        sym, c[0], name[5:45], ' (use)' if c in used else ''),
                    cell=[float(x) for x in c[3:9]], centring=sym[0],
                    source='index'))
            except Exception:
                pass
    for pname, pdata in tree.get('Phases', {}).items():
        if not pname or not isinstance(pdata, dict):
            continue
        try:
            gen = pdata['General']
            sg = gen['SGData']['SpGrp']
            cands.append(dict(label="phase %s (%s)" % (pname, sg),
                              cell=[float(x) for x in gen['Cell'][1:7]],
                              centring=sg.strip()[0], source='phase'))
        except Exception:
            pass
    return cands


# ---------------------------------------------------------------------------
# hits
# ---------------------------------------------------------------------------
def make_hit(query, db, ident, cell, centring, sg='', formula='', name='',
             ref='', cif=''):
    try:
        lat = Lattice(cell, centring)
    except Exception:
        return None
    dl, dv = query.compare(lat)
    return dict(db=db, id=str(ident), cell=[float(x) for x in cell[:6]],
                centring=lat.centring, sg=sg or '', formula=formula or '',
                name=name or '', ref=ref or '', cif=cif or '',
                dlen=100 * dl, dvol=100 * dv, volume=lat.conv_volume,
                ok=None)


def _clean(s):
    if s is None:
        return ''
    return re.sub(r'\s+', ' ', str(s)).strip().strip("'\"- ")


def formula_elements(formula):
    return set(re.findall(r'[A-Z][a-z]?', formula or ''))


def element_filter(hits, elements, strict=False):
    """Keep hits whose formula contains all the elements (and, with
    strict, nothing else). Hits without formula are kept."""
    want = set(e.capitalize() for e in elements)
    if not want:
        return hits
    out = []
    for h in hits:
        have = formula_elements(h['formula'])
        if not have:
            out.append(h)
        elif want <= have and (not strict or have <= want):
            out.append(h)
    return out


# ---------------------------------------------------------------------------
# COD
# ---------------------------------------------------------------------------
def _cod_get(params, timeout):
    url = COD_URL + '?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={
        'User-Agent': 'g2_cellsearch/%s (GSAS-II user script)' % __version__})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        raw = r.read()
    txt = raw.decode('utf-8', 'replace').strip()
    if not txt or txt[0] not in '[{':
        return []
    data = json.loads(txt)
    if isinstance(data, dict):
        data = list(data.values())
    return [d for d in data if isinstance(d, dict)]


CACHE_DIR = os.path.join(os.path.expanduser("~"), ".GSASII", "cellsearch_cache")
CACHE_DAYS = 30
BIG_ANSWER = 4000     # entries; larger answers are re-queried in halves
_KEEP = ('file', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'sg', 'sgHall',
         'formula', 'mineral', 'commonname', 'chemname', 'journal', 'year')


def _cache_file(params):
    key = json.dumps(sorted(params.items()))
    return os.path.join(CACHE_DIR, hashlib.md5(key.encode()).hexdigest()
                        + '.json')


def _cache_get(params):
    path = _cache_file(params)
    try:
        if time.time() - os.path.getmtime(path) < CACHE_DAYS * 86400:
            with open(path) as fp:
                return json.load(fp)
    except Exception:
        pass
    return None


def _cache_put(params, rows):
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        slim = [{k: r.get(k) for k in _KEEP if r.get(k) is not None}
                for r in rows]
        with open(_cache_file(params), 'w') as fp:
            json.dump(slim, fp)
    except Exception:
        pass


def _f(x):
    try:
        return float(str(x).split('(')[0])
    except Exception:
        return None


COD_SQL = dict(host='sql.crystallography.net', user='cod_reader',
               password='', database='cod')
_SQLCOLS = ("file,a,b,c,alpha,beta,gamma,sg,formula,chemname,commonname,"
            "mineral,journal,year")


def _cod_sql_rows(windows, log=print):
    """Read the light columns straight from the public COD MySQL server.
    windows = [(m, vlo, vhi)]. Returns {m: rows} or raises."""
    import pymysql                      # optional dependency
    conn = pymysql.connect(connect_timeout=12, read_timeout=180,
                           charset='utf8mb4', **COD_SQL)
    out = {}
    try:
        cur = conn.cursor(pymysql.cursors.DictCursor)
        extra = " AND (status IS NULL OR status NOT IN ('errors','retracted'))"
        for m, lo, hi in windows:
            q = ("SELECT " + _SQLCOLS + " FROM data WHERE vol BETWEEN %s "
                 "AND %s")
            try:
                cur.execute(q + extra, (lo, hi))
            except pymysql.err.Error:
                extra = ''
                cur.execute(q, (lo, hi))
            out[m] = list(cur.fetchall())
    finally:
        conn.close()
    return out


def _cod_rest_rows(windows, base, elements, timeout, nslices, log=print):
    """REST fallback: many small volume slices, one at a time, cached."""
    def get(p):
        """one slice, with patience when the server refuses connections"""
        wait = 4
        for attempt in range(4):
            try:
                return _cod_get(p, timeout)
            except urllib.error.URLError as exc:
                if attempt == 3 or 'refused' not in str(exc).lower():
                    raise
                time.sleep(wait)
                wait *= 3
        return []

    def fetch(vlo, vhi, depth=0):
        p = dict(base)
        p['vmin'] = '%.3f' % vlo
        p['vmax'] = '%.3f' % vhi
        rows = _cache_get(p)
        if rows is not None:
            return rows, []
        try:
            rows = get(p)
            time.sleep(0.3)             # be gentle with the server
            # very long answers may be capped by the server: an answer
            # this large is never trusted, the slice is halved instead
            if len(rows) < BIG_ANSWER or (vhi - vlo) < 0.5:
                _cache_put(p, rows)
                return rows, []
            last = "answer too long (%d entries)" % len(rows)
        except urllib.error.URLError as exc:
            if 'refused' in str(exc).lower():     # halving would not help
                return [], ["%.1f-%.1f A3 : %s" % (vlo, vhi, exc)]
            last = exc
        except Exception as exc:       # timeout, truncated JSON...
            last = exc
        if depth >= 4 or (vhi - vlo) < 0.5:
            return [], ["%.1f-%.1f A3 : %s" % (vlo, vhi, last)]
        mid = 0.5 * (vlo + vhi)
        r1, e1 = fetch(vlo, mid, depth + 1)
        r2, e2 = fetch(mid + 0.001, vhi, depth + 1)
        return r1 + r2, e1 + e2

    nsl = 1 if elements else nslices
    jobs = []
    for m, lo, hi in windows:
        w = (hi - lo) / nsl
        for k in range(nsl):
            jobs.append((m, lo + k * w + (0.001 if k else 0), lo + (k + 1) * w))
    log("  COD (web): %d requests%s ..." % (len(jobs), '' if elements else
        " (no element filter: this may take several minutes)"))
    out, errors, n = {m: [] for m, _lo, _hi in windows}, [], 0
    for i, (m, lo, hi) in enumerate(jobs):
        rows, errs = fetch(lo, hi)
        out[m] += rows
        n += len(rows)
        errors += ["V x%d, %s" % (m, e) for e in errs]
        if (i + 1) % 4 == 0 or i + 1 == len(jobs):
            log("  COD (web): %d/%d requests, %d entries read" % (
                i + 1, len(jobs), n))
    if errors:
        errors.append("run the search again: only the failed slices "
                      "will be requested again (cache)")
    return out, errors


def search_cod(query, tol, elements=(), strict=False, timeout=90,
               nslices=8, mode='auto', log=print):
    """Find COD entries with the same lattice. Entries are fetched by
    volume windows (x1,2,3,4 the primitive volume) and filtered locally
    on the lattice fingerprint.

    mode 'sql'  : public COD MySQL server (fast, needs pymysql + port 3306)
         'rest' : web API (slow without an element filter)
         'auto' : sql when possible, otherwise rest
    """
    base = {'format': 'json'}
    for i, el in enumerate(list(elements)[:8]):
        base['el%d' % (i + 1)] = el
    if strict and elements:
        base['strictmax'] = len(elements)
    windows = [(m, m * query.volume * (1 - 3 * tol),
                m * query.volume * (1 + 3 * tol)) for m in (1, 2, 3, 4)]
    bym, errors = None, []
    if mode in ('auto', 'sql'):
        try:
            t0 = time.time()
            bym = _cod_sql_rows(windows, log=log)
            log("  COD (SQL): read in %.1f s" % (time.time() - t0))
        except ImportError:
            msg = ("pymysql module missing (install it with: "
                   "%s -m pip install pymysql)" % sys.executable)
            log("  COD (SQL): " + msg)
            if mode == 'sql':
                return [], [msg]
        except Exception as exc:
            log("  COD (SQL): unavailable (%s)" % str(exc)[:150])
            if mode == 'sql':
                return [], [str(exc)]
    if bym is None:
        bym, errors = _cod_rest_rows(windows, base, elements, timeout,
                                     nslices, log=log)
    hits, nrows = {}, 0
    for m, rows in bym.items():
        nrows += len(rows)
        for r in rows:
            cell = [_f(r.get(k)) for k in
                    ('a', 'b', 'c', 'alpha', 'beta', 'gamma')]
            if None in cell or min(cell) <= 0:
                continue
            sg = _clean(r.get('sg') or r.get('sgHall'))
            cent = sg[:1] if sg else 'P'
            if _MULT.get(normalise_centring(cent, cell), 1) != m:
                continue
            cid = r.get('file') or r.get('id')
            h = make_hit(
                query, 'COD', cid, cell, cent, sg=sg,
                formula=_clean(r.get('formula')),
                name=_clean(r.get('mineral') or r.get('commonname')
                            or r.get('chemname')),
                ref=_clean('%s %s' % (r.get('journal') or '',
                                      r.get('year') or '')),
                cif=COD_CIF.format(cid))
            if h and h['dlen'] <= 100 * tol and h['dvol'] <= 300 * tol:
                hits[h['id']] = h
    log("  COD: entries per volume window  " + "  ".join(
        "x%d: %d" % (m, len(bym[m])) for m in sorted(bym)))
    log("  COD: %d entries examined, %d match(es)%s" % (
        nrows, len(hits), ''.join('\n    ERROR ' + e for e in errors)))
    return list(hits.values()), errors


def download_cod_cifs(hits, outdir, nmax=10, timeout=60, log=print):
    n = 0
    for h in hits:
        if h['db'] != 'COD' or n >= nmax:
            continue
        dest = os.path.join(outdir, 'cod_%s.cif' % h['id'])
        try:
            if not os.path.exists(dest):
                req = urllib.request.Request(COD_CIF.format(h['id']), headers={
                    'User-Agent': 'g2_cellsearch/%s' % __version__})
                with urllib.request.urlopen(req, timeout=timeout) as r, \
                        open(dest, 'wb') as fp:
                    fp.write(r.read())
            h['cif'] = dest
            n += 1
        except Exception as exc:
            log("  COD %s: CIF not downloaded (%s)" % (h['id'], exc))
    return n


# ---------------------------------------------------------------------------
# local CIF folder (e.g. ICSD exports)
# ---------------------------------------------------------------------------
_TAGS = {
    'a': ('_cell_length_a',), 'b': ('_cell_length_b',),
    'c': ('_cell_length_c',), 'alpha': ('_cell_angle_alpha',),
    'beta': ('_cell_angle_beta',), 'gamma': ('_cell_angle_gamma',),
    'sg': ('_space_group_name_h-m_alt', '_symmetry_space_group_name_h-m',
           '_space_group_name_h-m'),
    'formula': ('_chemical_formula_sum', '_chemical_formula_structural'),
    'name': ('_chemical_name_mineral', '_chemical_name_common',
             '_chemical_name_systematic', '_chemical_name_structure_type'),
    'code': ('_database_code_icsd', '_database_code_depnum_ccdc_archive',
             '_cod_database_code', '_database_code_amcsd'),
}
_TAGRE = re.compile(r"^\s*(_[^\s]+)[ \t]+(?:'([^']*)'|\"([^\"]*)\"|(\S+))",
                    re.M)


def parse_cif_blocks(text):
    """Very light CIF reader: yields (blockname, dict) with cell/sg/formula."""
    parts = re.split(r'(?im)^\s*data_', text)
    for part in parts[1:]:
        bname = part.split(None, 1)[0] if part.strip() else ''
        found = {}
        for m in _TAGRE.finditer(part):
            found.setdefault(m.group(1).lower(),
                             m.group(2) or m.group(3) or m.group(4))
        d = {}
        for key, tags in _TAGS.items():
            for t in tags:
                if t in found and found[t] not in ('?', '.'):
                    d[key] = found[t]
                    break
        yield bname, d


def search_local(query, tol, folder, log=print):
    hits, nfiles, nblocks = [], 0, 0
    folder = os.path.expanduser(folder)
    for root, _dirs, files in os.walk(folder):
        for fn in sorted(files):
            if not fn.lower().endswith('.cif'):
                continue
            path = os.path.join(root, fn)
            nfiles += 1
            try:
                with io.open(path, 'r', encoding='utf-8',
                             errors='replace') as fp:
                    text = fp.read()
            except Exception:
                continue
            for bname, d in parse_cif_blocks(text):
                cell = [_f(d.get(k)) for k in
                        ('a', 'b', 'c', 'alpha', 'beta', 'gamma')]
                if None in cell:
                    continue
                nblocks += 1
                sg = _clean(d.get('sg'))
                ident = d.get('code') or bname or fn
                h = make_hit(query, 'local', ident, cell, sg[:1] or 'P',
                             sg=sg, formula=_clean(d.get('formula')),
                             name=_clean(d.get('name')),
                             ref=os.path.relpath(path, folder), cif=path)
                if h and h['dlen'] <= 100 * tol and h['dvol'] <= 300 * tol:
                    hits.append(h)
    log("  local folder: %d CIF file(s), %d structure(s), "
        "%d match(es)" % (nfiles, nblocks, len(hits)))
    return hits


# ---------------------------------------------------------------------------
# CSD (through the CCDC's own Python, as a sub-process)
# ---------------------------------------------------------------------------
def _load_config():
    try:
        with open(CONFIG) as fp:
            return json.load(fp)
    except Exception:
        return {}


def _save_config(cfg):
    try:
        os.makedirs(os.path.dirname(CONFIG), exist_ok=True)
        with open(CONFIG, 'w') as fp:
            json.dump(cfg, fp, indent=1)
    except Exception:
        pass


def _python_has_ccdc(cmd, timeout=120):
    try:
        r = subprocess.run(list(cmd) + ['-c', 'import ccdc;print(ccdc.__version__)'],
                           capture_output=True, text=True, timeout=timeout)
        if r.returncode == 0:
            return True, (r.stdout.strip().splitlines() or ['?'])[-1]
        return False, (r.stdout + r.stderr).strip()[-300:]
    except Exception as exc:
        return False, str(exc)


def find_csd_python(csd_home=None, explicit=None, log=print):
    """Locate a Python that can 'import ccdc'. Returns a command list."""
    cfg = _load_config()
    tried = []
    cands = []
    if explicit:
        cands.append(os.path.expanduser(explicit))
    if os.environ.get('G2_CSD_PYTHON'):
        cands.append(os.environ['G2_CSD_PYTHON'])
    if cfg.get('csd_python'):
        cands.append(cfg['csd_python'])
    homes = [csd_home] if csd_home else []
    homes += [cfg.get('csd_home'), '~/CCDC', '/Applications/CCDC',
              '/opt/CCDC', 'C:/Program Files/CCDC']
    launch, pythons = [], []
    for home in homes:
        if not home:
            continue
        home = os.path.expanduser(home)
        if not os.path.isdir(home):
            continue
        depth0 = home.rstrip(os.sep).count(os.sep)
        for root, dirs, files in os.walk(home):
            if root.count(os.sep) - depth0 >= 7:
                dirs[:] = []
                continue
            dirs[:] = [d for d in dirs if d not in (
                'pkgs', 'site-packages', 'node_modules', 'Frameworks',
                'Resources', 'share', 'include', 'doc', 'docs')]
            for fn in files:
                low = fn.lower()
                full = os.path.join(root, fn)
                if low.startswith('run_csd_python_api'):
                    launch.append(full)
                elif (re.match(r'python(3(\.\d+)?)?(\.exe)?$', low)
                      and os.path.basename(root) in ('bin', 'Scripts',
                                                     'miniconda', 'python')
                      and ('api' in root.lower() or 'conda' in root.lower()
                           or 'ccdc' in root.lower())):
                    pythons.append(full)
    cands += sorted(launch) + sorted(pythons, key=len)
    seen = set()
    for c in cands:
        if c in seen or not os.path.exists(c):
            continue
        seen.add(c)
        ok, msg = _python_has_ccdc([c])
        tried.append((c, msg))
        if ok:
            cfg['csd_python'] = c
            _save_config(cfg)
            log("  CSD: CCDC Python = %s (ccdc %s)" % (
                c, msg.splitlines()[-1] if msg else '?'))
            return [c], tried
    return None, tried


def search_csd(query, tol, outdir, csd_home=None, csd_python=None,
               angle_tol=None, max_cifs=25, timeout=900, log=print):
    worker = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'csd_cellsearch_worker.py')
    if not os.path.exists(worker):
        return [], ["csd_cellsearch_worker.py not found next to the script"]
    cmd, tried = find_csd_python(csd_home, csd_python, log=log)
    if cmd is None:
        msg = ["CCDC Python not found (use --csd-python PATH)."]
        msg += ["   tried: %s -> %s" % t for t in tried[-4:]]
        log("  CSD : " + "\n".join(msg))
        return [], msg
    cifdir = os.path.join(outdir, 'csd')
    os.makedirs(cifdir, exist_ok=True)
    job = dict(cell=[float(x) for x in query.reduced_cell],
               percent_length_tolerance=100 * tol,
               absolute_angle_tolerance=angle_tol,
               cifdir=cifdir, max_cifs=max_cifs)
    try:
        r = subprocess.run(cmd + [worker], input=json.dumps(job),
                           capture_output=True, text=True, timeout=timeout)
    except Exception as exc:
        log("  CSD: could not start (%s)" % exc)
        return [], [str(exc)]
    out = None
    for line in reversed(r.stdout.splitlines()):
        if line.startswith('G2CSD_JSON:'):
            out = json.loads(line[len('G2CSD_JSON:'):])
            break
    if out is None or out.get('error'):
        err = (out or {}).get('error') or (r.stderr or r.stdout)[-1500:]
        log("  CSD: error\n%s" % err)
        return [], [err]
    hits = []
    for e in out.get('hits', []):
        h = make_hit(query, 'CSD', e['id'], e['cell'], e.get('centring', 'P'),
                     sg=e.get('sg', ''), formula=e.get('formula', ''),
                     name=e.get('name', ''), ref=e.get('ref', ''),
                     cif=e.get('cif', ''))
        if h is None:
            continue
        # the CSD has its own tolerance definition: keep everything it
        # returned, but flag entries outside our own criterion
        h['ok'] = h['dlen'] <= 100 * tol and h['dvol'] <= 300 * tol
        hits.append(h)
    log("  CSD %s: %d match(es) returned by the "
        "reduced-cell search" % (out.get('csd_version', ''), len(hits)))
    return hits, []


# ---------------------------------------------------------------------------
# ICSD helper (manual search in the browser -- no automated access)
# ---------------------------------------------------------------------------
def open_icsd(query_cell, url=ICSD_CELL_URL, log=print):
    txt = "%.4f %.4f %.4f %.3f %.3f %.3f" % tuple(query_cell[:6])
    copied = False
    for cmd in (['pbcopy'], ['clip'], ['xclip', '-selection', 'clipboard']):
        try:
            subprocess.run(cmd, input=txt, text=True, check=True, timeout=5)
            copied = True
            break
        except Exception:
            continue
    try:
        webbrowser.open(url)
    except Exception:
        pass
    log("  ICSD: unit-cell search page opened in the browser.")
    log("         cell %s: %s" % ("copied to the clipboard"
                                     if copied else "to be typed", txt))
    log("         Export the CIF files to a folder, then run again with "
        "--local FOLDER")
    return txt


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run_search(cell, centring='P', tol=0.015, use_cod=True, use_csd=True,
               local=None, elements=(), strict=False, outdir=None,
               csd_home=None, csd_python=None, csd_angle_tol=None,
               max_cifs=10, cod_mode='auto', log=print):
    """Run all searches. Returns (query Lattice, sorted hits, errors)."""
    query = Lattice(cell, centring)
    log("g2_cellsearch version %s" % __version__)
    log("Cell searched: %s  (%s)" % (fmt_cell(query.cell),
                                           query.centring))
    log("Reduced primitive cell: %s   V = %.2f A3" % (
        fmt_cell(query.reduced_cell), query.volume))
    log("Tolerance: %.2f %% on lengths, %.2f %% on volume" % (
        100 * tol, 300 * tol))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    hits, errors = [], []
    if use_cod:
        h, e = search_cod(query, tol, elements, strict, mode=cod_mode,
                          log=log)
        hits += h
        errors += ['COD: ' + x for x in e]
    if use_csd:
        if not outdir:
            raise ValueError("outdir needed for the CSD search")
        h, e = search_csd(query, tol, outdir, csd_home, csd_python,
                          csd_angle_tol, max_cifs=500, log=log)
        hits += h
        errors += ['CSD: ' + str(x) for x in e]
    for folder in (local or []):
        if os.path.isdir(os.path.expanduser(folder)):
            hits += search_local(query, tol, folder, log=log)
        else:
            errors.append("local folder not found: %s" % folder)
            log("  " + errors[-1])
    if elements:
        n0 = len(hits)
        hits = element_filter(hits, elements, strict)
        log("  element filter %s%s: %d -> %d" % (
            ' '.join(elements), ' (strict)' if strict else '', n0, len(hits)))
    hits.sort(key=lambda h: (h['dlen'], h['dvol']))
    if outdir and use_cod:
        download_cod_cifs(hits, outdir, nmax=max_cifs, log=log)
    # summary of the search, written by write_csv next to the CSV (cellsearch_summary.json)
    bases = [b for b, used in (('COD', use_cod), ('CSD', use_csd)) if used]
    query.search_info = {
        'format': 1, 'kind': 'cell', 'date': time.strftime('%Y-%m-%d %H:%M'),
        'cell': [float(x) for x in query.cell], 'centring': query.centring, 'tol_percent': 100 * tol,
        'elements': list(elements), 'strict': bool(strict), 'local': list(local or []),
        'databases': {b: {'hits': sum(1 for h in hits if h['db'] == b),
                          'error': '; '.join(e[5:] for e in errors if e.startswith(b + ': '))} for b in bases},
        'n_hits': len(hits)}
    return query, hits, errors


def format_table(hits, nmax=40):
    if not hits:
        return "No match."
    lines = ["%-3s %-6s %-12s %6s %6s  %-11s %-44s %s" % (
        '#', 'db', 'id', 'dL %', 'dV %', 'group', 'cell', 'formula / name')]
    for i, h in enumerate(hits[:nmax]):
        desc = h['formula']
        if h['name']:
            desc = (desc + ' | ' if desc else '') + h['name']
        flag = '' if h['ok'] is not False else ' (outside tol.)'
        lines.append("%-3d %-6s %-12s %6.2f %6.2f  %-11s %-44s %s%s" % (
            i + 1, h['db'], h['id'][:12], h['dlen'], h['dvol'], h['sg'][:11],
            fmt_cell(h['cell']), desc[:70], flag))
    if len(hits) > nmax:
        lines.append("... %d more (see the CSV file)" % (len(hits) - nmax))
    return "\n".join(lines)


def write_csv(hits, path, query):
    with io.open(path, 'w', newline='', encoding='utf-8') as fp:
        w = csv.writer(fp)
        w.writerow(['# query cell'] + ['%.5f' % x for x in query.cell]
                   + [query.centring])
        w.writerow(['rank', 'database', 'id', 'dlen_percent', 'dvol_percent',
                    'space_group', 'a', 'b', 'c', 'alpha', 'beta', 'gamma',
                    'volume', 'formula', 'name', 'reference', 'cif'])
        for i, h in enumerate(hits):
            w.writerow([i + 1, h['db'], h['id'], '%.3f' % h['dlen'],
                        '%.3f' % h['dvol'], h['sg']]
                       + ['%.5f' % x for x in h['cell']]
                       + ['%.3f' % h['volume'], h['formula'], h['name'],
                          h['ref'], h['cif']])
    info = getattr(query, 'search_info', None)
    if info:
        try:
            with io.open(os.path.join(os.path.dirname(path), 'cellsearch_summary.json'), 'w', encoding='utf-8') as fp:
                fp.write(json.dumps(info, ensure_ascii=False, indent=1))
        except (OSError, TypeError, ValueError):
            pass


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Search for known structures by unit cell "
                    "(COD, CSD, local CIF folder).")
    ap.add_argument('gpx', nargs='?', help="GSAS-II project (.gpx)")
    ap.add_argument('--cell', nargs=6, type=float,
                    metavar=('a', 'b', 'c', 'al', 'be', 'ga'))
    ap.add_argument('--centring', default='P',
                    help="P I F A B C R (with --cell; default P)")
    ap.add_argument('--list', action='store_true',
                    help="list the cells found in the project")
    ap.add_argument('--pick', type=int, default=None,
                    help="number of the cell to use (see --list)")
    ap.add_argument('--tol', type=float, default=1.5,
                    help="tolerance on lengths in %% (default 1.5)")
    ap.add_argument('--elements', default='',
                    help="required elements (all databases), e.g. 'Zn N C H'")
    ap.add_argument('--strict', action='store_true',
                    help="only these elements, no others")
    ap.add_argument('--no-cod', action='store_true')
    ap.add_argument('--cod-mode', default='auto',
                    choices=('auto', 'sql', 'rest'),
                    help="COD access: sql (fast, pymysql) or rest (web)")
    ap.add_argument('--no-csd', action='store_true')
    ap.add_argument('--local', action='append', default=[],
                    help="folder of CIF files (repeatable), e.g. ICSD exports")
    ap.add_argument('--icsd', action='store_true',
                    help="open the ICSD page and copy the cell")
    ap.add_argument('--icsd-url', default=ICSD_CELL_URL)
    ap.add_argument('--csd-home', default=None, help="default ~/CCDC")
    ap.add_argument('--csd-python', default=None)
    ap.add_argument('--csd-angle-tol', type=float, default=None)
    ap.add_argument('--max-cifs', type=int, default=10,
                    help="number of COD CIF files downloaded (default 10)")
    ap.add_argument('--out', default=None, help="output folder")
    a = ap.parse_args(argv)

    if a.cell:
        cell, cent, label = a.cell, a.centring, 'command line'
        stem = os.path.abspath('cellsearch')
    elif a.gpx:
        cands = cells_from_gpx(a.gpx)
        if a.list or not cands:
            if not cands:
                print("No cell (indexing or phase) in this project.")
                return 1
            for i, c in enumerate(cands):
                print("%2d  %s (%s)  %s" % (i + 1, fmt_cell(c['cell']),
                                           c['centring'], c['label']))
            return 0
        k = (a.pick - 1) if a.pick else 0
        if not 0 <= k < len(cands):
            print("--pick out of range (1..%d)" % len(cands))
            return 1
        cell, cent, label = cands[k]['cell'], cands[k]['centring'], \
            cands[k]['label']
        stem = os.path.splitext(os.path.abspath(a.gpx))[0]
    else:
        ap.error("give a .gpx or --cell")
    outdir = a.out or (stem + '_cellsearch')
    print("Cell source: %s" % label)
    query, hits, errors = run_search(
        cell, cent, tol=a.tol / 100., use_cod=not a.no_cod,
        use_csd=not a.no_csd, local=a.local,
        elements=a.elements.replace(',', ' ').split(), strict=a.strict,
        outdir=outdir, csd_home=a.csd_home, csd_python=a.csd_python,
        csd_angle_tol=a.csd_angle_tol, max_cifs=a.max_cifs,
        cod_mode=a.cod_mode)
    print()
    print(format_table(hits))
    csvpath = os.path.join(outdir, 'cellsearch_results.csv')
    write_csv(hits, csvpath, query)
    print("\nResults: %s\nCIF    : %s" % (csvpath, outdir))
    if a.icsd:
        open_icsd(query.cell, a.icsd_url)
    return 0


if __name__ == '__main__':
    sys.exit(main())
