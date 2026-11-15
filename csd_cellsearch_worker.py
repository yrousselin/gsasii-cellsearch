# -*- coding: utf-8 -*-
"""Worker run by g2_cellsearch.py inside the CCDC's own Python.

Reads a JSON job on stdin, runs a CSD reduced-cell search and prints one
line 'G2CSD_JSON:{...}' on stdout. Not meant to be run by hand.
"""
import json
import os
import sys


def _attr(obj, *names):
    for n in names:
        try:
            v = getattr(obj, n)
            if v is not None:
                return v
        except Exception:
            continue
    return None


def main():
    job = json.loads(sys.stdin.read())
    out = {'hits': []}
    try:
        import ccdc
        from ccdc.search import ReducedCellSearch
        out['csd_version'] = getattr(ccdc, '__version__', '')
        a, b, c, al, be, ga = job['cell']
        lengths, angles = (a, b, c), (al, be, ga)
        try:        # named tuples when available (same as crystal.cell_*)
            from ccdc.crystal import CellLengths, CellAngles
            lengths, angles = CellLengths(a, b, c), CellAngles(al, be, ga)
        except Exception:
            pass
        query = ReducedCellSearch.Query(lengths, angles, 'primitive')
        searcher = ReducedCellSearch(query)
        st = searcher.settings
        if job.get('percent_length_tolerance') is not None:
            try:
                st.percent_length_tolerance = float(
                    job['percent_length_tolerance'])
            except Exception as exc:
                out['warning'] = 'length tolerance not set: %s' % exc
        if job.get('absolute_angle_tolerance') is not None:
            try:
                st.absolute_angle_tolerance = float(
                    job['absolute_angle_tolerance'])
            except Exception as exc:
                out['warning'] = 'angle tolerance not set: %s' % exc
        hits = searcher.search()
        cifdir = job.get('cifdir')
        nmax = int(job.get('max_cifs') or 0)
        seen = set()
        for h in hits:
            ident = h.identifier
            if ident in seen:
                continue
            seen.add(ident)
            cr = h.crystal
            cl, ca = cr.cell_lengths, cr.cell_angles
            e = {'id': ident,
                 'cell': [float(x) for x in tuple(cl)[:3] + tuple(ca)[:3]]}
            sg = _attr(cr, 'spacegroup_symbol') or ''
            e['sg'] = str(sg)
            cent = str(_attr(cr, 'lattice_centring') or '')
            low = cent.lower()
            letter = 'P'
            if low.startswith('prim'):
                letter = 'P'
            elif low[:1] in 'abcfi' and low[:1]:
                letter = low[:1].upper()
            elif low.startswith('r'):
                letter = 'R'
            elif sg:
                letter = str(sg)[0]
            e['centring'] = letter
            try:
                ent = h.entry
                e['formula'] = str(_attr(ent, 'formula') or '')
                e['name'] = str(_attr(ent, 'chemical_name') or '')
                pub = _attr(ent, 'publication')
                if pub is not None:
                    e['ref'] = ' '.join(str(x) for x in (
                        _attr(pub, 'journal_name', 'journal') or '',
                        _attr(pub, 'year') or '') if x)
            except Exception:
                pass
            if cifdir and len(seen) <= nmax:
                path = os.path.join(cifdir, 'csd_%s.cif' % ident)
                try:
                    try:
                        txt = cr.to_string('cif')
                        with open(path, 'w') as fp:
                            fp.write(txt)
                    except Exception:
                        from ccdc.io import CrystalWriter
                        with CrystalWriter(path) as w:
                            w.write(cr)
                    e['cif'] = path
                except Exception:
                    pass
            out['hits'].append(e)
    except Exception:
        import traceback
        out['error'] = traceback.format_exc()
    sys.stdout.write('\nG2CSD_JSON:' + json.dumps(out) + '\n')


if __name__ == '__main__':
    main()
