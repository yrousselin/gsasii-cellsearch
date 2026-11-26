# -*- coding: utf-8 -*-
'''User "exporter" that adds to GSAS-II the menu entry

   Export -> Entire project as -> (user) Search databases by unit cell

It takes a unit cell of the current project (indexing result flagged "use"
in Unit Cells List, or a phase), looks for known structures with the same
lattice in COD, CSD and local CIF folders (e.g. ICSD exports), and lists
the matches. The CIFs are saved next to the project so that they can be
imported with Import -> Phase -> from CIF file.

All the work is done by g2_cellsearch.py (same folder).
'''
from __future__ import division, print_function
import os
import subprocess
import sys

from .. import GSASIIfiles as G2fil

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)


class ExportCellSearch(G2fil.ExportBaseclass):
    '''Search structure databases by unit cell

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self, G2frame):
        super(self.__class__, self).__init__(
            G2frame=G2frame,
            formatName='Search databases by unit cell',
            extension='.csv',
            longFormatName='Search COD / CSD / local CIFs by unit cell')
        self.exporttype = ['project']
        self.multiple = False

    # ------------------------------------------------------------------
    def _options(self, cs, cands, cfg):
        '''Dialog asking for the cell and the search options.'''
        import wx
        dlg = wx.Dialog(self.G2frame, title='Search databases by unit cell',
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        sz = wx.BoxSizer(wx.VERTICAL)
        sz.Add(wx.StaticText(dlg, label='Unit cell to search for:'), 0,
               wx.ALL, 5)
        labels = ['%s (%s)   %s' % (cs.fmt_cell(c['cell']), c['centring'],
                                    c['label']) for c in cands]
        cell = wx.Choice(dlg, choices=labels)
        cell.SetSelection(0)
        sz.Add(cell, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)

        grid = wx.FlexGridSizer(cols=2, hgap=8, vgap=6)
        grid.AddGrowableCol(1, 1)
        grid.Add(wx.StaticText(dlg, label='Length tolerance (%)'), 0,
                 wx.ALIGN_CENTER_VERTICAL)
        tol = wx.TextCtrl(dlg, value=str(cfg.get('tol', 1.5)))
        grid.Add(tol, 1, wx.EXPAND)
        grid.Add(wx.StaticText(dlg, label='Required elements (optional)'),
                 0, wx.ALIGN_CENTER_VERTICAL)
        elem = wx.TextCtrl(dlg, value='')
        elem.SetToolTip('e.g.  Zn N C H   -- filters all databases and speeds up the COD query')
        grid.Add(elem, 1, wx.EXPAND)
        grid.Add(wx.StaticText(dlg, label='Local CIF folder (ICSD exports)'),
                 0, wx.ALIGN_CENTER_VERTICAL)
        row = wx.BoxSizer(wx.HORIZONTAL)
        local = wx.TextCtrl(dlg, value=cfg.get('local', ''))
        row.Add(local, 1, wx.EXPAND)
        btn = wx.Button(dlg, label='Browse...')

        def onBrowse(evt):
            d = wx.DirDialog(dlg, 'Folder of CIF files',
                             local.GetValue() or os.path.expanduser('~'))
            if d.ShowModal() == wx.ID_OK:
                local.SetValue(d.GetPath())
            d.Destroy()
        btn.Bind(wx.EVT_BUTTON, onBrowse)
        row.Add(btn, 0, wx.LEFT, 4)
        grid.Add(row, 1, wx.EXPAND)
        sz.Add(grid, 0, wx.EXPAND | wx.ALL, 8)

        cod = wx.CheckBox(dlg, label='COD (online)')
        cod.SetValue(cfg.get('cod', True))
        csd = wx.CheckBox(dlg, label='CSD (local CCDC installation)')
        csd.SetValue(cfg.get('csd', True))
        icsd = wx.CheckBox(dlg, label='Also open the ICSD cell-search page '
                           'and copy the cell to the clipboard')
        icsd.SetValue(False)
        for w in (cod, csd, icsd):
            sz.Add(w, 0, wx.LEFT | wx.BOTTOM, 8)
        sz.Add(dlg.CreateButtonSizer(wx.OK | wx.CANCEL), 0,
               wx.EXPAND | wx.ALL, 8)
        dlg.SetSizerAndFit(sz)
        dlg.SetSize((max(620, dlg.GetSize()[0]), -1))
        dlg.CenterOnParent()
        res = None
        if dlg.ShowModal() == wx.ID_OK:
            try:
                t = float(tol.GetValue().replace(',', '.'))
            except ValueError:
                t = 1.5
            res = dict(cand=cands[cell.GetSelection()], tol=t,
                       elements=elem.GetValue().replace(',', ' ').split(),
                       local=local.GetValue().strip(), cod=cod.GetValue(),
                       csd=csd.GetValue(), icsd=icsd.GetValue())
        dlg.Destroy()
        return res

    def _results(self, cs, query, hits, errors, outdir, csvpath, logtxt):
        '''Dialog listing the matches.'''
        import wx
        dlg = wx.Dialog(self.G2frame, title='Unit-cell search results',
                        size=(1050, 560),
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        sz = wx.BoxSizer(wx.VERTICAL)
        head = 'Query: %s (%s)   primitive V = %.2f A^3   --   %d match(es)' % (
            cs.fmt_cell(query.cell), query.centring, query.volume, len(hits))
        sz.Add(wx.StaticText(dlg, label=head), 0, wx.ALL, 6)
        lc = wx.ListCtrl(dlg, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        cols = [('#', 35), ('database', 55), ('id', 95), ('dL %', 55),
                ('dV %', 55), ('space group', 95), ('cell', 300),
                ('formula', 170), ('name', 220), ('CIF', 45)]
        for i, (lbl, w) in enumerate(cols):
            lc.InsertColumn(i, lbl, width=w)
        for i, h in enumerate(hits):
            vals = [str(i + 1), h['db'], h['id'], '%.2f' % h['dlen'],
                    '%.2f' % h['dvol'], h['sg'], cs.fmt_cell(h['cell']),
                    h['formula'], h['name'],
                    'yes' if os.path.isfile(h['cif']) else '']
            idx = lc.InsertItem(i, vals[0])
            for j, v in enumerate(vals[1:]):
                lc.SetItem(idx, j + 1, v)
        sz.Add(lc, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        msg = ('CIFs are in %s\nImport one with  Import -> Phase -> '
               'from CIF file.' % outdir)
        if errors:
            msg += '\nProblems: ' + ' ; '.join(str(e)[:160] for e in errors)
        sz.Add(wx.StaticText(dlg, label=msg), 0, wx.ALL, 6)
        row = wx.BoxSizer(wx.HORIZONTAL)

        def reveal(path):
            try:
                isfile = os.path.isfile(path)
                folder = path if os.path.isdir(path) else os.path.dirname(path)
                if sys.platform == 'darwin':
                    subprocess.Popen(['open', '-R', path] if isfile
                                     else ['open', path])
                elif sys.platform.startswith('win'):
                    os.startfile(folder)
                else:
                    subprocess.Popen(['xdg-open', folder])
            except Exception as exc:
                print('could not open %s: %s' % (path, exc))

        def onShow(evt):
            i = lc.GetFirstSelected()
            if i < 0 or not os.path.isfile(hits[i]['cif']):
                reveal(outdir)
            else:
                reveal(hits[i]['cif'])

        def onLog(evt):
            d = wx.MessageDialog(dlg, logtxt, 'Search log', wx.OK)
            d.ShowModal()
            d.Destroy()

        def onICSD(evt):
            cs.open_icsd(query.cell)

        for lbl, fn in (('Show selected CIF / folder', onShow),
                        ('Open ICSD cell search', onICSD),
                        ('Log', onLog)):
            b = wx.Button(dlg, label=lbl)
            b.Bind(wx.EVT_BUTTON, fn)
            row.Add(b, 0, wx.RIGHT, 6)
        row.AddStretchSpacer()
        b = wx.Button(dlg, wx.ID_OK, 'Close')
        row.Add(b, 0)
        sz.Add(row, 0, wx.EXPAND | wx.ALL, 6)
        dlg.SetSizer(sz)
        dlg.CenterOnParent()
        dlg.ShowModal()
        dlg.Destroy()

    # ------------------------------------------------------------------
    def Exporter(self, event=None):
        '''Run the search (called from the Export menu)'''
        import wx
        import importlib
        import g2_cellsearch as cs
        importlib.reload(cs)        # pick up edits without restarting
        self.InitExport(event)
        G2frame = self.G2frame
        if not G2frame.GSASprojectfile:
            wx.MessageBox('Save the project first.', 'Cell search')
            return
        G2frame.OnFileSave(event)   # the search reads the .gpx file
        gpx = G2frame.GSASprojectfile
        try:
            cands = cs.cells_from_gpx(gpx)
        except Exception as exc:
            wx.MessageBox('Could not read %s\n%s' % (gpx, exc), 'Cell search')
            return
        if not cands:
            wx.MessageBox('No unit cell found: index the pattern (and flag '
                          'a cell with "use") or define a phase first.',
                          'Cell search')
            return
        cfg = cs._load_config()
        opt = self._options(cs, cands, cfg.get('gui', {}))
        if opt is None:
            return
        cfg = cs._load_config()
        cfg['gui'] = dict(tol=opt['tol'], local=opt['local'], cod=opt['cod'],
                          csd=opt['csd'])
        cs._save_config(cfg)
        outdir = os.path.splitext(gpx)[0] + '_cellsearch'
        log = []

        def logger(*a):
            s = ' '.join(str(x) for x in a)
            log.append(s)
            print(s)
        busy = wx.BusyInfo('Searching databases, please wait...\n(the COD '
                           'query can take a minute or two)')
        wx.BeginBusyCursor()
        try:
            wx.SafeYield()
            query, hits, errors = cs.run_search(
                opt['cand']['cell'], opt['cand']['centring'],
                tol=opt['tol'] / 100., use_cod=opt['cod'],
                use_csd=opt['csd'],
                local=[opt['local']] if opt['local'] else [],
                elements=opt['elements'], outdir=outdir, log=logger)
            csvpath = os.path.join(outdir, 'cellsearch_results.csv')
            cs.write_csv(hits, csvpath, query)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            wx.EndBusyCursor()
            del busy
            wx.MessageBox('Search failed:\n%s' % exc, 'Cell search')
            return
        wx.EndBusyCursor()
        del busy
        print(cs.format_table(hits))
        if opt['icsd']:
            cs.open_icsd(query.cell)
        self._results(cs, query, hits, errors, outdir, csvpath,
                      '\n'.join(log))
