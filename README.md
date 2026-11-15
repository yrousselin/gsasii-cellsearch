# Unit-cell search in COD, CSD and local CIFs for GSAS-II

*[Version française : LISEZMOI.md](LISEZMOI.md)*

A plugin for **GSAS-II** that takes a unit cell from your project, for example the result of indexing a powder pattern, and looks for **known structures with the same lattice** in:

- **COD**, the Crystallography Open Database (online, free);
- **CSD**, the Cambridge Structural Database (if the CCDC software is installed on your computer);
- **any folder of CIF files**, for example structures exported by hand from the ICSD.

The matching CIF files are saved next to your project, so you can import them straight away with *Import → Phase → from CIF file*.

Version 1.5.

## What it does

- Adds the menu entry **Export → Entire project as → (user) Search databases by unit cell** to GSAS-II.
- Lets you choose the cell: an indexing result marked *use* in the *Unit Cells List*, or a phase of the project. Then set the tolerance, optional required elements, and which databases to search.
- **How cells are compared.** Each cell is converted to its primitive lattice and described by a fingerprint that does not depend on the choice of axes: the cell volume and the lengths of the 7 shortest lattice vectors. Two lattices match when all 7 lengths agree within the tolerance (1.5 % by default). This avoids the instabilities of the Niggli reduced cell. Supercells and subcells are not searched.
- Lists the matches, sorted by closeness, in a table and a `.csv` file, and downloads the CIF files.
- Can also open the ICSD web cell-search page with the cell ready to paste.

## Requirements

| What | Why | Where to get it |
|---|---|---|
| **GSAS-II** | hosts the plugin | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| numpy | calculations | included with GSAS-II |
| Internet connection | COD search | |
| **Conograph** *(recommended)* | indexing the powder pattern to obtain the unit cell before the search | [Z-Code downloads](https://z-code-software.com/downloads/) (Windows and macOS) |
| pymysql *(optional)* | much faster COD queries | `pip` (the macOS installer adds it) |
| CCDC software with the CSD Python API *(optional)* | CSD search | [CCDC](https://www.ccdc.cam.ac.uk/) (requires a CSD licence) |

Without the CSD, the COD and local-folder searches still work.

## Installation

### macOS (and Linux)

1. Download this repository: green **Code** button → **Download ZIP**, then double-click the ZIP to unpack it.
2. Open the **Terminal** application. Type `bash `, with a space after it, then drag the file `install_cellsearch_plugin.sh` from the Finder into the Terminal window, and press Enter.
3. The installer:
   - copies the files into `~/.GSASII/exports/`;
   - installs pymysql;
   - looks for the CCDC's Python and says whether the CSD search will work.
4. Restart GSAS-II.

### Windows

1. Download this repository: green **Code** button → **Download ZIP**. Right-click the ZIP → **Extract All**.
2. Double-click **`install_windows.bat`**. It copies the three `.py` files into `C:\Users\<your name>\.GSASII\exports\`, which is the folder where GSAS-II looks for user export plugins.
3. Optional, for faster COD queries: install pymysql with the Python that runs GSAS-II:
   `"<GSAS-II folder>\python.exe" -m pip install pymysql`, where the GSAS-II folder is the one chosen in the GSAS-II installer.
4. Restart GSAS-II.

**If the CSD is not found automatically**, give the path of the CCDC's Python, the one that can `import ccdc`, either:
- with the environment variable `G2_CSD_PYTHON`, or
- with the option `--csd-python PATH` on the command line.

The path found is then remembered in `~/.GSASII/cellsearch.json`.

## Usage

**In GSAS-II**

1. Index your powder pattern in GSAS-II or with Conograph, and mark the chosen cell as *use* in the *Unit Cells List*. You can also open a project that already contains a phase.
2. Open **Export → Entire project as → (user) Search databases by unit cell**.
3. Choose the cell, the tolerance, the required elements (for example `Zn N C H`) and the databases, then confirm.
4. The results table opens. The CIF files are in a folder next to your project.

**From the command line** (with GSAS-II's Python, for example `~/g2main/bin/python` on macOS):

```
python g2_cellsearch.py project.gpx --list                 # list the cells found in the project
python g2_cellsearch.py project.gpx                        # search with the first cell
python g2_cellsearch.py --cell 8.86 10.45 13.42 101.3 99.9 92.7
python g2_cellsearch.py project.gpx --local ~/ICSD_exports --elements "Zn N C H"
python g2_cellsearch.py --help                             # all options
```

Command-line messages are in French.

## Limitations

- Supercells and subcells are not searched.
- The ICSD option only opens the ICSD web page. The address is set for the CNRS access; change it with `--icsd-url`.
- Developed and used on macOS. Windows should work, but has not been fully tested yet. Feedback is welcome.

## Credits and references

- **GSAS-II:** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549.
- **COD:** Crystallography Open Database, <https://www.crystallography.net/>. **CSD:** Cambridge Crystallographic Data Centre.
- **Conograph** (powder indexing): « New functions and graphical user interface attached to powder indexing software », *J. Appl. Cryst.* **2017**, 50, 651–659.
- Written by Yoann Rousselin with the help of Claude (Anthropic) as a coding assistant.

## Licence

MIT. See [LICENSE](LICENSE).

## Author

**Yoann Rousselin**, CNRS research engineer (X-ray crystallography), ICMUB, Dijon (France).
Questions or bug reports: please open an *Issue* on this page.
