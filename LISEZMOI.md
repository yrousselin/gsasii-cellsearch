# Recherche par maille dans la COD, la CSD et des CIF locaux pour GSAS-II

*[English version: README.md](README.md)*

Extension pour **GSAS-II**. Elle prend une maille de votre projet, par exemple le résultat de l'indexation d'un diagramme de poudre, et cherche des **structures connues ayant le même réseau** dans :

- la **COD**, Crystallography Open Database (en ligne, gratuite) ;
- la **CSD**, Cambridge Structural Database (si les logiciels du CCDC sont installés sur votre ordinateur) ;
- **n'importe quel dossier de fichiers CIF**, par exemple des structures exportées à la main depuis l'ICSD.

Les CIF trouvés sont enregistrés à côté de votre projet : vous pouvez les importer directement avec *Import → Phase → from CIF file*.

Version 1.5.

## Ce que fait l'extension

- Ajoute à GSAS-II le menu **Export → Entire project as → (user) Search databases by unit cell**. La fenêtre et les messages sont en anglais.
- Vous laisse choisir la maille : un résultat d'indexation marqué *use* dans la *Unit Cells List*, ou une phase du projet. Vous réglez ensuite la tolérance, les éléments exigés (facultatif) et les bases à interroger.
- **Comment les mailles sont comparées.** Chaque maille est ramenée à son réseau primitif et décrite par une empreinte qui ne dépend pas du choix des axes : le volume de la maille et la longueur des 7 plus courts vecteurs du réseau. Deux réseaux correspondent quand les 7 longueurs concordent à la tolérance près (1,5 % par défaut). Cette méthode évite les instabilités de la maille réduite de Niggli. Les surstructures et sous-mailles ne sont pas recherchées.
- Liste les correspondances, classées de la plus proche à la plus éloignée, dans un tableau et un fichier `.csv`, puis télécharge les CIF.
- Peut aussi ouvrir la page de recherche par maille de l'ICSD, avec la maille prête à coller.

## Prérequis

| Quoi | Pourquoi | Où le trouver |
|---|---|---|
| **GSAS-II** | accueille l'extension | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| numpy | calculs | fourni avec GSAS-II |
| Connexion Internet | recherche dans la COD | |
| **Conograph** *(recommandé)* | indexation du diagramme de poudre pour obtenir la maille avant la recherche, idéalement via l'extension [gsasii-autoindex](https://github.com/yrousselin/gsasii-autoindex) | [téléchargements Z-Code](https://z-code-software.com/downloads/) (Windows et macOS) |
| pymysql *(facultatif)* | interrogation de la COD beaucoup plus rapide | `pip` (l'installeur macOS l'ajoute) |
| Logiciels du CCDC avec la CSD Python API *(facultatif)* | recherche dans la CSD | [CCDC](https://www.ccdc.cam.ac.uk/) (licence CSD nécessaire) |

Sans la CSD, les recherches dans la COD et dans les dossiers locaux fonctionnent quand même.

## Installation

### macOS (et Linux)

1. Téléchargez ce dépôt : bouton vert **Code** → **Download ZIP**, puis double-cliquez sur le ZIP pour le décompresser.
2. Ouvrez l'application **Terminal**. Tapez `bash ` (avec un espace après), faites glisser le fichier `install_cellsearch_plugin.sh` depuis le Finder dans la fenêtre du Terminal, puis appuyez sur Entrée.
3. L'installeur :
   - copie les fichiers dans `~/.GSASII/exports/` ;
   - installe pymysql ;
   - cherche le Python du CCDC et indique si la recherche dans la CSD fonctionnera.
4. Redémarrez GSAS-II.

### Windows

1. Téléchargez ce dépôt : bouton vert **Code** → **Download ZIP**. Clic droit sur le ZIP → **Extraire tout**.
2. Double-cliquez sur **`install_windows.bat`**. Il copie les trois fichiers `.py` dans `C:\Users\<votre nom>\.GSASII\exports\`, le dossier où GSAS-II cherche les extensions d'export de l'utilisateur.
3. Facultatif, pour accélérer la COD : installez pymysql avec le Python de GSAS-II :
   `"<dossier de GSAS-II>\python.exe" -m pip install pymysql`, le dossier de GSAS-II étant celui choisi lors de son installation.
4. Redémarrez GSAS-II.

**Si la CSD n'est pas trouvée automatiquement**, indiquez le chemin du Python du CCDC, celui qui sait faire `import ccdc`, au choix :
- avec la variable d'environnement `G2_CSD_PYTHON` ;
- avec l'option `--csd-python CHEMIN` en ligne de commande.

Le chemin trouvé est ensuite mémorisé dans `~/.GSASII/cellsearch.json`.

## Utilisation

**Dans GSAS-II**

1. Indexez votre diagramme de poudre, par exemple avec l'extension [gsasii-autoindex](https://github.com/yrousselin/gsasii-autoindex) (Conograph et autres méthodes), puis marquez la maille retenue comme *use* dans la *Unit Cells List*. Vous pouvez aussi ouvrir un projet qui contient déjà une phase.
2. Ouvrez **Export → Entire project as → (user) Search databases by unit cell**.
3. Choisissez la maille, la tolérance, les éléments exigés (par exemple `Zn N C H`) et les bases, puis validez.
4. Le tableau des résultats s'ouvre. Les CIF sont dans un dossier à côté de votre projet.

**En ligne de commande** (avec le Python de GSAS-II, par exemple `~/g2main/bin/python` sur macOS) :

```
python g2_cellsearch.py projet.gpx --list                  # liste les mailles du projet
python g2_cellsearch.py projet.gpx                         # recherche avec la première maille
python g2_cellsearch.py --cell 8.86 10.45 13.42 101.3 99.9 92.7
python g2_cellsearch.py projet.gpx --local ~/ICSD_exports --elements "Zn N C H"
python g2_cellsearch.py --help                             # toutes les options
```

## Limites

- Les surstructures et sous-mailles ne sont pas recherchées.
- L'option ICSD ouvre seulement la page web de l'ICSD. L'adresse est réglée pour l'accès CNRS ; modifiez-la avec `--icsd-url`.
- Développé et utilisé sur macOS. Windows devrait fonctionner, mais n'a pas encore été testé complètement. Vos retours sont les bienvenus.

## Crédits et références

- **GSAS-II :** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549.
- **COD :** Crystallography Open Database, <https://www.crystallography.net/>. **CSD :** Cambridge Crystallographic Data Centre.
- **Conograph** (indexation de poudre) : « New functions and graphical user interface attached to powder indexing software », *J. Appl. Cryst.* **2017**, 50, 651–659.
- Écrit par Yoann Rousselin avec l'aide de Claude (Anthropic) comme assistant de programmation.

## Licence

MIT. Voir [LICENSE](LICENSE).

## Auteur

**Yoann Rousselin**, ingénieur de recherche CNRS (cristallographie), ICMUB, Dijon.
Questions ou signalement d'un problème : ouvrez une *Issue* sur cette page.
