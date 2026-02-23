# jLibSearch

## Co-authors

Denice van Herwerden, https://github.com/DvHerwerden
Andriy Rebryk, https://github.com/r3bryk

## Description

`jLibSearch` is used for compound identification using any in-house databases in an appropriate format (see below). 

## Prerequisites

Before using the script, several applications/tools have to be installed:
1. Visual Studio Code; https://code.visualstudio.com/download.
2. Julia Programming Language; https://julialang.org/downloads/.
3. Julia Extension in Visual Studio Code > Extensions > Search “julia” > Press `Install`.

Then, the packages and functions must be loaded as follows:
1. Open the script file, e.g., `jLibSearch.jl`, with Visual Studio Code and wait until Julia environment and extension are loaded.
2. Enable using packages by highlighting `using Pkg` in **line 2** and pressing `“Ctrl + Enter”`.
3. Install packages by highlighting **lines 5-12** and pressing `“Ctrl + Enter”`. This procedure will take some time. When done, you will see a message `julia>` in the terminal field; this message will appear after any operation is completed and Julia is ready to proceed.
4. Mute installation of packages by highlighting **lines 5-12** and pressing `“Ctrl + /”`. This step is needed to avoid repeated installation of the packages every time the script is used.
5. Load/import packages for use by highlighting **lines 15-21** and pressing `“Ctrl + Enter”`.
6. Load all functions by highlighting **lines 31-522** (everything between the headers `Functions for library search & visualization` and `Execution of library search & visualization`) and pressing `“Ctrl + Enter”`.

**NB!** All the steps, except for **steps 3-4**, have to be repeated every time you open the script file.

## How to use the script

### Search library/database:
1. Specify library/database and aligned files as: 
- `pathDB = "C:\\Documents\\Dust_database.xlsx"` in **line 535**.
- `pathFile = "C:\\Documents\\Input_file.xlsx"` in **line 539**.
- Specify the files by highlighting **lines 534-540** and pressing `“Ctrl + Enter”`.

2. Specify in **lines 543-546** the input parameters, such as: 
- RI windows for semi-standard non-polar RI, standard non-polar RI, and AI generated RI (`RIwin = [30, 50, 100]`, respectively). If you do not have any of the values for a particular compound, leave it blank; all the RI values can be left blank, but then the match will be based only on spectral similarity.
- m/z tolerance for matching spectra (`mz_tol = 0.1`)
- Minimum number of highest intensity fragments used from each spectrum; default = 15 (`numfrags = 50`).
- Spectra similarity method: `"DISCO"` (DIstance & Spectrum Correlation Optimization) or `"NDP"` (Normalized Dot Product) (`similarity_method = "DISCO"`)
- Load the parameters by highlighting **lines 543-546** and pressing `“Ctrl + Enter”`.

3. Run library search by highlighting **line 547** (`librarySearch(pathDB, pathFile, RIwin, mz_tol, numfrags, similarity_method)`) and pressing `“Ctrl + Enter”`.

### Create head-to-tail graphs for library vs. experimental (user) spectra:
1. Specify in **lines 553 and 555** the input parameters, such as:
- Path to file with library search results, e.g., `pathFile = "C:\\Documents\\ Input_file_LibSearch.csv”` in **line 553**.
- Which features should be used for library matching in **line 555**: `"all"` uses every entry, a vector of numbers `[1, 2, 6]` uses only those indices, and `collect(1:100)` uses only indices/entries from 1 to 100 (`index = "all"`).
- Load the parameters by highlighting **lines 552-555** and pressing `“Ctrl + Enter”`.

2. Run visualization by highlighting **ine 556** (`libraryVisualization(pathFile, index)`) and pressing `“Ctrl + Enter”`.

## Notes and recommendations

The script takes `CSV` or `Excel` files as input.

The input files must contain at least the following columns to be processed: 
`"Name"`

`"R.I. calc"` or `"Retention Index"` or `"RTI"`

`"ID"`or `"Nr"`

`"Spectrum"`

**NB!** `"Spectrum"` values must be in LECO ChromaTOF, ChromaTOF Sync 2D, or Guineu result file formats: 
- ChromaTOF: `39:4500 52:220 67:9999`

- ChromaTOF Sync 2D: `(39|450.0)(52|22.0)(67:1000.0)`

- Guineu: `[ 39:4500 , 52:220 , 67:1000 ]`

The library/database file must contain at least the following columns to be processed:
`"#", "Name", "Spectrum", “RI”, “RI StdNP”, “RI AI Semi-StdNP”`,

where `“RI”` is semi-standard non-polar RI, `“RI StdNP”` is standard non-polar RI, and `“RI AI Semi-StdNP”` is AI-generated semi-standard non-polar RI. If you do not have any of the RI values for a particular compound, leave the respective column cell blank; all the RI values can be left blank, but then the match will be based only on spectral similarity.

**NB!** `"Spectrum"` values must be in LECO ChromaTOF result file format, i.e. `39:4500 52:220 67:9999`.

## License
[![MIT License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)

Intended for academic and research use.