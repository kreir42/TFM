## Measurements of (alpha,n) reactions in CNA/HISPANOS
Master's thesis for the [Inter-University Master's Degree in Nuclear Physics](https://master.us.es/fisicanuclear/), carried out in the [CNA](https://cna.us.es/index.php/es) in Seville, Spain, during the 2022-2023 academic year.
The thesis can be compiled from latex/TFM.tex or downloaded from the [latest release](https://github.com/kreir42/TFM/releases).

All files are original, except for "process_ToF.cxx" and "process_filelist.cxx", which were provided by my supervisor and later modified.

## Usage
The data from the experiments, in the form of .root files, are included in this project using [git annex](https://git-annex.branchable.com/).
To analyze the .root files present in the "input" folder and listed in "files", run:
```sh
$ root process_filelist.cxx
```
The result will be several .txt files, intermediary .root files in the "output" folder, and "output.root" with most graphs.

To run with other files with the same format, the code would have to be modified.
