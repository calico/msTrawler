# msTrawler
A toolkit for analyzing multiplexed proteomics data
with an emphasis on finding a reliability threshold and accounting
for signal dependent heteroskedactity.  

# Installation

## Option 1: Using `devtools`

```
devtools::install_github("https://github.com/calico/msTrawler")
```

## Option 2: from source 

1. clone or pull repo

```
git clone git@github.com:calico/msTrawler.git
```
2. cd to parent directory
```
cd ~/workspace
```

3. install with `--no-staged-install` flag
```
R CMD INSTALL msTrawler --no-staged-install
```

the `--no-staged-install` option is essential, otherwise will not build.

4. test in RStudio
```
library(msTrawler)

4. test converting proteome discoverer ouput to R dataframe 
```
setwd("/PATH/TO/LIBRARY/FOLDER")
df <- convertPdFile("data/pd_example_export_PSMs.txt", "data/pd_example_export_Proteins.txt")

```

5.  For help using the package please see the tutorial here
...
https://github.com/calico/msTrawler/blob/main/msTrawler_tutorial.pdf
...
