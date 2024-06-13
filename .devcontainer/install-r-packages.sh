#!/bin/bash
# Get the number of CPU cores
Ncpus=$(cat /proc/cpuinfo | grep processor | wc -l)

# R -e 'install.packages(c("packageName"), version = "1.4.2", repos = "https://cran.r-project.org")', Ncpus = $Ncpus)'