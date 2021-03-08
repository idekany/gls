# Generalized Lomb-Scargle Periodogram (GLS)
---

This is a Fortran90 parallelized implementation of the GLS algorithm of [Zechmeister & Kürster](https://www.aanda.org/articles/aa/pdf/2009/11/aa11296-08.pdf) using OpenMP.

## Installation
---

Simply issue the 'make' command in the directory containing all the files of this package to compile the source code (a fortran compiler such as gfortran is required).

## Usage
---

The program can be run with the following syntax:

`gls <inputfile>`

where inputfile is an ASCII text file containing the time series to be transformed, and must have 2 or 3 columns separated by whitespace(s): 
`time, measurement, [measurement error]`

See also the output
`gls -help`
for more information.

An example input file called `test.lc` is included in this package.

The algorithm's parameters can be specified in the `gls.par` parameter file, and must be located in the directory where the `gls` exacutable is invoked.
A default parameter file is included in this package, please keep its format when editing. Please see the comments in the parameter file for more details.

Upon execution, the program writes the following results to the STDOUT:

`inputfile N freq power fmax nfreq FAP [FAPMC] 1/T`

where:
`N`: number of input data points
`freq`: frequency of the highest GLS `power` in the periodogram
`fmax`: maximum sampled frequency
`nfreq`: number of frequency steps
`FAP`: analytical estimate of the false-alarm pobability
`FAPMC`: estimate of the false-alarm pobability by Monte-Carlo simulation (if set in `gls.par`)
`1/T`: reciprocal of the total length of the time series.


