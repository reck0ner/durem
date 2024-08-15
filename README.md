# drem

drem is a R library to generate and estimate Duration Relational Event Models (DREM).

### Programming Languages
The package contains code written in:
* R (>= 4.0.0)
* Rcpp (>= 1.0.4.6) and RcppArmadillo (>= 0.9.860.2.0)
* C++11 (Compiler Version: GCC-8.1.0)
	
## Installation

To install the package in R using `devtools`:

```R
library(devtools)
install.package(remify)
install.package(remstats)
install_github("reck0ner/drem")

#load the package
library(drem)
```

## Usage
```R
effects <- ~ baseline(-4) + inertia(0.01) + reciprocity(-0.04) + itp(0.01,scaling="std")

dremulateTie(effects, actors = 1:25, time = 20, events = 500, initial = 200)

```
## Support
```R
#To view all help files in the remulate package
help(package='drem')

#To view available effects
help('dremulateTieEffects')

```
## Contributing
Pull requests and bug reports are welcome. For major changes, please open an issue first to discuss what you would like to change.

### Acknowledgement
The development of this package was supported by a Vidi Grant (452-17-006) awarded by the Netherlands Organization for Scientific Research (NWO) Grant and an ERC Starting Grant  (758791).
