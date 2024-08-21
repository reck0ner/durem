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
start_effects <- ~ 1 + remstats::inertia(scaling="std") + remstats::reciprocity(scaling="std")
stop_effects <- ~ 1 + remstats::outdegreeSender(scaling="std")

start_params <- c(-7,0.2,0.1)
stop_params <-  c(-4, -0.2)

drem::dremulateTie(start_effects, stop_effects, start_params, stop_params, num_actors = 10, num_events = 1000, stop_threshold = 1500)

```
## Support
```R
#To view all help files in the remulate package
help(package='drem')

```
## Contributing
Pull requests and bug reports are welcome. For major changes, please open an issue first to discuss what you would like to change.

### Acknowledgement
The development of this package was supported by a Vidi Grant (452-17-006) awarded by the Netherlands Organization for Scientific Research (NWO) Grant and an ERC Starting Grant  (758791).
