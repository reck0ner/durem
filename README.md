# durem

durem is a R library to generate and estimate Duration Relational Event Models (durem).

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
install_github("reck0ner/durem")

#load the package
library(durem)
```

## Usage
```R
# Define the effects for the start model
start_effects <- ~ 1 + 
    remstats::inertia(scaling="std") + 
    remstats::reciprocity(scaling="std")

# Define the effects for the end model
end_effects   <- ~ 1 + 
    remstats::outdegreeSender(scaling="std")

# Define the parameters for each effect
# (order must match the respective formula)
start_params <- c(-7, 0.2, 0.1)
end_params   <- c(-4, -0.2)

# Simulate data
sim <- durem::duremulate(
    start_effects, 
    end_effects, 
    start_params, 
    end_params, 
    num_actors  = 10,  
    num_events  = 500,
    psi_start = 0.5,
    psi_end = 0.5,    
    event_threshold = 1500
)

# Using the previously simulated sim$edgelist
# Estimate the durem model
est <- durem::duremstimate(start_effects = start_effects, 
                           end_effects = end_effects,
                           psi_start = 0.5,
                           psi_end = 0.5,
                           edgelist = sim$edgelist)

summary(est)

```
## Support
```R
#To view all help files in the remulate package
help(package='durem')

```
## Contributing
Pull requests and bug reports are welcome. For major changes, please open an issue first to discuss what you would like to change.

### Acknowledgement
The development of this package was supported by a Vidi Grant (452-17-006) awarded by the Netherlands Organization for Scientific Research (NWO) Grant and an ERC Starting Grant  (758791).
