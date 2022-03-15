# ctmc3: CTMC MCMC Pseudo-Marginal Inference

An R package for performing inference on discretely observed continuous-time Markov chains using a pseudo-marginal Markov chain Monte Carlo approach.


## Installation

### Step 1: system dependencies

#### FFTW (only needed if building dependencies from source)

The package [`mcmcse`](https://cran.r-project.org/package=mcmcse) -- used to compute ESSs -- indirectly depends on the [FFTW library](http://www.fftw.org/). Follow these steps to install FFTW in

- Ubuntu: `apt-get install libfftw3-dev`
- macOS: `brew install fftw`
- Windows: instructions [here](http://www.fftw.org/install/windows.html)

### Step 2: install `ctmc3`

```r
if (!require(remotes) || packageVersion("remotes") < package_version("2.4.2")) {
    install.packages("remotes")
}
Sys.setenv(GITHUB_PAT = "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")
remotes::install_github("UBC-Stat-ML/ctmc3@main")
```


## Usage example

Build a pre-tuned sampler object for a particular experiment
```r
sampler = ctmc3::get_sampler(
  exp_name   = "SG2019_Sch_log", # experiment: Schloegl data with sampler in log-space
  reg_ts     = TRUE,             # exploit regularity of time series by using RA method
  gtp_solver = "skeletoid"       # matrix exponential approximation
)
```
Run the sampler and measure ESS per billion matrix operations
```r
pske::reset_ops_counter() # set ops counter to 0
res=sampler$run_chain(S = 10000L, print_every = 100L) # run for 10000 iter, print every 100

# get ess for each parameter and take the mean across all
# report ESS/GMOs
ess_mean = mean(mcmcse::ess(exp(res$theta))) # compute average ESSs (need to invert log transform)
n_ops    = pske::get_ops_counter()
cat(sprintf("\nEfficiency: %.3f ESS/GMOs\n", ess_mean/(1E-9*n_ops)))
```
Compare this number to the results in the Experiments section of the paper. 

Graph traceplots and densities using the coda utility
```r
coda_ob = coda::mcmc(exp(res$theta))
plot(coda_ob)
```

### Available options for `get_sampler`

1. `exp_name`:
    - Sampler in theta space:
        - `"SG2019_Sch"`
        - `"SG2019_LV20"`
        - `"GHS2017_LV"`
        - `"GHS2017_SIR"`
    - Sampler in log(theta) space:
        - `"SG2019_Sch_log"`
        - `"SG2019_LV20_log"`
        - `"GHS2017_LV_log"`
2. `reg_ts`:
    - `TRUE` : uses RA
    - `FALSE`: uses IA
3. `gtp_solver` (matrix exponential algorithm):
    - `"skeletoid"`
    - `"unif"`

### Tuning (currently unsupported on Windows)

One may also tune the sampler from scratch (takes a while)
```r
sampler = ctmc3::get_sampler(
  exp_name   = "SG2019_Sch_log", # experiment: Schloegl data with sampler in log-space
  reg_ts     = TRUE,             # exploit regularity of time series by using RA method
  gtp_solver = "skeletoid",      # matrix exponential approximation
  tuned_pars = FALSE             # do not used pre-tuned parameters
)
ctmc3::tune_sampler(sampler,n_cores = 4L)
```
It will write the results to the current working directory.





