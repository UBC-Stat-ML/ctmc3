# ctmc3: CTMC MCMC Pseudo-Marginal Inference

## Installation

### Step 1: system dependencies

The package [`mcmcse`](https://cran.r-project.org/package=mcmcse) -- used to compute ESSs -- indirectly depends on the [FFTW library](http://www.fftw.org/). Follow these steps to install FFTW in

- Ubuntu: `apt-get install libfftw3-dev`
- macOS: `brew install fftw`
- Windows: instructions [here](http://www.fftw.org/install/windows.html)

### Step 2: install `ctmc3`

```r
if (!require(remotes)) {
    install.packages('remotes')
}
Sys.setenv(GITHUB_PAT = "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")
remotes::install_github("UBC-Stat-ML/ctmc3")
```


## Usage example

```r
# build a pre-tuned sampler object for a particular experiment
sampler = ctmc3::get_sampler(
  exp_name   = "SG2019_Sch", # name of the experiment
  reg_ts     = TRUE,         # exploit regularity of time series by using RA method
  gtp_solver = "skeletoid"   # matrix exponential approximation
)

# run sampler
assign("FLOPS_COUNTER",0,envir = pske:::glovars) # reset flops counter in pske
res=sampler$run_chain(S=10000L,print_every=100L) # run for 10000 iter, print every 100

# get ess for each parameter and take the mean across all
# report ESS/GFLOP
ess_mean = mean(mcmcse::ess(res$theta))
cat(sprintf("\nEfficiency: %.3f ESS/GFLOP\n",
            ess_mean/(1E-9*pske:::glovars$FLOPS_COUNTER)))

# traceplots and densities using coda utility
coda_ob = coda::mcmc(res$theta)
plot(coda_ob)
```

## Available options

1. `exp_name`: `"SG2019_Sch"` , `"SG2019_LV20"`, `"GHS2017_LV"`, `"GHS2017_SIR"`
2. `reg_ts`: `TRUE`, `FALSE` (but note that `"GHS2017_SIR"` is irregularly sampled)
3. `solver`: `"skeletoid"`, `"unif"`

