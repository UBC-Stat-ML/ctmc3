# ctmc3: CTMC MCMC Pseudo-Marginal Inference

## Installation

**Step 0**: the package [`mcmcse`](https://cran.r-project.org/package=mcmcse)---used to compute ESSs---indirectly depends on the FFTW library. Follow these steps to install FFTW in

- Ubuntu: `apt-get install libfftw3-dev`
- macOS: `brew install fftw`
- Other: view instructions [here](http://www.fftw.org/install/windows.html)

**Step 1**: install required package `pske` using the provided tarball

``` r
install.packages("pske_0.0.1.tar.gz", repos = NULL, type="source")
```

**Step 2**: install `ctmc3` using the provided tarball

``` r
install.packages("ctmc3_0.0.1.tar.gz", repos = NULL, type="source")
```

## Usage example

``` r
# build a pre-tuned sampler object for a particular experiment
sampler = ctmc3::get_sampler(
  exp_name="SG2019_Sch", # name of the experiment
  reg_ts=TRUE,           # exploit regularity of time series
  gtp_solver="skeletoid" # matrix exponential approximation
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

