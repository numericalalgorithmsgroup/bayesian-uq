This repository contains code for reproducing the results in the following paper.

```
MCMC for Bayesian uncertainty quantification from time-series data
LNCS, ICCS 2020 conference proceedings (accepted)
Philip Maybank, Patrick Peltzer, Uwe Naumann, Ingo Bojak
https://arxiv.org/abs/2005.14281
```

## Install dependencies for spectral analysis using smMALA

#### Eigen

The linear algebra uses the Eigen library.  For an AD friendly variant of Eigen clone Eigen-AD from here.

[https://gitlab.stce.rwth-aachen.de/stce/eigen-ad](https://gitlab.stce.rwth-aachen.de/stce/eigen-ad)

If you have any questions regarding Eigen-AD, email info@stce.rwth-aachen.de.

#### FFTs

The spectral analysis uses the FFTW library.  Download available from here.

[http://www.fftw.org/download.html](http://www.fftw.org/download.html)

#### dco/c++

Version 3.4.3 or later is needed.  As of early June 2020, this is still not publicly available.  The examples that use finite differences instead of AD for derivative computation should still work without dco/c++ installed.

## Customize makefile

Create a file called

```
makefiles/GNUmakefile.[USER].inc
```

where USER is the string returned by `whoami` at the command-line.

Add include paths to dependencies in the makefile.  For an example see, `makefiles/GNUmakefile.philipm.inc`.

## Generate smMALA results

Results for the harmonic oscillator example with the smMALA MCMC sampler can be generated as follows,

```
make smMALA-harmonic-oscillator.r
```

To compare test output with the expected standard output, add the argument **STORERES=1**.  E.g.,

```
make STORERES=1 smMALA-harmonic-oscillator.r
```

To use the basic C scalar type (i.e. double) instead of dco/c++ types, add the argument **scalar_type=basic**.

To produce a table that summarizes the MCMC samples (similar to Table 1 in the paper) do,

```
make scalar_type=basic pysummary_smMALA
```

This should reproduce the following table,

```
+-----------+---------+---------+---------+
|           |   0.025 |     0.5 |   0.975 |
|-----------+---------+---------+---------|
| omega0_c1 |  77.864 |  80.001 |  82.422 |
| omega0_c2 |  37.075 |  38.86  |  40.821 |
| sd_in_c1  |  92.937 | 100.826 | 110.169 |
| sd_in_c2  |   9.41  |  10.552 |  11.888 |
| zeta      |   0.167 |   0.194 |   0.225 |
+-----------+---------+---------+---------+

```

This requires the pandas and tabulate Python modules to be installed, e.g.,

```
python -m pip install -U pandas tabulate
```

## Install CmdStan for use of NUTS

The CmdStan interface to Stan is documented here,

[https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan](https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan)

The Stan software is decribed in this paper,

Carpenter, B, Gelman, A, Hoffman, MD, Lee, D, Goodrich, B, Betancourt, M, Brubaker, MA, Li, P, & Riddell, A. (2017). Stan : A Probabilistic Programming Language.
[https://www.jstatsoft.org/article/view/v076i01](https://www.jstatsoft.org/article/view/v076i01)

The repository comes with a patch for Stan that is needed to make Stan compatible with Eigen-AD.  After installing Stan and adding the Stan path to your makefile, you can  apply this patch as follows,

```
make stan_patch
```

## Generate NUTS results

Results for the harmonic oscillator example with the NUTS MCMC sampler (with dco types used in derivative computation) can be generated as follows,

```
make stansummary_nuts
```

For analagous results where finite differences are used for the derivative computation do,

```
make scalar_type=basic stansummary_nuts
```

This should approximately reproduce the following table,

```
                 Mean     MCSE   StdDev     5%    50%    95%    N_Eff  N_Eff/s    R_hat
lp__             6155  9.0e-02  1.7e+00   6152   6155   6157  3.6e+02  7.7e+01  1.0e+00
accept_stat__    0.91  2.8e-03  9.9e-02   0.71   0.95    1.0  1.2e+03  2.6e+02  1.0e+00
stepsize__       0.46      nan  2.8e-16   0.46   0.46   0.46      nan      nan      nan
treedepth__       2.8  1.9e-02  5.6e-01    2.0    3.0    4.0  8.1e+02  1.7e+02  1.0e+00
n_leapfrog__      8.3  1.3e-01  3.8e+00    3.0    7.0     15  8.1e+02  1.7e+02  1.0e+00
divergent__      0.00      nan  0.0e+00   0.00   0.00   0.00      nan      nan      nan
energy__        -6152  1.3e-01  2.4e+00  -6155  -6153  -6148  3.1e+02  6.6e+01  1.0e+00
theta[1]          4.4  5.7e-04  1.5e-02    4.4    4.4    4.4  6.6e+02  1.4e+02  1.0e+00
theta[2]          3.7  9.1e-04  2.2e-02    3.6    3.7    3.7  5.8e+02  1.2e+02  1.0e+00
theta[3]          4.6  2.0e-03  4.2e-02    4.6    4.6    4.7  4.4e+02  9.3e+01  1.0e+00
theta[4]          2.4  2.7e-03  6.0e-02    2.3    2.4    2.5  4.9e+02  1.0e+02  1.0e+00
theta[5]         -1.6  3.7e-03  7.5e-02   -1.8   -1.6   -1.5  4.1e+02  8.6e+01  1.0e+00
```

Note that N_Eff/s will be machine-dependent because it depends on the sampling execution time.