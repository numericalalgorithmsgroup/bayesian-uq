This repository contains code for reproducing the results in the following paper.

```
MCMC for Bayesian uncertainty quantification from time-series data
LNCS, volume 12143, Coputational Science - ICCS 2020, https://doi.org/10.1007/978-3-030-50436-6
https://arxiv.org/abs/2005.14281
```

The code was developed in a Linux environment using the gcc compiler.  It has also been built in a Windows environment using the MSVC cl compiler.  The code tends to compile faster on Linux.  And the Stan examples may not work under MSVC.  On a Windows machine, use of gcc in the Windows Subsystem for Linux environment is recommended.

## Install dependencies for spectral analysis using smMALA

#### Eigen

The linear algebra uses the Eigen library.  For an AD friendly variant of Eigen clone Eigen-AD from here.

[https://gitlab.stce.rwth-aachen.de/stce/eigen-ad](https://gitlab.stce.rwth-aachen.de/stce/eigen-ad)

If you have any questions regarding Eigen-AD, email info@stce.rwth-aachen.de.

#### FFTs

The spectral analysis uses the FFTW library.  Download available from here.

[http://www.fftw.org/download.html](http://www.fftw.org/download.html)

#### dco/c++

Version 3.4.3 or later is needed.  The latest version of dco/c++ 3.4.x for Linux is available here.

[https://www.nag.com/content/downloads-dco-c-dcl6i34ngl](https://www.nag.com/content/downloads-dco-c-dcl6i34ngl)

The examples that use finite differences instead of AD for derivative computation should still work without dco/c++ installed.

## Customize makefile

Copy `makefiles/GNUmakefile.default.inc` to a new file called

```
makefiles/GNUmakefile.[USER].inc
```

where USER is the string returned by `whoami` at the command-line.

Modify the dependency paths according to your setup.

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

By default a central finite difference formula is used to evaluate derivatives for basic scalar types.  To use a forward finite difference formula add the argument **finite_difference_type=ffd**.

To produce a table that summarizes the MCMC samples and includes the N Eff. (Effective Sample Size) diagnostics do,

```
make scalar_type=basic stansummary_smMALA
```

This should reproduce the following table.  Note that the stansummary tool from CmdStan needs to be installed - see below.  Also the last column will be different due to machine and runtime dependent differences in elapsed time.

```
+-----------+----------+----------+----------+---------+-----------+
| name      |       5% |      50% |      95% |   N_Eff |   N_Eff/s |
|-----------+----------+----------+----------+---------+-----------|
| lp__      | 6163.33  | 6166.76  | 6168.51  |     417 |    67.499 |
| omega0_c1 |   78.263 |   80.001 |   82.052 |     282 |    45.737 |
| omega0_c2 |   37.413 |   38.86  |   40.395 |     150 |    24.332 |
| sd_in_c1  |   93.99  |  100.826 |  108.607 |     296 |    47.839 |
| sd_in_c2  |    9.527 |   10.552 |   11.68  |     241 |    38.99  |
| zeta      |    0.169 |    0.194 |    0.219 |     265 |    42.98  |
+-----------+----------+----------+----------+---------+-----------+

```

Producing the table above requires the pandas and tabulate Python modules to be installed, e.g.,

```
python -m pip install -U pandas tabulate
```

## Install CmdStan for use of NUTS

The CmdStan interface to Stan is documented here,

[https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan](https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan)

The Stan software is decribed in this paper,

Carpenter, B, Gelman, A, Hoffman, MD, Lee, D, Goodrich, B, Betancourt, M, Brubaker, MA, Li, P, & Riddell, A. (2017). Stan : A Probabilistic Programming Language.
[https://www.jstatsoft.org/article/view/v076i01](https://www.jstatsoft.org/article/view/v076i01)

This repository comes with a patch for Stan that is needed to make Stan compatible with Eigen-AD.  After installing Stan and adding the Stan path to your makefile, you can  apply this patch as follows,

```
make stan_patch
```

The Stan code in this repository was developed using the v2.20-v2.24 releases of CmdStan (July 2019 - July 2020).  The results table below was generated using v2.24.

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
+-----------+----------+----------+----------+---------+-----------+
| name      |       5% |      50% |      95% |   N_Eff |   N_Eff/s |
|-----------+----------+----------+----------+---------+-----------|
| lp__      | 6163.22  | 6166.7   | 6168.39  |     497 |    58.82  |
| omega0_c1 |   78.308 |   80.306 |   82.487 |     643 |    76.085 |
| omega0_c2 |   37.597 |   39.138 |   40.567 |     658 |    77.881 |
| sd_in_c1  |   94.65  |  101.152 |  109.488 |     384 |    45.537 |
| sd_in_c2  |    9.574 |   10.595 |   11.806 |     546 |    64.645 |
| zeta      |    0.169 |    0.193 |    0.219 |     494 |    58.477 |
+-----------+----------+----------+----------+---------+-----------+
```

Note that N_Eff/s will be machine-dependent because it depends on the sampling execution time.
