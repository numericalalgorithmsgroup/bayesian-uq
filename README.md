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

Version 3.4.3 or later is needed.  As of mid July 2020, this is not publicly available, but early access can be arranged on request.  The examples that use finite differences instead of AD for derivative computation should still work without dco/c++ installed.

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

To produce a table that summarizes the MCMC samples and includes the N Eff. (Effective Sample Size) diagnostics do,

```
make scalar_type=basic stansummary_smMALA
```

This should reproduce the following table.  Note the last column will be different due to machine and runtime dependent differences in elapsed time.

```
+----+-----------+----------+----------+----------+---------+-----------+
|    | name      |       5% |      50% |      95% |   N_Eff |   N_Eff/s |
|----+-----------+----------+----------+----------+---------+-----------|
|  0 | lp__      | 6163.33  | 6166.76  | 6168.51  |     417 |    79.546 |
|  1 | omega0_c1 |   78.263 |   80.001 |   82.052 |     282 |    53.904 |
|  2 | omega0_c2 |   37.413 |   38.86  |   40.395 |     150 |    28.678 |
|  3 | sd_in_c1  |   93.99  |  100.825 |  108.607 |     296 |    56.381 |
|  4 | sd_in_c2  |    9.527 |   10.552 |   11.68  |     241 |    45.953 |
|  5 | zeta      |    0.169 |    0.194 |    0.219 |     265 |    50.655 |
+----+-----------+----------+----------+----------+---------+-----------+

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

The Stan code in this repository was developed using the v2.20.0 release of CmdStan (July 2019).

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
| lp__      | 6163.32  | 6166.69  | 6168.4   |     448 |    93.091 |
| omega0_c1 |   78.338 |   80.288 |   82.385 |     563 |   117.019 |
| omega0_c2 |   37.702 |   39.083 |   40.642 |     605 |   125.843 |
| sd_in_c1  |   94.265 |  101.422 |  109.488 |     366 |    76.177 |
| sd_in_c2  |    9.596 |   10.598 |   11.729 |     520 |   108.164 |
| zeta      |    0.168 |    0.194 |    0.219 |     373 |    77.531 |
+-----------+----------+----------+----------+---------+-----------+
```

Note that N_Eff/s will be machine-dependent because it depends on the sampling execution time.
