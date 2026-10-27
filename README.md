# BVAR_ — Empirical Macro Toolbox (EMT)

by F. Ferroni and F. Canova

The empirical macro toolbox (EMT) is a free, open-source, end-to-end macro empirics kit, written in MATLAB and compatible with GNU Octave. At its core are Bayesian VAR (BVAR) models—one of the most versatile ways to structure and analyze macroeconomic data, spanning forecasting, scenario design, nowcasting, and shock identification within a single disciplined framework. The package also contains functions and routines to estimate panel VARs, FAVARs, local projections, trend/cycles decomposition and other models. We describe the methodology employed and implementation of the functions. We illustrate the main features with a number of practical examples.

## Getting started

Download or clone the repository and add the two toolbox folders to the MATLAB path:

```matlab
addpath C:/your_folder/BVAR_/bvartools
addpath C:/your_folder/BVAR_/cmintools   % optimisers used to choose the prior hyperparameters
```

Requirements: MATLAB with the Statistics and Optimization toolboxes; the trend-cycle tutorial also uses the Signal Processing toolbox. GNU Octave works too, more slowly and with simpler graphics; Octave versions of the main examples are in `examples/BVAR tutorial Octave`.

Use `bvar_` to estimate a BVAR: since R2020a MATLAB has a built-in function called `bvar`. The toolbox file `bvar.m` is kept only for older MATLAB releases.

## A first estimation

A Bayesian VAR with a Minnesota prior on the monthly US data of Gertler and Karadi (2015), and the responses to a monetary policy shock. Run it from the folder `examples/BVAR tutorial`, where the data are stored.

```matlab
addpath ../../bvartools/ ../../cmintools/
load DataGK                            % US monthly data, 1979m7-2012m6
y    = [logip logcpi gs1 ebp];        % industrial production, CPI, 1-year rate, excess bond premium
lags = 12;

options.priors.name = 'Minnesota';    % Minnesota prior with default hyperparameters
options.K           = 1000;           % posterior draws
options.hor         = 24;             % horizon of the impulse responses
BVAR = bvar_(y, lags, options);

% responses to a monetary policy shock: recursive (Cholesky) identification,
% the 1-year rate is the third variable
opt.varnames    = {'IP', 'CPI', '1-year rate', 'EBP'};
opt.shocksnames = {'Monetary policy'};
plot_irfs_(BVAR.ir_draws(:,:,3,:), opt)
```

`BVAR` holds the posterior draws of the coefficients (`BVAR.Phi_draws`) and of the residual covariance (`BVAR.Sigma_draws`), the impulse responses (`BVAR.ir_draws`, variables x horizons x shocks x draws), the forecasts (`BVAR.forecasts`) and the log marginal likelihood (`BVAR.logmlike`). Other priors, identification schemes and forecasts are switched on through fields of `options`; type `help bvar_` for the list.

## Tutorials

The scripts in `examples/` reproduce the examples of the Hitchhiker's Guide. Run each from its own folder.

`examples/BVAR tutorial`

| Script                        | What it shows                                                                                                                    |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `example_1_classical`         | VAR with a flat prior and bootstrap: lag selection, IRFs and variance decomposition                                              |
| `example_2_minn`              | Minnesota prior: calibrated and estimated hyperparameters, Pandemic Priors, prior for the long run                               |
| `example_3_irf`               | Identification: recursive, sign, narrative, zero, long-run, proxy and heteroskedasticity; variance and historical decompositions |
| `example_4_mfvar`             | Mixed-frequency (monthly-quarterly) VAR and nowcasting                                                                           |
| `example_5_favar`             | Factor-augmented VAR                                                                                                             |
| `example_6_VARX`              | VAR with exogenous variables                                                                                                     |
| `example_7_LP`                | Local projections, classical (OLS and IV) and Bayesian                                                                           |
| `example_8_panels`            | Panels of VARs: average of individual responses and pooled estimator                                                             |
| `example_9_prediction`        | Unconditional and conditional forecasts                                                                                          |
| `example_10_VAR_heterosked`   | Forecasting with heteroskedasticity weights (COVID-19 period)                                                                    |
| `example_11_connectedness`    | Connectedness and spillovers in the cryptocurrency market                                                                        |
| `example_12_bdfm`             | Bayesian static and dynamic factor models                                                                                        |
| `example_13_nonGaussian_bvar` | Non-Gaussian VAR and identification with higher-order moments                                                                    |
| `example_14_set_comparison`   | Set identification with `options.set`, compared with the older options                                                           |
| `dataloading`                 | Loading and ordering a monthly euro area data set                                                                                |

`examples/Trend-Cycle-Dating tutorial`

| Script                     | What it shows                                                |
| -------------------------- | ------------------------------------------------------------ |
| `example_1_cyles`          | Trend-cycle decompositions with different filters            |
| `example_2_dating`         | Dating business-cycle turning points and cyclical statistics |
| `example_3_recession_prob` | Predicting recessions with a probit model                    |

## Documentation

- **The Hitchhiker's Guide** (`HitchhikerGuide_.pdf`) explains the methods, every option and the examples; slides in `HitchhikerGuide_Slides.pdf`.
- **Help in MATLAB**: every function has a help text with inputs, outputs and references, e.g. `help bvar_`, `help plot_irfs_`.
- **Function map**: [dependency graph and list of all functions](https://naffe15.github.io/_pages/function_map.html) (also in [`docs/function_map.html`](docs/function_map.html)). After changing the toolbox, rebuild it in MATLAB with `cd docs; build_function_map`.
- **Video tutorials**: [YouTube channel](https://www.youtube.com/channel/UCDepPX4wbdkIqfg438J0h_g).

## License and citation

The toolbox is free software, released under the GNU General Public License version 3 (see `LICENSE`). If you use it, please cite the Hitchhiker's Guide. Problems and suggestions are welcome as GitHub issues.
