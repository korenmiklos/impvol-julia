# 2018-10-13

## Exploring nominal anchor

```
13-Oct 09:57:50:INFO:root:Nominal world expenditure: 3.726180697148625e15
13-Oct 09:57:50:INFO:root:In the data: [8.33469e6]
```

Data is about final expenditure? No, calculated from total revenue per sector per country.

The two do not match up even in free trade equilibrium.

```
13-Oct 10:21:33:INFO:root:Nominal expenditure in free trade: 1.024434905807753e7
13-Oct 10:21:33:INFO:root:----------------------in the data: 8.334687984343763e6
```

In period 1, they are still the same.

```
13-Oct 10:24:16:INFO:root:Nominal expenditure in free trade: 8.334687984343763e6
13-Oct 10:24:16:INFO:root:----------------------in the data: 8.334687984343763e6
```

Even in period 1, it is lower, moving slowly. (It may be because we calibrate to free-trade revenue, which is higher in this model.)

```
3-Oct 10:30:51:INFO:root:--- Period 1 ---
13-Oct 10:30:53:INFO:root:Total expenditure in the model: 3.1928874924874576e6
13-Oct 10:30:54:INFO:root:Total expenditure in the model: 3.1566047958590924e6
13-Oct 10:30:55:INFO:root:Total expenditure in the model: 3.1615631102681905e6
13-Oct 10:30:56:INFO:root:Total expenditure in the model: 3.1733005969951646e6
13-Oct 10:30:57:INFO:root:Total expenditure in the model: 3.1847511885460145e6
13-Oct 10:30:58:INFO:root:Total expenditure in the model: 3.1942779582839618e6
13-Oct 10:30:59:INFO:root:Total expenditure in the model: 3.2016520691218963e6
13-Oct 10:31:00:INFO:root:Total expenditure in the model: 3.2070644740706724e6
```

Emjs and Rnjs do not match up in the model:

```
13-Oct 10:34:09:INFO:root:Total expenditure in the model: 3.721319483124286e6
13-Oct 10:34:09:INFO:root:------ Middle 1: 0.0002450388818035214
13-Oct 10:34:09:INFO:root:Starting from 1.1561598775110739e-17
13-Oct 10:34:09:INFO:root:-- Outer 12: 0.0008046753470789342
13-Oct 10:34:09:INFO:root:Nominal world expenditure: 3.726180697148625e15
```

No, second line is REAL output. 

# 2018-10-15

## Checking price normalizatiom

Already under free trade eq, US prices are 10 orders of magnitude smaller than we assume in normalization:
```
15-Oct 10:03:20:INFO:root:--- Period 1 ---
15-Oct 10:03:21:INFO:root:US prices: [6.95608e-10, 1.00439e-10, 3.27178e-10, 1.52279e-10, 1.60343e-10]
15-Oct 10:03:21:INFO:root:in the data: [1.0, 1.0, 1.0, 1.0, 1.0]
```

Fixed a bug on line 380 of `calibrate_params`, but still similar magnitudes.

Wages are also off, but not by the same magnitude:
```
15-Oct 10:12:32:INFO:root:Data wages: [1.12031e6, 1.12031e6, 1.12031e6, 1.12031e6, 1.12031e6]
WARNING: using ImpvolEquilibrium.array_transpose in module Main conflicts with an existing identifier.
15-Oct 10:12:37:INFO:root:--- Period 1 ---
15-Oct 10:12:40:INFO:root:US wages: [2.02839e6, 3.48099e5, 6.17292e5, 6.92271e5, 740738.0]
```

Apparently, `relative_A` is not 1.0 for US:

```
15-Oct 10:58:35:INFO:root:[7165.76 7131.23 6239.92 6359.79 6353.39 6465.62 6446.54 6460.09 5606.42 
```

Bug: `calculate_A` already calculates level of A, not relative to US. It was only the price that needed a level pin down. My calculation and Balazs's calculation of US productivity is the same.
```
15-Oct 11:20:04:INFO:root:[7165.76, 580.413, 4490.71, 2305.41, 2757.17, 1644.78, 1512.43, 1492.24, 8984.6, 112.126, 1893.54, 3111.52, 5174.74, 1665.38, 3379.38, 3511.17, 1390.29, 3924.83, 2698.91, 7473.85,1340.21, 3155.8, 3718.4]
15-Oct 11:20:04:INFO:root:[7165.76, 580.413, 4490.71, 2305.41, 2757.17, 1766.57, 1522.09, 1492.24, 8984.6, 112.126, 1893.54, 3150.18, 5174.74, 1665.38, 3379.38, 3618.48, 1464.01, 3924.83, 2839.14, 7598.12,1399.67, 3309.09, 3866.28]
```

Magnitude of US prices now ok, but still not exact level. Why?
```
15-Oct 11:25:01:INFO:root:US prices: [0.451145, 0.451288, 0.451569, 0.450765, 0.450904]
15-Oct 11:25:01:INFO:root:in the data: [1.0, 1.0, 1.0, 1.0, 1.0]
```

Model captures trade shares pretty well, so this is not the cause for price level differences.
```
15-Oct 11:58:23:INFO:root:-- Outer 12: 0.0008042939716096958
15-Oct 11:58:23:INFO:root:Model trade shares: [0.97146, 0.942892, 0.995727, 0.95153, 0.97311]
15-Oct 11:58:23:INFO:root:Data trade shares: [0.970499, 0.941882, 0.994321, 0.946233, 0.966882]
15-Oct 11:58:33:INFO:root:US prices: [0.451145, 0.451288, 0.451569, 0.450765, 0.450904]
15-Oct 11:58:33:INFO:root:in the data: [1.0, 1.0, 1.0, 1.0, 1.0]
```

# 2018-10-16

## Simulate model in sterile environment

Opened an [issue](https://github.com/korenmiklos/impvol-julia/issues/18). The formula y/psi is only a good approximation of wages if labor cannot readjust.

Calibrated a minimal, 2x3 economy.
```
WARNING: replacing module CalibrateParameters
WARNING: using ImpvolEquilibrium.array_transpose in module CalibrateParameters conflicts with an existing identifier.
2×3 Array{Float64,2}:
 0.899504  1.11172  1.03593
 1.0       1.0      1.0
2×3 Array{Float64,2}:
 13.0891  10.2124  11.8764
 11.1547  11.878   12.195
16-Oct 11:21:45:INFO:root:US wage rate: 3.0
Test Passed
```

Prices have diverged from original 1.0:
```
16-Oct 16:44:33:INFO:root:Model trade shares: [0.725587, 0.910351, 1.0]
16-Oct 16:44:33:INFO:root:Data trade shares: [0.7, 0.9, 1.0]
16-Oct 16:44:43:INFO:root:US prices: [2.61974, 2.60394, 2.59458]
16-Oct 16:44:43:INFO:root:in the data: [1.0, 1.0, 1.0]
16-Oct 16:44:43:INFO:root:Nominal world expenditure: 10.549620676278591
```

Equation (15) in "paper November 8 2017" has Pnt, not p_sectoral. This seems wrong in line 323. Redoing now using eq (15) in "algorithm.pdf".

Fixed calculation of `A`. We are closer to the prices, not not fully there:
```
16-Oct 18:05:16:INFO:root:Model trade shares: [0.725587, 0.910351, 1.0]
16-Oct 18:05:16:INFO:root:Data trade shares: [0.7, 0.9, 1.0]
16-Oct 18:05:26:INFO:root:US prices: [0.873248, 0.86798, 0.86486]
16-Oct 18:05:26:INFO:root:in the data: [1.0, 1.0, 1.0]
16-Oct 18:05:26:INFO:root:Nominal world expenditure: 10.549620676278582
16-Oct 18:05:26:INFO:root:--------------In the data: [12.0]
```

# 2018-10-18

## Code choice of numeraire

Our current choice of numeraire is that we start the free-trade initialization from a given nominal world expenditure. However, after trade costs are reclibrated, world nominal expenditure is not an exact match. (See above.) Propose: Set US aggregate price index to 1.0 in each time period.

```
18-Oct 10:13:54:INFO:root:Model trade shares: [0.725587, 0.910351, 1.0]
18-Oct 10:13:54:INFO:root:Data trade shares: [0.7, 0.9, 1.0]
18-Oct 10:14:04:INFO:root:US prices: [1.00642, 1.00035, 0.99675]
18-Oct 10:14:04:INFO:root:in the data: [1.0, 1.0, 1.0]
18-Oct 10:14:04:INFO:root:Nominal world expenditure: 12.158414061090557
18-Oct 10:14:04:INFO:root:--------------In the data: [12.0]
``` 

## Dependencies and name space

`scenario.jl` loads `init_parameters.jl`. This defines a module `Environment` and within that module, loads `config.jl`. This loads `calibrate_params.jl`. This defines a module `CalibrateParameters`. (The absolute name of this module seems to be `Environment.CalibrateParameters`.) This loads `equilibriu.jl`, defining the module `ImpvolEquilibrium` (`Environment.CalibrateParameters.ImpvolEquilibrium`).

New workflow:
`scenario.jl` loads `equilibrium.jl`, defining module `ImpvolEquilibrium`. It also loads `init_parameters.jl`, which loads `config.jl`. The latter loads `calibrate_params.jl`, which defines module `CalibrateParameters`. init and config define no other modules, `parameters` is simply a global. At this point, there are two modules defined, `ImpvolEquilibrium` and `CalibrateParameters`. `scenario.jl` further loads `change_parameters.jl`, which has to change the global `parameters` and can access the functions exported by the two modules.


## Nominal prices in simulation

Prices move around too much over time. Already in 1973, US sectoral prices are:
```
24-element Array{Float64,1}:
 1.00258
 1.24319
 1.39467
 1.19452
 1.23191
 1.15468
 1.30421
 1.21179
 1.25413
 1.27171
 1.12897
 1.15808
 1.19688
 1.14249
 1.2054
 1.19571
 1.08363
 1.1775
 1.08566
 1.16241
 1.11747
 1.12605
 1.18582
 0.859164
 ```
 (Remember that 1972 prices are normalized to 1.0). In 2007, they are
 ```
 24-element Array{Float64,1}:
  3.62405
 13.8792
 41.2767
 20.1324
 27.1493
 27.3631
 25.1278
 34.7121
 46.1475
 22.5707
 15.9228
 36.1945
 41.1135
 23.7509
 43.9461
 46.6444
  5.08522
  7.71534
  6.52232
  7.26517
 44.5359
 45.7468
 38.4647
  0.280734
 ```

 There is something off with nontradable productivity/price dynamics.