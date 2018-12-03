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

Nontradable productivity is calibrate to go up by a huge amount (see last row)
```
 julia> [A1[1,end,:,1] A36[1,end,:,1]]
24×2 Array{Float64,2}:
 . . .
      3.80674e5  871490.0
      4.80477e5       8.99647e5
      2.86229e8       3.59988e11
```

Nominal US value added increases gradually:
```
julia> sum(data["va"][1,end,:,:], 1)
1×36 Array{Float64,2}:
 1.12031e6  1.25371e6  1.35986e6  1.4886e6  1.6613e6  …  1.0998e7  1.1712e7  1.2409e7  1.2997e7
```
Estimated sectoral prices are also smooth:
```
julia> parameters[:p_sectoral][1,end,end-2:end,:]
3×36 Array{Float64,2}:
 1.0  1.02114  1.11355  1.25773  1.34228  1.43614  …  2.88145  2.93181  2.97346  2.99878  3.03941
 1.0  1.06523  1.2077   1.33567  1.39121  1.42458     3.08848  3.13912  3.23314  3.32482  3.4055
 1.0  1.02928  1.07685  1.21017  1.30383  1.38651     4.67841  4.79941  4.90868  4.56207  4.00687
```

Fixed an apparent bug on line 309, but now prices and trade shares are both off
```
18-Oct 13:07:50:INFO:root:Model trade shares: [0.839914, 0.81171, 0.979848]
18-Oct 13:07:50:INFO:root:Data trade shares: [0.970499, 0.941882, 0.994321]
18-Oct 13:08:00:INFO:root:US prices: [3.46751, 13.2969, 41.1655]
```
Nontradable productiviy still blows up
```
  3866.28    3098.45
 59285.3        1.67605e7
```
Input prices are off!
```
18-Oct 13:16:31:INFO:root:Calibrated input prices:
18-Oct 13:16:31:INFO:root:[1.0, 1.0, 1.0]
18-Oct 13:16:31:INFO:root:[1.26707, 1.18481, 254.632]
```

It is `gamma'` to calculate the sectoral price index, not `gamma`. Now productivity growth looks normal
```
julia> (A36 ./ A1)[1,end,:,1]
24-element Array{Float64,1}:
 1.50131
 1.48263
 1.63589
 1.74895
 1.53929
 1.17424
 1.54995
 1.24335
 1.42291
 0.650903
 1.15389
 ⋮
 1.27075
 1.36602
 1.25941
 3.98871
 8.16857
 5.21963
 7.53213
 1.30153
 1.63239
 1.58981
 1.87969
```

## Back to adjustment loop

As long as wage gap is small, algo seems to converge:
```
18-Oct 17:02:16:INFO:root:Mean absolute wage gap: 0.33467752623816616
18-Oct 17:02:16:INFO:root:Gradient: 0.0008632666889728881
18-Oct 17:02:16:INFO:root:Difference: -1.7866160602952874e-5
18-Oct 17:02:16:INFO:root:Proportional increase: -0.3995673719477588
18-Oct 17:02:16:INFO:root:---- Adjustment 11: 0.0008632666889728881
```

But then wage gap blows up:
```
18-Oct 17:02:20:INFO:root:Mean absolute wage gap: 0.06710612235301563
18-Oct 17:02:20:INFO:root:Starting from 0.0011696159608411702
18-Oct 17:02:20:INFO:root:------ Middle 1: 0.0018415561136925984
18-Oct 17:02:20:INFO:root:------ Middle 2: 0.0004696180190675084
18-Oct 17:02:20:INFO:root:Mean absolute wage gap: 5.14264715222008
18-Oct 17:02:20:INFO:root:Gradient: 1.2187102768125433
18-Oct 17:02:20:INFO:root:Difference: -0.0005269412387152173
18-Oct 17:02:20:INFO:root:Proportional increase: -5.9130287546716004e-6
18-Oct 17:02:20:INFO:root:---- Adjustment 1: 1.2187102768125433
18-Oct 17:02:21:INFO:root:------ Middle 1: 0.2154776298323131
18-Oct 17:02:22:INFO:root:------ Middle 2: 0.040436994183319114
18-Oct 17:02:22:INFO:root:------ Middle 3: 0.010356126020335758
18-Oct 17:02:22:INFO:root:------ Middle 4: 0.006033075579390297
18-Oct 17:02:22:INFO:root:------ Middle 5: 0.001665137520795596
18-Oct 17:02:22:INFO:root:------ Middle 6: 0.0012098033085186535
18-Oct 17:02:22:INFO:root:------ Middle 7: 0.0005212181457145956
18-Oct 17:02:22:INFO:root:Mean absolute wage gap: 70825.42096895793
18-Oct 17:02:22:INFO:root:Gradient: 8957.330061518034
```

# 2018-10-19

## Testing the adjustment loop in the test environment

In `test.jl`, even with `one_over_rho=0.1`, the adjustment loop converges. Curiously, the gradient decreases even after a negative difference.
```
19-Oct 09:30:52:INFO:root:---- Adjustment 1: 0.0023589792249189443
19-Oct 09:30:53:INFO:root:------ Middle 1: 0.00011191856902462655
19-Oct 09:30:53:INFO:root:Mean absolute wage gap: 0.020973740549026548
19-Oct 09:30:53:INFO:root:Gradient: 0.002074178005398289
19-Oct 09:30:53:INFO:root:Difference: -2.618292490889863e-5
19-Oct 09:30:53:INFO:root:Proportional increase: -10.143197001291403
19-Oct 09:30:53:INFO:root:---- Adjustment 2: 0.002074178005398289
19-Oct 09:30:54:INFO:root:------ Middle 1: 9.194374672958078e-5
19-Oct 09:30:54:INFO:root:Mean absolute wage gap: 0.020559640891600024
19-Oct 09:30:54:INFO:root:Gradient: 0.0018238077156721327
19-Oct 09:30:54:INFO:root:Difference: -2.4434405946142414e-5
19-Oct 09:30:54:INFO:root:Proportional increase: -12.243129328302281
19-Oct 09:30:54:INFO:root:---- Adjustment 3: 0.0018238077156721327
```

Trying to blow up the loop to reproduce the error. With `step_size=0.5` the convergence is even faster. `step_size=0.01` also converged.

Going in the wrong direction (`step_size=-0.1`) managed to blow it up:
```
19-Oct 09:44:51:INFO:root:Mean absolute wage gap: 1.3396099092278523e11
19-Oct 09:44:51:INFO:root:Gradient: 1.9055211199800068e10
19-Oct 09:44:51:INFO:root:Difference: 5.298608375392533e-5
19-Oct 09:44:51:INFO:root:Proportional increase: -2.432109017164174e-25
19-Oct 09:44:51:INFO:root:---- Adjustment 10: 1.9055211199800068e10
```

This is associated with a sector totally disappearing:
```
19-Oct 09:48:24:INFO:root:Wage gap:
2×3 Array{Float64,2}:
 0.673761  3.5846e11   0.967779
 0.582594  4.45306e11  0.9721
19-Oct 09:48:24:INFO:root:Labor shares:
2×3 Array{Float64,2}:
 0.5  5.0e-13  0.5
 0.5  5.0e-13  0.5
```

Back to normal step size. Low adjustment cost ensures constant wages:

```
19-Oct 09:51:12:INFO:root:Wage gap:
2×3 Array{Float64,2}:
 1.0004    1.00015  0.999654
 0.999067  1.00071  1.00013
19-Oct 09:51:12:INFO:root:Labor shares:
2×3 Array{Float64,2}:
 0.366299  0.148274  0.485427
 0.260195  0.254087  0.485718
19-Oct 09:51:12:INFO:root:Deviations:
2×3 Array{Float64,2}:
 5.55112e-17  2.77556e-17  5.55112e-17
 5.55112e-17  5.55112e-17  5.55112e-17
```

This cannot be tested with `S=1`. Expected labor will by definition be equal to actual labor.

There is a deviation from `L_njs_star` even with no adjustment costs. I set log productivity shocks with standard deviation of 0.2 and a low adjustment cost (`one_over_rho=2.0`) to enforce equal wages across sectors. Still, labor shares do not vary by more than 1.3 percentage point.
```
 2×3 Array{Float64,2}:
 0.883334  0.844817  1.00436
 1.16359   1.11103   1.54041
19-Oct 10:49:11:INFO:root:Wage gap:
2×3 Array{Float64,2}:
 1.00441  0.992284  0.998666
 1.0065   0.992997  0.999957
19-Oct 10:49:11:INFO:root:Labor shares:
2×3 Array{Float64,2}:
 0.380304  0.133281  0.486416
 0.268211  0.246088  0.4857
19-Oct 10:49:11:INFO:root:Deviations:
2×3 Array{Float64,2}:
 0.0119963  -0.0124865  0.000490195
 0.0132756  -0.0133582  8.2596e-5
```

It seems important to stay within the range of admissible labor shares and not wander towards the edges of the simplex.

# 2018-10-25

## Explore trajectory within simplex of labor shares

Two constraints on how labor shares can deviation from `L_njs_star`. They have to remain in the simples (sum to 1 and be positive). They also cannot deviation in expectation, because `L_njs_star` is the expected labor share. So maybe enforce `sum(deviation, 4)=0`?

Labor share seems to explode after a too-large positive step. Clearly we should not step outside the simplex!
```
3.627430347125549e-5
0.0010337950041388172

0.0010700693076100727
-0.0013675088049718562

9.997026489455864e-13
7.9343637929069875

0.6654293739274104
-1600.971241191903

4.999999999944995e-13
-5.825394528522916e7
```

With small-enough steps (`step_size=0.01`), there seems to be convergence even at the edge of the simplex.
```
3.627430347125549e-5
0.00010337950041388173

0.00013965380388513722
-3.4553603074329935e-5

0.00010510020081080728
-7.394455903240657e-6

9.770574490756662e-5
1.0119315041760835e-6

9.871767641174268e-5
1.1521386920515173e-6

9.986981510379422e-5
8.281200668252414e-7

0.00010069793517061947
5.879211552113975e-7

0.00010128585632583086
3.999139321885642e-7

0.00010168577025801942
2.8904775995975544e-7

0.00010197481801797918
1.9367923071517828e-7

0.00010216849724869436
1.8980039479827433e-7

0.00010235829764349263
1.3124200337770453e-7
```

Where do we enforce labor shares are in the simplex? In line 286 in `evaluate_utility`. This normalization does not respect direction of gradient. If, for example, gradient step asks for (0.9, 0.5, -0.4), this becomes (0.64, 0.36, 0.00), irrespective of where we started from. Better to ask: what is the biggest step we can make in the direction of the gradient while remaining strictly inside simplex.

`max_step_size` calculates maximum permissible step (for each n and s). This stops explosion, but labor share seems to oscillate around the edge of the simplex, hitting different boundaries at each step (6, 16, 6, 19, 16, 23, etc).
```
25-Oct 13:29:52:INFO:root:Where in simplex: (3.627430347125549e-5, (6,))
25-Oct 13:29:52:INFO:root:Maximum step size: [0.912477 0.0219675 0.559824 2.1986 3.78518 0.0320764 1.18061 0.503222 8.18811 2.00244 0.0339305 0.286422 0.306572 3.03256 4.70981 0.0730031 1.46702 0.146459 0.0602662 4.17549 0.0301096 0.205036 1.26943 1.65259 5.39229]
3.627430347125549e-5
0.0010337950041388172
25-Oct 13:29:52:INFO:root:---- Adjustment 1: 14.999579545220561
25-Oct 13:29:52:INFO:root:Where in simplex: (0.0004132746001339241, (16,))
25-Oct 13:29:52:INFO:root:Maximum step size: [1.13564 6.43288e-5 0.498125 2.60678 4.28432 2.87321e-5 1.71346 0.662054 11.1217 2.34357 0.000129213 0.38611 0.0784774 5.34376 6.44664 5.36728e-5 1.44645 0.242897 0.000657214 4.79343 4.66434e-6 0.475801 1.52652 2.11122 6.78365]
0.0010700693076100727
-0.001070069307623708
25-Oct 13:29:52:INFO:root:---- Adjustment 2: 34.18805444515999
25-Oct 13:29:52:INFO:root:Where in simplex: (9.999999999989864e-13, (6,))
25-Oct 13:29:52:INFO:root:Maximum step size: [1.43271 6.91924e-6 0.336059 3.10313 4.92082 4.15498e-5 2.89727 0.931369 14.3463 2.8237 8.89739e-6 0.669791 4.90698e-5 7.14535 8.11303 0.000302023 2.67731 0.0797327 2.94993e-5 5.69387 9.90778e-5 1.44791 2.01423 2.82139 8.57219]
9.999999999989864e-13
0.004007971803441376
25-Oct 13:29:53:INFO:root:---- Adjustment 3: 68.82348464904753
25-Oct 13:29:53:INFO:root:Where in simplex: (9.999999999954521e-13, (19,))
25-Oct 13:29:53:INFO:root:Maximum step size: [1.78991 7.85522e-5 2.48438 3.6911 5.92534 0.000563179 4.33828 1.32201 14.6831 3.53532 5.4916e-6 1.8334 9.3014e-5 9.99991 10.0688 1.40842e-5 2.49591 0.000268681 1.00802e-5 6.94456 2.01203e-6 11.5965 2.85202 4.74535 10.5718]
0.004007971804423149
-0.0003210212998963792
25-Oct 13:29:53:INFO:root:---- Adjustment 4: 109.98786787148647
25-Oct 13:29:53:INFO:root:Where in simplex: (9.999999999955717e-13, (16,))
25-Oct 13:29:53:INFO:root:Maximum step size: [2.4945 0.000459214 0.279526 4.45695 7.06039 7.64905e-6 7.05752 2.81443 18.5571 4.40844 1.30194e-6 6.80329 5.03052e-5 14.5207 12.3439 1.87535e-6 5.30868 0.000309475 7.77409e-6 8.39652 5.97379e-6 2.28573 4.40017 5.49427 12.386]
0.0036869505045104426
-0.00016963090201660032
25-Oct 13:29:53:INFO:root:---- Adjustment 5: 188.72744655494225
25-Oct 13:29:53:INFO:root:Where in simplex: (9.999999999956621e-13, (23,))
25-Oct 13:29:53:INFO:root:Maximum step size: [3.69338 5.43921e-6 0.739311 4.93555 9.12966 3.58504e-6 4.05642 3.17064 19.8829 5.61125 1.24486e-6 4.38812 0.000134543 16.6265 15.1173 4.43949e-7 7.79618 0.000205369 7.087e-7 11.5988 2.25969e-5 5.23239 20.1943 8.24677 15.0126]
0.0035173196024785845
-0.0005589923539492487
25-Oct 13:29:53:INFO:root:---- Adjustment 6: 101.95161434022357
25-Oct 13:29:53:INFO:root:Where in simplex: (9.999999999948747e-13, (7,))
25-Oct 13:29:53:INFO:root:Maximum step size: [5.07632 4.68541e-6 0.252393 5.96449 9.90934 2.09553e-6 4.37903 2.20496 20.7567 7.65102 1.44674e-5 2.85565 2.75242e-5 31.3843 19.0103 2.83032e-7 3.06008 0.000213629 3.48407e-5 11.0196 9.6257e-6 2.1224 2.60975 14.3702 17.0703]
0.002958327248514174
-0.00035498582148351613
25-Oct 13:29:54:INFO:root:---- Adjustment 7: 8812.935450379218
25-Oct 13:29:54:INFO:root:Where in simplex: (9.999999999861419e-13, (22,))
25-Oct 13:29:54:INFO:root:Maximum step size: [8.11066 0.000243481 2.84349 7.27191 11.4649 8.12365e-7 2.67319 0.713664 29.1681 8.72877 1.34046e-6 31.3033 2.00093e-5 37.0849 23.7209 8.04633e-6 9.00131 2.97692e-6 5.70308e-6 13.2912 3.46378e-9 11.011 2.82791 14.6179 19.4417]
```

Introducing additional dampening (`0.33*max_step_size`) seems to help.

# 2018-10-27

## Resolve remaining issues

Issue #18 closed.

No IO linkages. Recalculated A. Outer loop does not seem to converge.

Can outer loop leave simplex? No, because wage share is by defintion (0,1). Reducing step size helps. 

Prices and trad shares are ok in Period 1, but off in Period 2. 
```
27-Oct 18:00:14:INFO:root:Model trade shares: [0.984066, 0.681995, 0.994389]
27-Oct 18:00:14:INFO:root:Data trade shares: [0.970394, 0.938128, 0.99445]
27-Oct 18:00:14:INFO:root:US prices: [0.247416, 1.19084, 0.858142]
27-Oct 18:00:14:INFO:root:in the data: [1.21921, 1.17859, 1.17859]
```

Puzzle: some US prices go down in the model
```
27-Oct 18:05:50:INFO:root:US prices: [0.546038, 1.02751, 0.888753]
27-Oct 18:05:50:INFO:root:in the data: [1.21921, 1.17859, 1.17859]
27-Oct 18:05:50:INFO:root:Nominal world expenditure: 6.693994194508271e6
27-Oct 18:05:50:INFO:root:--------------In the data: [1.02443e7]
```
This may be because of deflation and expressing prices relative to US price index? Looks excessive for that, would require world inflation of 50% or so.

# 2018-1030
## Check output
There is a non-positive expenditure share in one of the scenarios:
```
29-Oct 11:31:57:INFO:root:--------------In the data: [2.31824e7]
ERROR: LoadError: On worker 3:
DomainError:
log will only return a complex result if called with a complex argument. Try log(complex(x)).
nan_dom_err at ./math.jl:300 [inlined]
log at ./math.jl:419 [inlined]
...
distance at /home/koren/projects/impvol-julia/experiments/trade_imbalance/kappa1972/../../../equilibrium.jl:49
middle_loop! at /home/koren/projects/impvol-julia/experiments/trade_imbalance/kappa1972/../../../equilibrium.jl:271
```

I added a lower bound for numerical zero, but then middle loop oscillates
```
30-Oct 17:38:27:INFO:root:------ Middle 44: 0.028858571230919443
30-Oct 17:38:27:DEBUG:root:------ BEGIN Inner loop
30-Oct 17:38:27:DEBUG:root:-------- Inner 1: 0.006595752733840549
30-Oct 17:38:27:DEBUG:root:-------- Inner 2: 0.00534547637915506
30-Oct 17:38:27:DEBUG:root:-------- Inner 3: 0.0034377046310015936
30-Oct 17:38:28:DEBUG:root:-------- Inner 4: 0.002160947679340457
30-Oct 17:38:28:DEBUG:root:-------- Inner 5: 0.0013469907958464184
30-Oct 17:38:28:DEBUG:root:-------- Inner 6: 0.0008475161457900738
30-Oct 17:38:28:DEBUG:root:------ END Inner loop
30-Oct 17:38:28:INFO:root:------ Middle 45: 0.7411762538111589
30-Oct 17:38:28:DEBUG:root:------ BEGIN Inner loop
30-Oct 17:38:28:DEBUG:root:-------- Inner 1: 0.0005468679924894514
30-Oct 17:38:28:DEBUG:root:------ END Inner loop
30-Oct 17:38:28:INFO:root:------ Middle 46: 0.7087367860040116
30-Oct 17:38:28:DEBUG:root:------ BEGIN Inner loop
30-Oct 17:38:28:DEBUG:root:-------- Inner 1: 0.003475341346266561
30-Oct 17:38:29:DEBUG:root:-------- Inner 2: 0.0055330339781147425
30-Oct 17:38:29:DEBUG:root:-------- Inner 3: 0.004031365631101054
30-Oct 17:38:29:DEBUG:root:-------- Inner 4: 0.0028935406483291707
30-Oct 17:38:29:DEBUG:root:-------- Inner 5: 0.0021044124450564734
30-Oct 17:38:29:DEBUG:root:-------- Inner 6: 0.001535315468523144
30-Oct 17:38:29:DEBUG:root:-------- Inner 7: 0.0011230838106543845
30-Oct 17:38:30:DEBUG:root:-------- Inner 8: 0.0008229403237599034
30-Oct 17:38:30:DEBUG:root:------ END Inner loop
30-Oct 17:38:30:INFO:root:------ Middle 47: 0.17211336901095803
30-Oct 17:38:30:DEBUG:root:------ BEGIN Inner loop
30-Oct 17:38:30:DEBUG:root:-------- Inner 1: 0.0006039731566345187
30-Oct 17:38:30:DEBUG:root:------ END Inner loop
30-Oct 17:38:30:INFO:root:------ Middle 48: 0.028858410790122604
```

Middle loop is looking for expenditure shares, so should also remain in simplex. Both components are inside simplex, maybe reduce step size?

Around numerical zero, a log distance measure is very sensitive. For shares, I would switch to percentage point distance.

Inner loop gets stuck
```
30-Oct 18:44:12:DEBUG:root:-------- Inner 113: 0.0029659561893971872
30-Oct 18:44:12:DEBUG:root:-------- Inner 114: 0.002395913356269078
30-Oct 18:44:12:DEBUG:root:-------- Inner 115: 0.0029574539028234077
30-Oct 18:44:12:DEBUG:root:-------- Inner 116: 0.002382642240396397
30-Oct 18:44:12:DEBUG:root:-------- Inner 117: 0.0029511202077374702
30-Oct 18:44:13:DEBUG:root:-------- Inner 118: 0.0023721295794866238
30-Oct 18:44:13:DEBUG:root:-------- Inner 119: 0.0029465367901150377
30-Oct 18:44:13:DEBUG:root:-------- Inner 120: 0.0023638847092739496
30-Oct 18:44:13:DEBUG:root:-------- Inner 121: 0.002943349924611318
30-Oct 18:44:13:DEBUG:root:-------- Inner 122: 0.0023574894971424335
30-Oct 18:44:13:DEBUG:root:-------- Inner 123: 0.0029412650265131994
30-Oct 18:44:14:DEBUG:root:-------- Inner 124: 0.0023525918811992917
30-Oct 18:44:14:DEBUG:root:-------- Inner 125: 0.0029400399887136866
30-Oct 18:44:14:DEBUG:root:-------- Inner 126: 0.0023488985536536918
30-Oct 18:44:14:DEBUG:root:-------- Inner 127: 0.002939477364985182
30-Oct 18:44:14:DEBUG:root:-------- Inner 128: 0.0023461670243482267
30-Oct 18:44:15:DEBUG:root:-------- Inner 129: 0.002939417249673482
30-Oct 18:44:15:DEBUG:root:-------- Inner 130: 0.0023441981457784913
30-Oct 18:44:15:DEBUG:root:-------- Inner 131: 0.002939731173711331
30-Oct 18:44:15:DEBUG:root:-------- Inner 132: 0.0023428294201043474
```

Log distance is not a good measure for prices, either, use difference in predicted country shares. With price to (-theta), middle loop still oscillates
```
30-Oct 20:41:19:INFO:root:------ Middle 96: 0.0014719370573773188
30-Oct 20:41:19:INFO:root:------ Middle 97: 0.0008296230463568038
30-Oct 20:41:20:INFO:root:------ Middle 98: 0.0014006377369308487
30-Oct 20:41:20:INFO:root:------ Middle 99: 0.0007991132846278228
30-Oct 20:41:21:INFO:root:------ Middle 100: 0.0014624854880997877
```
# 2018-10-31

## Debug middle loop for trade imbalance case

`experiments/baseline/actual` running fine. In `trade_imbalance`, distance reduces up to a point, at which oscillation starts.
```
31-Oct 10:42:33:INFO:root:------ Middle 1: 0.020087978317916387
31-Oct 10:42:36:INFO:root:------ Middle 2: 0.011652563694581347
31-Oct 10:42:36:INFO:root:------ Middle 3: 0.006437027036450518
31-Oct 10:42:40:INFO:root:------ Middle 4: 0.004150244550192232
31-Oct 10:42:40:INFO:root:------ Middle 5: 0.0023019164014744034
31-Oct 10:42:42:INFO:root:------ Middle 6: 0.002043543946662378
31-Oct 10:42:42:INFO:root:------ Middle 7: 0.0011459547865934937
31-Oct 10:42:44:INFO:root:------ Middle 8: 0.0014189201742071077
31-Oct 10:42:44:INFO:root:------ Middle 9: 0.0007939137525280888
31-Oct 10:42:45:INFO:root:------ Middle 10: 0.0014162449268553495
31-Oct 10:42:45:INFO:root:------ Middle 11: 0.0008199946488115949
31-Oct 10:42:46:INFO:root:------ Middle 12: 0.0013745745781782964
31-Oct 10:42:47:INFO:root:------ Middle 13: 0.0007632439102759337
31-Oct 10:42:48:INFO:root:------ Middle 14: 0.001207056069492887
31-Oct 10:42:48:INFO:root:------ Middle 15: 0.000696214758000288
```
Reduce inner tolerance to 0.0001 to see if it affects where oscillation occurs. The answer is no.
```
31-Oct 10:46:10:INFO:root:------ Middle 1: 0.020087978317916387
31-Oct 10:46:15:INFO:root:------ Middle 2: 0.01181989324586644
31-Oct 10:46:15:INFO:root:------ Middle 3: 0.006516140887436318
31-Oct 10:46:20:INFO:root:------ Middle 4: 0.004190241006031715
31-Oct 10:46:21:INFO:root:------ Middle 5: 0.0023140470379854855
31-Oct 10:46:24:INFO:root:------ Middle 6: 0.002124708958196941
31-Oct 10:46:24:INFO:root:------ Middle 7: 0.0011761832113734727
31-Oct 10:46:26:INFO:root:------ Middle 8: 0.0014420729253306187
31-Oct 10:46:26:INFO:root:------ Middle 9: 0.0007999099736026548
31-Oct 10:46:28:INFO:root:------ Middle 10: 0.0014572418746174375
31-Oct 10:46:28:INFO:root:------ Middle 11: 0.0008216580045495448
31-Oct 10:46:29:INFO:root:------ Middle 12: 0.0014267370777544227
31-Oct 10:46:30:INFO:root:------ Middle 13: 0.0007886929002944327
31-Oct 10:46:31:INFO:root:------ Middle 14: 0.0013339401942805069
31-Oct 10:46:31:INFO:root:------ Middle 15: 0.0007497840098704164
```

CES works fine.

After deflating Snt, we again have convergence.

# 2018-11-01

## Compare model to data

Model GDP exhibits huge volatility. Is it some data artifact? No, real GDP summed in simulation results fluctuates wildly. For USA, for example,
```
36-element Array{Float64,1}:
 1.1344e6 
 8.55546e5
 7.81789e5
 2.28051e6
 ...
```

In the data, productivity growth is normal, +- 4-5% each year. However, from period 2 onwards, productivity in State 1 is not = productivity in the data.
```
[results[2][2][:A_njs][1,end,:,1] parameters[:A][1,end,:,2]]
24×2 Array{Float64,2}:
 11167.3     7131.23 
   468.8      549.146
  4153.15    4206.29 
  2457.03    2347.47 
  2678.18    2773.51 
  1786.46    1787.4  
  1550.43    1460.16 
  1435.31    1502.05 
  6761.0     9071.16 
   104.712    103.896
  1565.77    1932.88 
  3420.29    3282.29 
  5284.93    5389.49 
  1859.35    1726.08 
  3422.0     3449.72 
  3899.27    3651.91 
  1549.28    1525.91 
  4631.08    4155.96 
  3161.06    2991.72 
  7984.78    7968.36 
  1709.28    1468.59 
  3664.7     3475.8  
  3824.77    3901.97 
 41172.0    62202.5  
```
Productivities also match for period 36.
```
results[36][2][:A_njs][1,end,:,1] ≈ parameters[:A][1,end,:,36]
true
```
Can this be related to detrending? The match still holds for parameters, added unit tests for this. However, results seem to have a different Anjs.

In test run (`experiments/baseline/actual/test.jl`), US GDP seems to grow smoothly,
```
3-element Array{Float64,1}:
 1.1344e6 
 1.20329e6
 1.18642e6
```

This is inconsistent with results file from server run:
```
3-element Array{Float64,1}:
 1.1344e6 
 8.55546e5
 7.81789e5
```

Local run has good GDP numbers, see notebook.

# 2018-11-05
## Debug labor adjustment scenario
All scenarios not involving labor adjustment have regular volatilities [3451744ebea18560491614ca94e013cb5c53402c]. Those with labor adjustment have 10x, 30x fluctuationa. Likely wage calibration is at fault. Indeed, there are huge fluctuations in calibrate wages
```
 1.12031e6  1.1735e6   6.79931e5     1.13751e7  1.23575e7  1.2997e7
 1.12031e6  6.26857e5  6.79931e5     5.85601e6  2.09113e7  1.2997e7
```
Recovering the labor shares consistent with these wage ratios, some of them are greater than one. Last two sectors for US,
```
 0.00765139  0.00758658  0.0130356   0.0122568      0.0051943   0.0048989  
 0.71855     1.41978     1.44574     0.239047       0.504511    0.850966   
```

While wage ratio may be approximately 1 across all periods and countries, this is not true per period per country. For example, USA:
```
mean(wage_ratio[1,end,:,:], 1)
1×36 Array{Float64,2}:
 1.0  1.12633  0.901342  0.887479  …  1.02272  0.957281  0.971533  1.0
```

Expressing labor shares from trend deviations gives smooth fluctuations. See algorithm.tex.

# 2018-11-07
## Evaluate S=1000
Estimate volatilities differ by 1% on average. Maybe set S=100 by default and test agains S=10000?

# 2018-12-03
## New CES results
Removing sectoral shocks results in even higher volatility.
```
South Korea 0.0030050596714052257 0.002944000329299149  7.611564497864669 7.712889346155691 2.074026334114317 -3441.7403857812665 3443.8144121153773
Spain 0.00047496319368355895  0.0012767659599013084 0.024062947713052854  0.024488399211282888  -62.79950996498429  -33.32259095182335  -29.47691901316096
```