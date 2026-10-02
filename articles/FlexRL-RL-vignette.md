# FlexRL-RL-vignette

``` r

library(FlexRL)
#> If you are happy with FlexRL, please cite us!  Also, if you are unhappy, please cite us anyway.
```

`FlexRL` is a package proposing flexible probabilistic **R**ecord
**L**inkage with diagnostic tools for downstream inference on linked
data. The record linkage algorithm [`StEM()`](../reference/stEM.md) uses
a **St**ochastic **E**xpectation **M**aximisation approach to combine 2
data sources and outputs the set of records referring to the same
entities.

More details on the [record linkage
methodology](https://doi.org/10.1093/jrsssc/qlaf016); on the [false
discovery proportion estimation](https://doi.org/10.1002/sim.70292); and
on the potential risks of downstream inference on linked data (Robach et
al., in preparation).

This vignette uses example subsets from the
[SHIW](https://github.com/robachowyk/FlexRL-experiments/tree/main/SHIWApplication/SHIWData)
data, from the
[NLTCS](https://www.icpsr.umich.edu/web/NACDA/studies/9681/versions/V5)
data, and synthetic data to showcase how `FlexRL` works, how
[`StEM()`](../reference/stEM.md) adapts to dynamic PIVs, how
`RL_diagnostics` may guide downstream inference.

We use 4 **P**artially **I**dentifying **V**ariable**s** to link the
records: 2 variables such that `dynamics = "stable"` and 2 variables
that may change over time. Among these 2 dynamic PIVs, we will model the
changes of one (`dynamics = "structured"`) and we will let the algorithm
flexible about the other (`dynamics = "flexible"`).

``` r

PIVs_config <- list( V1 = list(dynamics = "stable",
                               bound_mistakes = c(0.10,0.10),
                               fix_mistakes = c(NA,NA)),
                    V2 = list(dynamics = "stable",
                              bound_mistakes = c(0.10,0.10),
                              fix_mistakes = c(NA,NA)),
                    V3 = list(dynamics = "flexible",
                              bound_mistakes = c(NA,NA),
                              fix_mistakes = c(NA,NA)),
                    V4 = list(dynamics = "structured",
                              bound_mistakes = c(NA,NA),
                              fix_mistakes = c(0.02,0.02),
                              cond_hazard_cov = list(cov1 = c("Xe", "Xf"),
                                                     cov2 = c()))
)
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$dynamics != "structured")
PIVs_type = list( V1 = FALSE, V2 = FALSE, V3 = FALSE, V4 = TRUE )
```

With PIVs that evolve over time, when modelling is possible, it is
natural to parametrise the probability that underlying true values of
the registered data for a pair of linked records coincide using a
survival function. It is possible to use an exponential model, whose
baseline hazard is constant in time; a Weibull model, whose baseline
hazard is monotone in time; a Gompertz model, whose baseline hazard
grows or decays exponentially in time; a piecewise-constant model, whose
baseline hazard takes one value on each interval defined by given cut
points; a custom model, supplied as a survival function S(X, alpha,
times) together with a number of parameters and their starting values.
These models may assume proportional hazards, accelerated failure time,
proportional odds, additive hazards, … If the covariates are
categorical, make sure to encode them into dummies and drop one category
(an intercept is included automatically). Since only the agreement of
the true values at the observed time gap is used, the likelihood is the
same for every model, -sum(Hequal \* log(S) + (1 - Hequal) \* log(1 -
S)), and the parameters are estimated with stats::nlminb() at each
M-step.

``` r

n_values  <- c(8, 9, 10, 15)
p_mistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
                  V3 = c(0.02, 0.02), V4 = c(0.02, 0.02))
p_missing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
                  V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
cond_hazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = log(c(0.7, 0.6, 0.5)))

gen_data <- simulate_data(
  PIVs_config = PIVs_config, n_values = n_values, n_records = c(450, 500), 
  n_links = 400, p_mistake = p_mistake, p_missing = p_missing, 
  cond_hazard_params = cond_hazard_params, enforce_estimability = TRUE,
  model_dynamics = survival_model("exponential")
)
```

We pre-process the data for the record linkage task: remove the records
for which the linking variables are outside of the common support, check
that PIVs are not strongly associated with each other (which may degrade
record linkage performance), map PIV values to natural numbers and
encode missing values as 0. The StEM algorithm denotes the smallest
source as `A` and the largest as `B`.

``` r

prep_data <- prepare_data(
  data1 = gen_data$data1, data2 = gen_data$data2, label1 = "1", label2 = "2", 
  PIVs_config = PIVs_config, same_mistakes = TRUE, uniq_id = "entity_id",
  restrict_support_intersection = TRUE
)
#> '2' is the larger source, saved as B; '1' saved as A.
```

``` r

head(prep_data$encodedA)
#>   V1 V2 V3 V4 change         date        Xe        Xf local_id entity_id source
#> 1  8  3  4  9  FALSE 9.358000e-03 0.8473464  1.826180        1         1      1
#> 2  3  6  3 12  FALSE 5.267033e-03 0.9779230  1.292910        2         2      1
#> 3  5  6  9 12  FALSE 9.798704e-03 1.2878413 -1.141012        3         3      1
#> 4  8  8  1 12  FALSE 3.133263e-03 2.8660074  1.338841        4         4      1
#> 5  0  8  3  5  FALSE 8.280712e-03 0.7181951  2.373188        5         5      1
#> 6  6  8  3  8  FALSE 2.433672e-05 1.2339214  1.105553        6         6      1
```

``` r

head(prep_data$encodedB)
#>   V1 V2 V3 V4 change         date local_id entity_id source
#> 1  7  6  4  9  FALSE 0.0033904564        1         1      2
#> 2  3  6  3 12  FALSE 0.0017498800        2         2      2
#> 3  5  6  9 12  FALSE 0.0017083757        3         3      2
#> 4  8  8  1 12  FALSE 0.0073588952        4         4      2
#> 5  8  8  3  5  FALSE 0.0090912612        5         5      2
#> 6  6  8  3  8  FALSE 0.0005722134        6         6      2
```

For this example we know the true linkage structure.

``` r

prep_data$true_pairs
#>       1   2
#> 1     1   1
#> 2     2   2
#> 3     3   3
#> 4     4   4
#> 5     5   5
#> 6     6   6
#> 7     7   7
#> 8     8   8
#> 9     9   9
#> 10   10  10
#> 11   11  11
#> 12   12  12
#> 13   13  13
#> 14   14  14
#> 15   15  15
#> 16   16  16
#> 17   17  17
#> 18   18  18
#> 19   19  19
#> 20   20  20
#> 21   21  21
#> 22   22  22
#> 23   23  23
#> 24   24  24
#> 25   25  25
#> 26   26  26
#> 27   27  27
#> 28   28  28
#> 29   29  29
#> 30   30  30
#> 31   31  31
#> 32   32  32
#> 33   33  33
#> 34   34  34
#> 35   35  35
#> 36   36  36
#> 37   37  37
#> 38   38  38
#> 39   39  39
#> 40   40  40
#> 41   41  41
#> 42   42  42
#> 43   43  43
#> 44   44  44
#> 45   45  45
#> 46   46  46
#> 47   47  47
#> 48   48  48
#> 49   49  49
#> 50   50  50
#> 51   51  51
#> 52   52  52
#> 53   53  53
#> 54   54  54
#> 55   55  55
#> 56   56  56
#> 57   57  57
#> 58   58  58
#> 59   59  59
#> 60   60  60
#> 61   61  61
#> 62   62  62
#> 63   63  63
#> 64   64  64
#> 65   65  65
#> 66   66  66
#> 67   67  67
#> 68   68  68
#> 69   69  69
#> 70   70  70
#> 71   71  71
#> 72   72  72
#> 73   73  73
#> 74   74  74
#> 75   75  75
#> 76   76  76
#> 77   77  77
#> 78   78  78
#> 79   79  79
#> 80   80  80
#> 81   81  81
#> 82   82  82
#> 83   83  83
#> 84   84  84
#> 85   85  85
#> 86   86  86
#> 87   87  87
#> 88   88  88
#> 89   89  89
#> 90   90  90
#> 91   91  91
#> 92   92  92
#> 93   93  93
#> 94   94  94
#> 95   95  95
#> 96   96  96
#> 97   97  97
#> 98   98  98
#> 99   99  99
#> 100 100 100
#> 101 101 101
#> 102 102 102
#> 103 103 103
#> 104 104 104
#> 105 105 105
#> 106 106 106
#> 107 107 107
#> 108 108 108
#> 109 109 109
#> 110 110 110
#> 111 111 111
#> 112 112 112
#> 113 113 113
#> 114 114 114
#> 115 115 115
#> 116 116 116
#> 117 117 117
#> 118 118 118
#> 119 119 119
#> 120 120 120
#> 121 121 121
#> 122 122 122
#> 123 123 123
#> 124 124 124
#> 125 125 125
#> 126 126 126
#> 127 127 127
#> 128 128 128
#> 129 129 129
#> 130 130 130
#> 131 131 131
#> 132 132 132
#> 133 133 133
#> 134 134 134
#> 135 135 135
#> 136 136 136
#> 137 137 137
#> 138 138 138
#> 139 139 139
#> 140 140 140
#> 141 141 141
#> 142 142 142
#> 143 143 143
#> 144 144 144
#> 145 145 145
#> 146 146 146
#> 147 147 147
#> 148 148 148
#> 149 149 149
#> 150 150 150
#> 151 151 151
#> 152 152 152
#> 153 153 153
#> 154 154 154
#> 155 155 155
#> 156 156 156
#> 157 157 157
#> 158 158 158
#> 159 159 159
#> 160 160 160
#> 161 161 161
#> 162 162 162
#> 163 163 163
#> 164 164 164
#> 165 165 165
#> 166 166 166
#> 167 167 167
#> 168 168 168
#> 169 169 169
#> 170 170 170
#> 171 171 171
#> 172 172 172
#> 173 173 173
#> 174 174 174
#> 175 175 175
#> 176 176 176
#> 177 177 177
#> 178 178 178
#> 179 179 179
#> 180 180 180
#> 181 181 181
#> 182 182 182
#> 183 183 183
#> 184 184 184
#> 185 185 185
#> 186 186 186
#> 187 187 187
#> 188 188 188
#> 189 189 189
#> 190 190 190
#> 191 191 191
#> 192 192 192
#> 193 193 193
#> 194 194 194
#> 195 195 195
#> 196 196 196
#> 197 197 197
#> 198 198 198
#> 199 199 199
#> 200 200 200
#> 201 201 201
#> 202 202 202
#> 203 203 203
#> 204 204 204
#> 205 205 205
#> 206 206 206
#> 207 207 207
#> 208 208 208
#> 209 209 209
#> 210 210 210
#> 211 211 211
#> 212 212 212
#> 213 213 213
#> 214 214 214
#> 215 215 215
#> 216 216 216
#> 217 217 217
#> 218 218 218
#> 219 219 219
#> 220 220 220
#> 221 221 221
#> 222 222 222
#> 223 223 223
#> 224 224 224
#> 225 225 225
#> 226 226 226
#> 227 227 227
#> 228 228 228
#> 229 229 229
#> 230 230 230
#> 231 231 231
#> 232 232 232
#> 233 233 233
#> 234 234 234
#> 235 235 235
#> 236 236 236
#> 237 237 237
#> 238 238 238
#> 239 239 239
#> 240 240 240
#> 241 241 241
#> 242 242 242
#> 243 243 243
#> 244 244 244
#> 245 245 245
#> 246 246 246
#> 247 247 247
#> 248 248 248
#> 249 249 249
#> 250 250 250
#> 251 251 251
#> 252 252 252
#> 253 253 253
#> 254 254 254
#> 255 255 255
#> 256 256 256
#> 257 257 257
#> 258 258 258
#> 259 259 259
#> 260 260 260
#> 261 261 261
#> 262 262 262
#> 263 263 263
#> 264 264 264
#> 265 265 265
#> 266 266 266
#> 267 267 267
#> 268 268 268
#> 269 269 269
#> 270 270 270
#> 271 271 271
#> 272 272 272
#> 273 273 273
#> 274 274 274
#> 275 275 275
#> 276 276 276
#> 277 277 277
#> 278 278 278
#> 279 279 279
#> 280 280 280
#> 281 281 281
#> 282 282 282
#> 283 283 283
#> 284 284 284
#> 285 285 285
#> 286 286 286
#> 287 287 287
#> 288 288 288
#> 289 289 289
#> 290 290 290
#> 291 291 291
#> 292 292 292
#> 293 293 293
#> 294 294 294
#> 295 295 295
#> 296 296 296
#> 297 297 297
#> 298 298 298
#> 299 299 299
#> 300 300 300
#> 301 301 301
#> 302 302 302
#> 303 303 303
#> 304 304 304
#> 305 305 305
#> 306 306 306
#> 307 307 307
#> 308 308 308
#> 309 309 309
#> 310 310 310
#> 311 311 311
#> 312 312 312
#> 313 313 313
#> 314 314 314
#> 315 315 315
#> 316 316 316
#> 317 317 317
#> 318 318 318
#> 319 319 319
#> 320 320 320
#> 321 321 321
#> 322 322 322
#> 323 323 323
#> 324 324 324
#> 325 325 325
#> 326 326 326
#> 327 327 327
#> 328 328 328
#> 329 329 329
#> 330 330 330
#> 331 331 331
#> 332 332 332
#> 333 333 333
#> 334 334 334
#> 335 335 335
#> 336 336 336
#> 337 337 337
#> 338 338 338
#> 339 339 339
#> 340 340 340
#> 341 341 341
#> 342 342 342
#> 343 343 343
#> 344 344 344
#> 345 345 345
#> 346 346 346
#> 347 347 347
#> 348 348 348
#> 349 349 349
#> 350 350 350
#> 351 351 351
#> 352 352 352
#> 353 353 353
#> 354 354 354
#> 355 355 355
#> 356 356 356
#> 357 357 357
#> 358 358 358
#> 359 359 359
#> 360 360 360
#> 361 361 361
#> 362 362 362
#> 363 363 363
#> 364 364 364
#> 365 365 365
#> 366 366 366
#> 367 367 367
#> 368 368 368
#> 369 369 369
#> 370 370 370
#> 371 371 371
#> 372 372 372
#> 373 373 373
#> 374 374 374
#> 375 375 375
#> 376 376 376
#> 377 377 377
#> 378 378 378
#> 379 379 379
#> 380 380 380
#> 381 381 381
#> 382 382 382
#> 383 383 383
#> 384 384 384
#> 385 385 385
#> 386 386 386
#> 387 387 387
#> 388 388 388
#> 389 389 389
#> 390 390 390
#> 391 391 391
#> 392 392 392
#> 393 393 393
#> 394 394 394
#> 395 395 395
#> 396 396 396
#> 397 397 397
#> 398 398 398
#> 399 399 399
#> 400 400 400
true_pairs <- do.call(paste, c(prep_data$true_pairs, list(sep="_")))
```

The number of unique values per **PIV** gives information on their
discriminating power:

``` r

prep_data$n_values
#> V1 V2 V3 V4 
#>  8  9 10 15
```

We gauge how hard the record linkage task is using a naive linkage
approach: matching records on exact agreement of identifying
information.

``` r

naive_fit <- naive_linkage(PIVs = PIVs, 
                           encodedA = prep_data$encodedA, 
                           encodedB = prep_data$encodedB)

df_results = data.frame( matrix(NA, nrow = 11, ncol = 0) )
rownames(df_results) = c("TP              ",
                         "FP              ",
                         "FN              ",
                         "sensitivity     ",
                         "FDP             ",
                         "hat FDP score   ",
                         "hat FDP score*  ",
                         "hat FDP synth*  ",
                         "max SMD         ",
                         "min support IoU ",
                         "max MMD         ")

linked_pairs    = do.call(paste, c(naive_fit, list(sep = "_")))
true_positive   = length( intersect(linked_pairs, true_pairs) ) 
false_positive  = length( setdiff(linked_pairs, true_pairs) ) 
false_negative  = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 
df_results[1:5,"Naive     "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

round(df_results, 2)
#>                  Naive     
#> TP                   258.00
#> FP                    94.00
#> FN                   142.00
#> sensitivity            0.64
#> FDP                    0.27
#> hat FDP score            NA
#> hat FDP score*           NA
#> hat FDP synth*           NA
#> max SMD                  NA
#> min support IoU          NA
#> max MMD                  NA
```

We run record linkage with `FlexRL`. The `data` should contain encodedA,
encodedB, n_values, PIVs_config, same_mistakes. `StEM_iter`,
`StEM_burnin` and `gibbs_iter`, `gibbs_burnin` are total number of
iterations (including burn-in) and iterations to be discarded as
burn-in. If `music_on` the algorithm opens a short tune in the browser
when it finishes. You can `save_info_iter` (save environment information
at each iteration) to `new_directory`. You can set starting values of
the algorithm for parameters gamma and phi with `gamma0`, `phiA0`,
`phiB0`, and change the number of samples drawn a posteriori to provide
a linkage estimate with `n_post_samp` (set to 1000 as default).

``` r

fit <- StEM(data = prep_data, StEM_iter = 30, StEM_burnin = 15,
           gibbs_iter = 20, gibbs_burnin = 10, music_on = FALSE)
#> FlexRL
#> Running StEM algorithm ■■■                                7% | iter 2/30 [1.4s]
#> Running StEM algorithm ■■■■                              10% | iter 3/30 [1.9s]
#> Running StEM algorithm ■■■■■■■■■■                        30% | iter 9/30 [4.8s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■                 53% | iter 16/30 [7.9s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■          77% | iter 23/30 [11s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% | iter 29/30 [13.8…
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 30/30 [14.2…
#> 
#> Drawing Delta          ■■                                 2% [2.4s]
#> Drawing Delta          ■■■                                5% [5.4s]
#> Drawing Delta          ■■■                                8% [8.4s]
#> Drawing Delta          ■■■■                              11% [11.4s]
#> Drawing Delta          ■■■■■                             14% [14.4s]
#> Drawing Delta          ■■■■■■                            16% [17.3s]
#> Drawing Delta          ■■■■■■■                           19% [20.4s]
#> Drawing Delta          ■■■■■■■■                          22% [23.4s]
#> Drawing Delta          ■■■■■■■■                          25% [26.4s]
#> Drawing Delta          ■■■■■■■■■                         28% [29.4s]
#> Drawing Delta          ■■■■■■■■■■                        30% [32.4s]
#> Drawing Delta          ■■■■■■■■■■■                       33% [35.4s]
#> Drawing Delta          ■■■■■■■■■■■■                      36% [38.4s]
#> Drawing Delta          ■■■■■■■■■■■■■                     39% [41.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    42% [44.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    45% [47.6s]
#> Drawing Delta          ■■■■■■■■■■■■■■■                   48% [50.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  50% [53.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■                 53% [56.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                56% [59.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               59% [1m 2.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              62% [1m 5.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              65% [1m 8.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■             68% [1m 11.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            70% [1m 14.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■           73% [1m 17.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■          76% [1m 20.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         79% [1m 23.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% [1m 26.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% [1m 29.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% [1m 32.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% [1m 35.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% [1m 38.6s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% [1m 41.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% [1m 44.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% [1m 45.5s]
#> 
```

The missing values and potential mistakes in the registration, the size
of the overlapping set of records between the files, the total number of
entities in the population, the low discriminative power of PIVs, PIVs
potential dynamics across time, PIVs distribution, dependencies among
PIVs, are obstacles to the linkage.

The algorithm returns: `Delta`, `gamma`, `eta`, `alpha`, `phi`.

``` r

delta_result = fit$Delta
delta_result = delta_result[delta_result$x>0.5, ]

linked_pairs    = do.call(paste, c(delta_result[,c("i","j")], list(sep = "_")))
true_positive   = length( intersect(linked_pairs, true_pairs) ) 
false_positive  = length( setdiff(linked_pairs, true_pairs) ) 
false_negative  = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive)  

df_results[1:5,"FlexRL    "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

round(df_results, 2)
#>                  Naive      FlexRL    
#> TP                   258.00     270.00
#> FP                    94.00      46.00
#> FN                   142.00     130.00
#> sensitivity            0.64       0.68
#> FDP                    0.27       0.15
#> hat FDP score            NA         NA
#> hat FDP score*           NA         NA
#> hat FDP synth*           NA         NA
#> max SMD                  NA         NA
#> min support IoU          NA         NA
#> max MMD                  NA         NA
```

Let us run record linkage diagnostics:

``` r

diag <- RL_diagnostics(fit = fit, encodedA = prep_data$encodedA, encodedB = prep_data$encodedB,
                      compare_vars = PIVs, vars_type_cont = PIVs_type, true_pairs = prep_data$true_pairs, 
                      FDP_estimation = TRUE, RL_method = "FlexRL", data = prep_data,
                      maxIter4CV = 1, n_repeats = 2)
#> Iteration: 0, Accuracy: 47.38%
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> Running StEM algorithm ■■■■                              10% | iter 3/30 [1.3s]
#> Running StEM algorithm ■■■■■                             13% | iter 4/30 [1.7s]
#> Running StEM algorithm ■■■■■■■■■■■■                      37% | iter 11/30 [5s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■                57% | iter 17/30 [7.7s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■         80% | iter 24/30 [10.9…
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 30/30 [13.5…
#> 
#> Drawing Delta          ■                                  1% [1.1s]
#> Drawing Delta          ■■                                 3% [3s]
#> Drawing Delta          ■■■                                6% [6.1s]
#> Drawing Delta          ■■■■                               8% [9.1s]
#> Drawing Delta          ■■■■                              11% [12.1s]
#> Drawing Delta          ■■■■■                             14% [15.1s]
#> Drawing Delta          ■■■■■■                            17% [18.1s]
#> Drawing Delta          ■■■■■■■                           20% [21s]
#> Drawing Delta          ■■■■■■■■                          23% [24s]
#> Drawing Delta          ■■■■■■■■■                         26% [27s]
#> Drawing Delta          ■■■■■■■■■■                        29% [30s]
#> Drawing Delta          ■■■■■■■■■■                        32% [33s]
#> Drawing Delta          ■■■■■■■■■■■                       34% [36s]
#> Drawing Delta          ■■■■■■■■■■■■                      37% [39.1s]
#> Drawing Delta          ■■■■■■■■■■■■■                     40% [42.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    43% [45.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■                   46% [48.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  49% [51.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  52% [54.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■                 54% [57s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                57% [1m]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               60% [1m 3.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              63% [1m 6.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■             66% [1m 9.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            69% [1m 12.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            72% [1m 15s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■           74% [1m 18s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■          77% [1m 21.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         80% [1m 24.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% [1m 27.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% [1m 30s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% [1m 33s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      92% [1m 36s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% [1m 39s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% [1m 42.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% [1m 44.9s]
#> 
#> Iteration: 0, Accuracy: 48.14%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> FlexRL
#> Running StEM algorithm ■■■■                              10% | iter 3/30 [1.4s]
#> Running StEM algorithm ■■■■■■■■                          23% | iter 7/30 [3.2s]
#> Running StEM algorithm ■■■■■■■■■■■■■■                    43% | iter 13/30 [5.9s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■             67% | iter 20/30 [9s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% | iter 27/30 [12.3…
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 30/30 [13.6…
#> 
#> Drawing Delta          ■                                  1% [1.4s]
#> Drawing Delta          ■■                                 4% [4.4s]
#> Drawing Delta          ■■■                                7% [7.4s]
#> Drawing Delta          ■■■■                              10% [10.3s]
#> Drawing Delta          ■■■■■                             13% [13.4s]
#> Drawing Delta          ■■■■■■                            16% [16.3s]
#> Drawing Delta          ■■■■■■■                           18% [19.4s]
#> Drawing Delta          ■■■■■■■                           21% [22.3s]
#> Drawing Delta          ■■■■■■■■                          24% [25.4s]
#> Drawing Delta          ■■■■■■■■■                         27% [28.4s]
#> Drawing Delta          ■■■■■■■■■■                        30% [31.3s]
#> Drawing Delta          ■■■■■■■■■■■                       33% [34.3s]
#> Drawing Delta          ■■■■■■■■■■■■                      36% [37.4s]
#> Drawing Delta          ■■■■■■■■■■■■■                     38% [40.4s]
#> Drawing Delta          ■■■■■■■■■■■■■                     41% [43.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    44% [46.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■                   47% [49.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  50% [52.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■                 53% [55.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                56% [58.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               58% [1m 1.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               61% [1m 4.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              64% [1m 7.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■             67% [1m 10.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            70% [1m 13.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■           73% [1m 16.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■          76% [1m 19.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         79% [1m 22.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         81% [1m 25.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% [1m 28.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% [1m 31.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% [1m 34.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% [1m 37.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% [1m 40.3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% [1m 43.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% [1m 44.7s]
#> 
#> FlexRL results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.13   0.12   0.12   0.11   0.11   0.10   0.10
#> FDP synth data estimator    0.24   0.22   0.21   0.17   0.14   0.14   0.14
#> Linked pairs (augm. RL)   294.00 291.00 288.00 286.00 283.00 280.00 277.00
#> Linked pairs (RL)         302.00 298.00 294.00 292.00 287.00 284.00 281.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.10   0.09   0.09   0.09   0.08   0.08   0.07
#> FDP synth data estimator    0.13   0.11   0.11   0.11   0.11   0.12   0.10
#> Linked pairs (augm. RL)   275.00 271.00 268.00 266.00 262.00 261.00 256.00
#> Linked pairs (RL)         278.00 274.00 271.00 268.00 265.00 264.00 258.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.07   0.07   0.07   0.07   0.07   0.06   0.06
#> FDP synth data estimator    0.08   0.08   0.06   0.06   0.06   0.06   0.06
#> Linked pairs (augm. RL)   254.00 252.00 251.00 250.00 250.00 248.00 248.00
#> Linked pairs (RL)         256.00 254.00 252.00 252.00 252.00 250.00 250.00
#>                             0.71   0.72   0.73   0.74   0.75   0.76   0.77
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.05   0.05   0.05
#> FDP synth data estimator    0.04   0.04   0.04   0.04   0.04   0.04   0.04
#> Linked pairs (augm. RL)   246.00 242.00 241.00 240.00 236.00 235.00 230.00
#> Linked pairs (RL)         247.00 244.00 242.00 241.00 237.00 236.00 232.00
#>                             0.78   0.79   0.80   0.81   0.82   0.83   0.84
#> FDP model score estimator   0.05   0.05   0.04   0.04   0.04   0.04   0.04
#> FDP synth data estimator    0.04   0.04   0.04   0.04   0.05   0.05   0.05
#> Linked pairs (augm. RL)   228.00 226.00 224.00 222.00 220.00 218.00 215.00
#> Linked pairs (RL)         230.00 228.00 226.00 224.00 221.00 219.00 216.00
#>                             0.85   0.86   0.87   0.88   0.89   0.90   0.91
#> FDP model score estimator   0.04   0.03   0.03   0.03   0.03   0.02   0.02
#> FDP synth data estimator    0.05   0.05   0.05   0.05   0.05   0.05   0.03
#> Linked pairs (augm. RL)   210.00 206.00 203.00 196.00 194.00 186.00 182.00
#> Linked pairs (RL)         212.00 208.00 204.00 198.00 194.00 187.00 182.00
#>                             0.92   0.93   0.94   0.95   0.96   0.97   0.98 0.99
#> FDP model score estimator   0.02   0.02   0.02   0.02   0.01   0.01   0.01    0
#> FDP synth data estimator    0.03   0.00   0.00   0.00   0.00   0.00   0.00    0
#> Linked pairs (augm. RL)   179.00 170.00 163.00 158.00 147.00 134.00 108.00   69
#> Linked pairs (RL)         180.00 170.00 163.00 158.00 147.00 134.00 108.00   69
#> FlexRL
#> Running StEM algorithm ■■■■■                             13% | iter 4/30 [1.6s]
#> Running StEM algorithm ■■■■■■■■■■■■                      37% | iter 11/30 [4.7s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■                57% | iter 17/30 [7.6s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■         80% | iter 24/30 [10.9…
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 30/30 [13.6…
#> 
#> Drawing Delta          ■■                                 3% [3s]
#> Drawing Delta          ■■■                                6% [6.1s]
#> Drawing Delta          ■■■■                               9% [9s]
#> Drawing Delta          ■■■■■                             12% [12.1s]
#> Drawing Delta          ■■■■■                             15% [15s]
#> Drawing Delta          ■■■■■■                            17% [18s]
#> Drawing Delta          ■■■■■■■                           20% [21s]
#> Drawing Delta          ■■■■■■■■                          23% [24s]
#> Drawing Delta          ■■■■■■■■■                         26% [27s]
#> Drawing Delta          ■■■■■■■■■■                        29% [30.1s]
#> Drawing Delta          ■■■■■■■■■■■                       32% [33s]
#> Drawing Delta          ■■■■■■■■■■■                       35% [36s]
#> Drawing Delta          ■■■■■■■■■■■■                      38% [39s]
#> Drawing Delta          ■■■■■■■■■■■■■                     41% [42.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    44% [45s]
#> Drawing Delta          ■■■■■■■■■■■■■■■                   46% [48s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  49% [51s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■                 52% [54.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                55% [57s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                58% [1m]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               61% [1m 3s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              64% [1m 6s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■             67% [1m 9s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            70% [1m 12s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■           73% [1m 15.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■          76% [1m 18s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         79% [1m 21.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         82% [1m 24s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% [1m 27s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% [1m 30s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% [1m 33s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% [1m 36.1s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% [1m 39s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% [1m 42s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% [1m 42.7s]
#> 
#> FlexRL
#> Running StEM algorithm ■■■■■■■                           20% | iter 6/30 [2.6s]
#> Running StEM algorithm ■■■■■■■■■■■■■■                    43% | iter 13/30 [5.6s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■              63% | iter 19/30 [8.3s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% | iter 27/30 [11.6…
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 30/30 [12.8…
#> 
#> Drawing Delta          ■                                  1% [1.4s]
#> Drawing Delta          ■■                                 4% [4.4s]
#> Drawing Delta          ■■■                                7% [7.5s]
#> Drawing Delta          ■■■■                              10% [10.5s]
#> Drawing Delta          ■■■■■                             13% [13.5s]
#> Drawing Delta          ■■■■■■                            16% [16.4s]
#> Drawing Delta          ■■■■■■■                           19% [19.4s]
#> Drawing Delta          ■■■■■■■■                          22% [22.4s]
#> Drawing Delta          ■■■■■■■■■                         25% [25.4s]
#> Drawing Delta          ■■■■■■■■■                         28% [28.4s]
#> Drawing Delta          ■■■■■■■■■■                        31% [31.5s]
#> Drawing Delta          ■■■■■■■■■■■                       34% [34.5s]
#> Drawing Delta          ■■■■■■■■■■■■                      37% [37.4s]
#> Drawing Delta          ■■■■■■■■■■■■■                     40% [40.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■                    43% [43.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■                   46% [46.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■                  50% [49.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■                 52% [52.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                55% [55.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■                58% [58.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■               61% [1m 1.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■              64% [1m 4.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■             67% [1m 7.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■            70% [1m 10.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■           73% [1m 13.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■          76% [1m 16.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■         79% [1m 19.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% [1m 22.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       85% [1m 25.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% [1m 28.5s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      91% [1m 31.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% [1m 34.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% [1m 37.4s]
#> Drawing Delta          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% [1m 40.2s]
#> 
#> FlexRL results (average over 2 iterations):
#>                            0.50  0.51  0.52   0.53   0.54   0.55   0.56   0.57
#> FDP model score estimator   0.1   0.1   0.1   0.09   0.09   0.09   0.09   0.09
#> Linked pairs (RL)         306.0 304.0 302.0 298.00 297.00 295.00 294.00 293.00
#>                             0.58   0.59   0.60   0.61   0.62   0.63   0.64
#> FDP model score estimator   0.09   0.08   0.08   0.07   0.07   0.07   0.07
#> Linked pairs (RL)         292.00 288.00 285.00 281.00 279.00 276.00 274.00
#>                             0.65   0.66   0.67   0.68   0.69   0.70   0.71
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.06   0.06   0.05
#> Linked pairs (RL)         273.00 272.00 269.00 268.00 266.00 265.00 262.00
#>                             0.72   0.73   0.74   0.75   0.76   0.77   0.78
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.05   0.05
#> Linked pairs (RL)         260.00 260.00 258.00 256.00 254.00 252.00 252.00
#>                             0.79   0.80   0.81   0.82   0.83   0.84   0.85
#> FDP model score estimator   0.04   0.04   0.04   0.04   0.04   0.03   0.03
#> Linked pairs (RL)         250.00 248.00 245.00 240.00 236.00 232.00 228.00
#>                             0.86   0.87   0.88   0.89   0.90   0.91   0.92
#> FDP model score estimator   0.03   0.03   0.03   0.02   0.02   0.02   0.02
#> Linked pairs (RL)         222.00 219.00 214.00 210.00 205.00 201.00 197.00
#>                             0.93   0.94   0.95   0.96   0.97   0.98 0.99
#> FDP model score estimator   0.02   0.02   0.01   0.01   0.01   0.01    0
#> Linked pairs (RL)         192.00 185.00 178.00 162.00 154.00 130.00   90

df_results[6,"FlexRL    "] = diag$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"FlexRL    "] = diag$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"FlexRL    "] = diag$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag$discrepancy_measures))
iou_idx = grep("iou", names(diag$discrepancy_measures))
mmd_idx = grep("mmd", names(diag$discrepancy_measures))

df_results[9, "FlexRL    "] = max(diag$discrepancy_measures[1,smd_idx])
df_results[10,"FlexRL    "] = min(diag$discrepancy_measures[1,iou_idx])
df_results[11,"FlexRL    "] = max(diag$discrepancy_measures[1,mmd_idx])
```

``` r

diag # print(diag)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  316
#> 
#>   FDP           (true pairs):    0.146
#>   Sensitivity   (true pairs):    0.675
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.103
#>                                                threshold 0.50: FDP ~ 0.103
#>                                            max threshold 0.99: FDP ~ 0.004
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.125
#>                                                threshold 0.50: FDP ~ 0.125
#>                                            max threshold 0.99: FDP ~ 0.004
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.238
#>                                                threshold 0.50: FDP ~ 0.238
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0058
#>   Multivariate MMD  (linked B vs. B):  0.0035
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  1.0000
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  0.8750
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9810
#>   Agreement V2  (linked A vs. linked B):  0.9905
#>   Agreement V3  (linked A vs. linked B):  0.9905
#>   Agreement V4  (linked A vs. linked B):  0.7943
#>   Agreement V1  (true pairs):  0.9425
#>   Agreement V2  (true pairs):  0.9450
#>   Agreement V3  (true pairs):  0.9525
#>   Agreement V4  (true pairs):  0.7150
#> 
#>   SMD V1 values: {2, 6, 8, ...}  (linked A vs. A):  0.0405, -0.0663, -0.0556, ...
#>   SMD V2 values: {2, 6, 9, ...}  (linked A vs. A):  0.0552, 0.0448, -0.1255, ...
#>   SMD V3 values: {3, 5, 8, ...}  (linked A vs. A):  -0.1198, 0.0787, 0.0358, ...
#>   SMD V4                         (linked A vs. A):  -0.0221
#>   SMD V1 values: {1, 6, 8, ...}   (linked B vs. B):  0.0538, -0.0627, -0.0657, ...
#>   SMD V2 values: {2, 5, 9, ...}   (linked B vs. B):  0.0551, 0.0932, -0.1197, ...
#>   SMD V3 values: {3, 9, 10, ...}  (linked B vs. B):  -0.0602, 0.0529, -0.0370, ...
#>   SMD V4                          (linked B vs. B):  0.0192

print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  253
#> 
#>   FDP           (true pairs):    0.059
#>   Sensitivity   (true pairs):    0.595
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.103
#>                                                threshold 0.75: FDP ~ 0.049
#>                                            max threshold 0.99: FDP ~ 0.004
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.125
#>                                                threshold 0.75: FDP ~ 0.053
#>                                            max threshold 0.99: FDP ~ 0.004
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.238
#>                                                threshold 0.75: FDP ~ 0.042
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0098
#>   Multivariate MMD  (linked B vs. B):  0.0063
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  0.8889
#>   IoU support V3  (linked A vs. A):  0.9000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  1.0000
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9842
#>   Agreement V2  (linked A vs. linked B):  0.9960
#>   Agreement V3  (linked A vs. linked B):  0.9960
#>   Agreement V4  (linked A vs. linked B):  0.8458
#>   Agreement V1  (true pairs):  0.9425
#>   Agreement V2  (true pairs):  0.9450
#>   Agreement V3  (true pairs):  0.9525
#>   Agreement V4  (true pairs):  0.7150
#> 
#>   SMD V1 values: {1, 4, 8, ...}  (linked A vs. A):  0.0755, 0.0595, -0.1132, ...
#>   SMD V2 values: {0, 2, 9, ...}  (linked A vs. A):  -0.0818, 0.1088, -0.1636, ...
#>   SMD V3 values: {3, 4, 5, ...}  (linked A vs. A):  -0.2055, 0.0784, 0.0972, ...
#>   SMD V4                         (linked A vs. A):  -0.0608
#>   SMD V1 values: {1, 2, 8, ...}  (linked B vs. B):  0.0968, 0.0562, -0.1299, ...
#>   SMD V2 values: {2, 5, 9, ...}  (linked B vs. B):  0.0933, 0.1332, -0.1674, ...
#>   SMD V3 values: {1, 3, 6, ...}  (linked B vs. B):  0.0614, -0.1520, 0.1039, ...
#>   SMD V4                         (linked B vs. B):  -0.0081
```

These summaries provide information on the linked data. We can also look
at: the convergence of the StEM algorithm parameters, the linkage scores
distribution, the data distributions between linked subset and original
data sources, the divergence metrics of the linked data from the data
sources, the FDP estimation.

``` r

plot(diag, "convergence")
```

![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-1.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-2.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-3.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-4.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-5.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-6.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-7.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-8.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-9.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-10.png)

``` r


plot(diag, "scores")
```

![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-11.png)

``` r

 
plot(diag, "distributions", threshold = 0.75)
```

![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-12.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-13.png)

``` r


plot(diag, "discrepancy")
```

![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-14.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-15.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-16.png)![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-17.png)

``` r


plot(diag, "FDP")
```

![](FlexRL-RL-vignette_files/figure-html/unnamed-chunk-14-18.png)

Several other algorithms for **R**ecord **L**inkage have been developed
so far in the literature. Each has its own qualities and flaws.
[`StEM()`](../reference/stEM.md) usually outperforms in scenarios where
the **PIVs** have low discriminative power with few expected
registration errors; it is particularly efficient when there is enough
information to model the dynamics of some **PIV** changing across time
(such as postal code). In 2026, it is the only package taking account of
dynamic PIVs. More can be found on this topic in the simulation setting
of the [methodology paper](https://doi.org/10.1093/jrsssc/qlaf016), code
is provided in
[`FlexRL-experiments`](https://github.com/robachowyk/FlexRL-experiments).
More can be found on the FDP estimation methods in [this
paper](https://doi.org/10.1002/sim.70292).

With \[link_with\_\*()\] wrappers we can link records with “multilink”,
“fastLink”, “BRL”, “reclin2”, “diyar”, “fedmatch”.

``` r

brl_args <- list(flds = PIVs, types = rep("bi",length(PIVs)))
fit_brl <- link_with_BRL(prep_data$encodedA, prep_data$encodedB, brl_args)

delta_result = as.data.frame(fit_brl)
delta_result = delta_result[delta_result$LinkScore>0.5, ]
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"BRL       "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

diag_brl <- do.call(RL_diagnostics, c(list(fit_brl, prep_data$encodedA, prep_data$encodedB,
                           PIVs, PIVs_type, true_pairs = prep_data$true_pairs,
                           FDP_estimation = TRUE, RL_method = "BRL", maxIter4CV = 1, n_repeats = 2),
                           brl_args))
#> Iteration: 0, Accuracy: 49.54%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 49.54%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> BRL results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.15   0.13   0.12   0.12   0.12   0.11   0.10
#> FDP synth data estimator    0.28   0.29   0.30   0.26   0.26   0.22   0.22
#> Linked pairs (augm. RL)   251.00 242.00 236.00 235.00 233.00 228.00 223.00
#> Linked pairs (RL)         258.00 249.00 243.00 241.00 239.00 233.00 228.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.09   0.09   0.08   0.07   0.07   0.07   0.07
#> FDP synth data estimator    0.19   0.19   0.10   0.10   0.10   0.10   0.10
#> Linked pairs (augm. RL)   216.00 215.00 210.00 207.00 207.00 204.00 204.00
#> Linked pairs (RL)         220.00 219.00 212.00 209.00 209.00 206.00 206.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.07   0.07   0.07   0.07   0.07   0.07   0.07
#> FDP synth data estimator    0.10   0.10   0.10   0.10   0.10   0.10   0.10
#> Linked pairs (augm. RL)   204.00 204.00 204.00 204.00 204.00 204.00 204.00
#> Linked pairs (RL)         206.00 206.00 206.00 206.00 206.00 206.00 206.00
#>                             0.71   0.72   0.73   0.74   0.75   0.76   0.77
#> FDP model score estimator   0.07   0.07   0.07   0.07   0.07   0.07   0.07
#> FDP synth data estimator    0.10   0.10   0.10   0.10   0.10   0.10   0.10
#> Linked pairs (augm. RL)   204.00 204.00 204.00 204.00 203.00 203.00 203.00
#> Linked pairs (RL)         206.00 206.00 206.00 206.00 205.00 205.00 205.00
#>                             0.78   0.79   0.80   0.81   0.82   0.83   0.84
#> FDP model score estimator   0.07   0.07   0.07   0.07   0.07   0.06   0.06
#> FDP synth data estimator    0.10   0.10   0.10   0.10   0.10   0.10   0.10
#> Linked pairs (augm. RL)   202.00 202.00 202.00 201.00 200.00 198.00 191.00
#> Linked pairs (RL)         204.00 204.00 204.00 203.00 202.00 200.00 193.00
#>                             0.85   0.86   0.87   0.88   0.89   0.90   0.91
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.05   0.05   0.05
#> FDP synth data estimator    0.05   0.05   0.05   0.06   0.06   0.00   0.00
#> Linked pairs (augm. RL)   190.00 187.00 183.00 181.00 173.00 167.00 151.00
#> Linked pairs (RL)         191.00 188.00 184.00 182.00 174.00 167.00 151.00
#>                             0.92   0.93   0.94  0.95  0.96  0.97 0.98 0.99
#> FDP model score estimator   0.04   0.04   0.04  0.04  0.03  0.03 0.02    0
#> FDP synth data estimator    0.00   0.00   0.00  0.00  0.00  0.00 0.00    0
#> Linked pairs (augm. RL)   142.00 128.00 114.00 98.00 69.00 15.00 1.00    0
#> Linked pairs (RL)         142.00 128.00 114.00 98.00 69.00 15.00 1.00    0
#> BRL results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.14   0.13   0.13   0.13   0.12   0.11   0.11
#> Linked pairs (RL)         267.00 265.00 262.00 258.00 254.00 250.00 245.00
#>                             0.57  0.58  0.59  0.60   0.61   0.62   0.63   0.64
#> FDP model score estimator   0.11   0.1   0.1   0.1   0.09   0.08   0.07   0.07
#> Linked pairs (RL)         243.00 242.0 240.0 236.0 230.00 222.00 218.00 213.00
#>                             0.65   0.66   0.67   0.68   0.69   0.70   0.71
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.06   0.06   0.06
#> Linked pairs (RL)         211.00 210.00 208.00 208.00 208.00 208.00 208.00
#>                             0.72   0.73   0.74   0.75   0.76   0.77   0.78
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.06   0.06   0.06
#> Linked pairs (RL)         208.00 208.00 208.00 208.00 208.00 208.00 208.00
#>                             0.79   0.80   0.81   0.82   0.83   0.84   0.85
#> FDP model score estimator   0.06   0.06   0.06   0.06   0.06   0.06   0.06
#> Linked pairs (RL)         208.00 207.00 207.00 206.00 206.00 204.00 201.00
#>                             0.86   0.87   0.88   0.89   0.90   0.91   0.92
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.04   0.04
#> Linked pairs (RL)         196.00 193.00 190.00 183.00 177.00 169.00 154.00
#>                             0.93   0.94   0.95   0.96  0.97 0.98 0.99
#> FDP model score estimator   0.04   0.03   0.03   0.03  0.02 0.02    0
#> Linked pairs (RL)         139.00 125.00 114.00 101.00 52.00 4.00    0

df_results[6,"BRL       "] = diag_brl$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"BRL       "] = diag_brl$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"BRL       "] = diag_brl$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_brl$discrepancy_measures))
iou_idx = grep("iou", names(diag_brl$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_brl$discrepancy_measures))

df_results[9, "BRL       "] = max(diag_brl$discrepancy_measures[1,smd_idx])
df_results[10,"BRL       "] = min(diag_brl$discrepancy_measures[1,iou_idx])
df_results[11,"BRL       "] = max(diag_brl$discrepancy_measures[1,mmd_idx])
```

``` r

multilink_args <- list(types = rep("bi",length(PIVs)), breaks = list(V1 = NA, V2 = NA, V3 = NA, V4 = NA),
                       duplicates = c(0,0), n_iter = 500, max_cc_size = 50)
fit_multilink <- link_with_multilink( prep_data$encodedA, prep_data$encodedB, multilink_args)
#> [1] "Creating comparison data for field V1"
#> [1] "Creating comparison data for field V2"
#> [1] "Creating comparison data for field V3"
#> [1] "Creating comparison data for field V4"
#> [1] "Running Gibbs sampler with Gibbs updates to partition."
#> [1] "Processing inputs"
#> [1] "Beginning sampling"
#> Beginning iteration 1/500
#> Beginning iteration 100/500
#> Beginning iteration 200/500
#> Beginning iteration 300/500
#> Beginning iteration 400/500
#> Beginning iteration 500/500
#> [1] "Finished sampling in 29.7413575649261 seconds"
#> [1] "Finding Bayes estimate with a threshold of 0.0375 and a maximum connected component of size 41"

delta_result = as.data.frame(fit_multilink[c("idxA", "idxB")])
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"multilink "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

diag_multilink <- do.call(RL_diagnostics, c(list(fit_multilink, prep_data$encodedA, prep_data$encodedB,
                                PIVs, PIVs_type, true_pairs = prep_data$true_pairs,
                                FDP_estimation = TRUE, RL_method = "multilink", maxIter4CV = 1, n_repeats = 2),
                                multilink_args))
#> Iteration: 0, Accuracy: 50.61%
#> Iteration: 1, Accuracy: 41.41%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> [1] "Creating comparison data for field V1"
#> [1] "Creating comparison data for field V2"
#> [1] "Creating comparison data for field V3"
#> [1] "Creating comparison data for field V4"
#> [1] "Running Gibbs sampler with Gibbs updates to partition."
#> [1] "Processing inputs"
#> [1] "Beginning sampling"
#> Beginning iteration 1/500
#> Beginning iteration 100/500
#> Beginning iteration 200/500
#> Beginning iteration 300/500
#> Beginning iteration 400/500
#> Beginning iteration 500/500
#> [1] "Finished sampling in 34.8073935508728 seconds"
#> [1] "Finding Bayes estimate with a threshold of 0.0375 and a maximum connected component of size 30"
#> Iteration: 0, Accuracy: 51.32%
#> Iteration: 1, Accuracy: 43.61%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> [1] "Creating comparison data for field V1"
#> [1] "Creating comparison data for field V2"
#> [1] "Creating comparison data for field V3"
#> [1] "Creating comparison data for field V4"
#> [1] "Running Gibbs sampler with Gibbs updates to partition."
#> [1] "Processing inputs"
#> [1] "Beginning sampling"
#> Beginning iteration 1/500
#> Beginning iteration 100/500
#> Beginning iteration 200/500
#> Beginning iteration 300/500
#> Beginning iteration 400/500
#> Beginning iteration 500/500
#> [1] "Finished sampling in 33.212367773056 seconds"
#> [1] "Finding Bayes estimate with a threshold of 0.035 and a maximum connected component of size 43"
#> multilink results (average over 2 iterations):
#>                             0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#> FDP model score estimator    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator    0.17  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)   266.00  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)            NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> [1] "Creating comparison data for field V1"
#> [1] "Creating comparison data for field V2"
#> [1] "Creating comparison data for field V3"
#> [1] "Creating comparison data for field V4"
#> [1] "Running Gibbs sampler with Gibbs updates to partition."
#> [1] "Processing inputs"
#> [1] "Beginning sampling"
#> Beginning iteration 1/500
#> Beginning iteration 100/500
#> Beginning iteration 200/500
#> Beginning iteration 300/500
#> Beginning iteration 400/500
#> Beginning iteration 500/500
#> [1] "Finished sampling in 29.5847890377045 seconds"
#> [1] "Finding Bayes estimate with a threshold of 0.0375 and a maximum connected component of size 41"
#> Warning: No valid FDP estimate after 1 attempts at iteration 1. Increase
#> `maxIter4CV`, or the estimator may be unreliable for this method/data.
#> multilink results (average over 2 iterations):
#>                           0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN

df_results[6,"multilink "] = diag_multilink$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"multilink "] = diag_multilink$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"multilink "] = diag_multilink$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_multilink$discrepancy_measures))
iou_idx = grep("iou", names(diag_multilink$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_multilink$discrepancy_measures))

df_results[9, "multilink "] = max(diag_multilink$discrepancy_measures[1,smd_idx])
df_results[10,"multilink "] = min(diag_multilink$discrepancy_measures[1,iou_idx])
df_results[11,"multilink "] = max(diag_multilink$discrepancy_measures[1,mmd_idx])
```

``` r

fastLink_args <- list(varnames = PIVs, tol.em = 1e-06, threshold.match = 0.5, return.all = TRUE, n.cores = 1 )
fit_fastLink <- link_with_fastLink( prep_data$encodedA, prep_data$encodedB, fastLink_args )
#> 
#> ==================== 
#> fastLink(): Fast Probabilistic Record Linkage
#> ==================== 
#> 
#> Calculating matches for each variable.
#> Getting counts for parameter estimation.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Running the EM algorithm.
#> Iteration number 100 
#> Maximum difference in log-likelihood = 0.0029 
#> Iteration number 200 
#> Maximum difference in log-likelihood = 4e-04 
#> Iteration number 300 
#> Maximum difference in log-likelihood = 1e-04 
#> Iteration number 400 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 500 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 600 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 700 
#> Maximum difference in log-likelihood = 0 
#> Getting the indices of estimated matches.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Deduping the estimated matches.
#> Getting the match patterns for each estimated match.

delta_result = as.data.frame(fit_fastLink)
delta_result = delta_result[delta_result$LinkScore>0.5, ]
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"fastLink  "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

diag_fastLink <- do.call(RL_diagnostics, c(list(fit_fastLink, prep_data$encodedA, prep_data$encodedB,
                                PIVs, PIVs_type, true_pairs = prep_data$true_pairs,
                                FDP_estimation = TRUE, RL_method = "fastLink", maxIter4CV = 1, n_repeats = 2),
                                fastLink_args))
#> Iteration: 0, Accuracy: 51.87%
#> Iteration: 1, Accuracy: 40.85%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> 
#> ==================== 
#> fastLink(): Fast Probabilistic Record Linkage
#> ==================== 
#> 
#> Calculating matches for each variable.
#> Getting counts for parameter estimation.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Running the EM algorithm.
#> Iteration number 100 
#> Maximum difference in log-likelihood = 0.0026 
#> Iteration number 200 
#> Maximum difference in log-likelihood = 5e-04 
#> Iteration number 300 
#> Maximum difference in log-likelihood = 1e-04 
#> Iteration number 400 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 500 
#> Maximum difference in log-likelihood = 0 
#> Getting the indices of estimated matches.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Deduping the estimated matches.
#> Getting the match patterns for each estimated match.
#> Iteration: 0, Accuracy: 46.6%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> 
#> ==================== 
#> fastLink(): Fast Probabilistic Record Linkage
#> ==================== 
#> 
#> Calculating matches for each variable.
#> Getting counts for parameter estimation.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Running the EM algorithm.
#> Iteration number 100 
#> Maximum difference in log-likelihood = 0.002 
#> Iteration number 200 
#> Maximum difference in log-likelihood = 3e-04 
#> Iteration number 300 
#> Maximum difference in log-likelihood = 1e-04 
#> Iteration number 400 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 500 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 600 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 700 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 800 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 900 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 1000 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 1100 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 1200 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 1300 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 1400 
#> Maximum difference in log-likelihood = 0 
#> Getting the indices of estimated matches.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Deduping the estimated matches.
#> Getting the match patterns for each estimated match.
#> fastLink results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.24   0.24   0.24   0.24   0.24   0.24   0.24
#> FDP synth data estimator    0.20   0.20   0.20   0.20   0.20   0.20   0.20
#> Linked pairs (augm. RL)   272.00 272.00 272.00 272.00 272.00 272.00 272.00
#> Linked pairs (RL)         277.00 277.00 277.00 277.00 277.00 277.00 277.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.24   0.24   0.24   0.24   0.24   0.24   0.24
#> FDP synth data estimator    0.20   0.20   0.20   0.20   0.20   0.20   0.20
#> Linked pairs (augm. RL)   272.00 272.00 272.00 272.00 272.00 272.00 272.00
#> Linked pairs (RL)         277.00 277.00 277.00 277.00 277.00 277.00 277.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.24   0.24   0.24   0.24   0.24   0.24   0.24
#> FDP synth data estimator    0.20   0.20   0.20   0.20   0.20   0.20   0.20
#> Linked pairs (augm. RL)   272.00 272.00 272.00 272.00 272.00 272.00 272.00
#> Linked pairs (RL)         277.00 277.00 277.00 277.00 277.00 277.00 277.00
#>                             0.71   0.72   0.73   0.74   0.75   0.76 0.77 0.78
#> FDP model score estimator   0.24   0.24   0.24   0.24   0.24   0.24    0    0
#> FDP synth data estimator    0.20   0.20   0.20   0.20   0.20   0.20    0    0
#> Linked pairs (augm. RL)   272.00 272.00 272.00 272.00 272.00 272.00    0    0
#> Linked pairs (RL)         277.00 277.00 277.00 277.00 277.00 277.00    0    0
#>                           0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.99
#> FDP model score estimator    0
#> FDP synth data estimator     0
#> Linked pairs (augm. RL)      0
#> Linked pairs (RL)            0
#> 
#> ==================== 
#> fastLink(): Fast Probabilistic Record Linkage
#> ==================== 
#> 
#> Calculating matches for each variable.
#> Getting counts for parameter estimation.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Running the EM algorithm.
#> Iteration number 100 
#> Maximum difference in log-likelihood = 0.0027 
#> Iteration number 200 
#> Maximum difference in log-likelihood = 4e-04 
#> Iteration number 300 
#> Maximum difference in log-likelihood = 1e-04 
#> Iteration number 400 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 500 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 600 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 700 
#> Maximum difference in log-likelihood = 0 
#> Getting the indices of estimated matches.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Deduping the estimated matches.
#> Getting the match patterns for each estimated match.
#> 
#> ==================== 
#> fastLink(): Fast Probabilistic Record Linkage
#> ==================== 
#> 
#> Calculating matches for each variable.
#> Getting counts for parameter estimation.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Running the EM algorithm.
#> Iteration number 100 
#> Maximum difference in log-likelihood = 0.0022 
#> Iteration number 200 
#> Maximum difference in log-likelihood = 3e-04 
#> Iteration number 300 
#> Maximum difference in log-likelihood = 1e-04 
#> Iteration number 400 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 500 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 600 
#> Maximum difference in log-likelihood = 0 
#> Iteration number 700 
#> Maximum difference in log-likelihood = 0 
#> Getting the indices of estimated matches.
#>     Parallelizing calculation using OpenMP. 1 threads out of 4 are used.
#> Deduping the estimated matches.
#> Getting the match patterns for each estimated match.
#> fastLink results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.23   0.23   0.23   0.23   0.23   0.23   0.23
#> Linked pairs (RL)         275.00 275.00 275.00 275.00 275.00 275.00 275.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.23   0.23   0.23   0.23   0.23   0.23   0.23
#> Linked pairs (RL)         275.00 275.00 275.00 275.00 275.00 275.00 275.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.23   0.23   0.23   0.23   0.23   0.23   0.23
#> Linked pairs (RL)         275.00 275.00 275.00 275.00 275.00 275.00 275.00
#>                             0.71   0.72   0.73   0.74   0.75   0.76   0.77 0.78
#> FDP model score estimator   0.23   0.23   0.23   0.23   0.23   0.23   0.23    0
#> Linked pairs (RL)         275.00 275.00 275.00 275.00 275.00 275.00 275.00    0
#>                           0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.99
#> FDP model score estimator    0
#> Linked pairs (RL)            0

df_results[6,"fastLink  "] = diag_fastLink$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"fastLink  "] = diag_fastLink$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"fastLink  "] = diag_fastLink$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_fastLink$discrepancy_measures))
iou_idx = grep("iou", names(diag_fastLink$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_fastLink$discrepancy_measures))

df_results[9, "fastLink  "] = max(diag_fastLink$discrepancy_measures[1,smd_idx])
df_results[10,"fastLink  "] = min(diag_fastLink$discrepancy_measures[1,iou_idx])
df_results[11,"fastLink  "] = max(diag_fastLink$discrepancy_measures[1,mmd_idx])
```

``` r

reclin2_args <- list(on = PIVs, formula = stats::as.formula(paste("~", paste(PIVs, collapse = "+"))),
                     type = "mpost", add = TRUE, variable = "selected", score = "mpost", threshold = 0.5)
fit_reclin2 <- link_with_reclin2( prep_data$encodedA, prep_data$encodedB, reclin2_args )

delta_result = as.data.frame(fit_reclin2)
delta_result = delta_result[delta_result$LinkScore>0.5, ]
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"reclin2   "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

diag_reclin2 <- do.call(RL_diagnostics, c(list(fit_reclin2, prep_data$encodedA, prep_data$encodedB, 
                               PIVs, PIVs_type, true_pairs = prep_data$true_pairs, 
                               FDP_estimation = TRUE, RL_method = "reclin2", maxIter4CV = 1, n_repeats = 2),
                               reclin2_args))
#> Iteration: 0, Accuracy: 51.62%
#> Iteration: 1, Accuracy: 40.68%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 48.28%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> reclin2 results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.31   0.31   0.31   0.31   0.31   0.31   0.31
#> FDP synth data estimator    0.30   0.30   0.30   0.30   0.30   0.30   0.30
#> Linked pairs (augm. RL)   312.00 312.00 312.00 312.00 312.00 312.00 312.00
#> Linked pairs (RL)         322.00 322.00 322.00 322.00 322.00 322.00 322.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.31   0.31   0.31   0.31   0.31   0.31   0.31
#> FDP synth data estimator    0.30   0.30   0.30   0.30   0.30   0.30   0.30
#> Linked pairs (augm. RL)   312.00 312.00 312.00 312.00 312.00 312.00 312.00
#> Linked pairs (RL)         322.00 322.00 322.00 322.00 322.00 322.00 322.00
#>                             0.64   0.65   0.66   0.67   0.68 0.69 0.70 0.71
#> FDP model score estimator   0.31   0.31   0.31   0.31   0.31    0    0    0
#> FDP synth data estimator    0.30   0.30   0.30   0.30   0.30    0    0    0
#> Linked pairs (augm. RL)   312.00 312.00 312.00 312.00 312.00    0    0    0
#> Linked pairs (RL)         322.00 322.00 322.00 322.00 322.00    0    0    0
#>                           0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0
#> reclin2 results (average over 2 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator   0.3   0.3   0.3   0.3   0.3   0.3   0.3   0.3   0.3
#> Linked pairs (RL)         312.0 312.0 312.0 312.0 312.0 312.0 312.0 312.0 312.0
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator   0.3   0.3   0.3   0.3   0.3   0.3   0.3   0.3   0.3
#> Linked pairs (RL)         312.0 312.0 312.0 312.0 312.0 312.0 312.0 312.0 312.0
#>                            0.68  0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77
#> FDP model score estimator   0.3   0.3    0    0    0    0    0    0    0    0
#> Linked pairs (RL)         312.0 312.0    0    0    0    0    0    0    0    0
#>                           0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.98 0.99
#> FDP model score estimator    0    0
#> Linked pairs (RL)            0    0

df_results[6,"reclin2   "] = diag_reclin2$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"reclin2   "] = diag_reclin2$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"reclin2   "] = diag_reclin2$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_reclin2$discrepancy_measures))
iou_idx = grep("iou", names(diag_reclin2$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_reclin2$discrepancy_measures))

df_results[9, "reclin2   "] = max(diag_reclin2$discrepancy_measures[1,smd_idx])
df_results[10,"reclin2   "] = min(diag_reclin2$discrepancy_measures[1,iou_idx])
df_results[11,"reclin2   "] = max(diag_reclin2$discrepancy_measures[1,mmd_idx])
```

``` r

diyar_args <- list(probabilistic = TRUE, return_weights = TRUE)
fit_diyar <- link_with_diyar(prep_data$encodedA[, PIVs], prep_data$encodedB[, PIVs], diyar_args)

delta_result = as.data.frame(fit_diyar[c("idxA", "idxB")])
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"diyar     "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

names4diag <- c(PIVs, "local_id", "source")
diag_diyar <- do.call(RL_diagnostics, c(list(fit_diyar, 
                             prep_data$encodedA[, names4diag], prep_data$encodedB[, names4diag], 
                             PIVs, PIVs_type, true_pairs = prep_data$true_pairs, 
                             FDP_estimation = TRUE, RL_method = "diyar", maxIter4CV = 1, n_repeats = 2),
                             diyar_args))
#> Iteration: 0, Accuracy: 48.93%
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 48.43%
#> diyar results (average over 2 iterations):
#>                           0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator     0  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)     86  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Warning: No valid FDP estimate after 1 attempts at iteration 1. Increase
#> `maxIter4CV`, or the estimator may be unreliable for this method/data.
#> diyar results (average over 2 iterations):
#>                           0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (RL)          NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN

df_results[6,"diyar     "] = diag_diyar$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"diyar     "] = diag_diyar$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"diyar     "] = diag_diyar$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_diyar$discrepancy_measures))
iou_idx = grep("iou", names(diag_diyar$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_diyar$discrepancy_measures))

df_results[9, "diyar     "] = max(diag_diyar$discrepancy_measures[1,smd_idx])
df_results[10,"diyar     "] = min(diag_diyar$discrepancy_measures[1,iou_idx])
df_results[11,"diyar     "] = max(diag_diyar$discrepancy_measures[1,mmd_idx])
```

``` r

fedmatch_args <- list(by = PIVs, match_type = "multivar", 
                      unique_key_1 = "unique_key_A", unique_key_2 = "unique_key_B",
                      multivar_settings = fedmatch::build_multivar_settings(
                        compare_type = rep("indicator", length(PIVs)),
                        wgts = rep(1/length(PIVs), length(PIVs))))
fit_fedmatch <- link_with_fedmatch( prep_data$encodedA, prep_data$encodedB, fedmatch_args )

delta_result = as.data.frame(fit_fedmatch)
delta_result = delta_result[delta_result$LinkScore>0.5, ]
linked_pairs = do.call(paste, c(delta_result[,c("idxA","idxB")], list(sep = "_")))
true_positive = length( intersect(linked_pairs, true_pairs) )
false_positive = length( setdiff(linked_pairs, true_pairs) )
false_negative = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive) 

df_results[1:5,"fedmatch  "] = c(true_positive,false_positive,false_negative,sensitivity,fdp)

diag_fedmatch <- do.call(RL_diagnostics, c(list(fit_fedmatch, prep_data$encodedA, prep_data$encodedB, 
                               PIVs, PIVs_type, true_pairs = prep_data$true_pairs, 
                               FDP_estimation = TRUE, RL_method = "fedmatch", maxIter4CV = 1, n_repeats = 2),
                               fedmatch_args))
#> Iteration: 0, Accuracy: 46.51%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Warning: No valid FDP estimate after 1 attempts at iteration 1. Increase
#> `maxIter4CV`, or the estimator may be unreliable for this method/data.
#> fedmatch results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> FDP synth data estimator     NaN    NaN    NaN    NaN    NaN    NaN    NaN
#> Linked pairs (augm. RL)   348.00 348.00 348.00 348.00 348.00 348.00 348.00
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> FDP synth data estimator     NaN    NaN    NaN    NaN    NaN    NaN    NaN
#> Linked pairs (augm. RL)   348.00 348.00 348.00 348.00 348.00 348.00 348.00
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> FDP synth data estimator     NaN    NaN    NaN    NaN    NaN    NaN    NaN
#> Linked pairs (augm. RL)   348.00 348.00 348.00 348.00 348.00 348.00 348.00
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.71   0.72   0.73   0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator   0.09   0.09   0.09   0.09    0    0    0    0    0
#> FDP synth data estimator     NaN    NaN    NaN    NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)   348.00 348.00 348.00 348.00  217  217  217  217  217
#> Linked pairs (RL)         443.00 443.00 443.00 443.00  276  276  276  276  276
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    217  217  217  217  217  217  217  217  217  217
#> Linked pairs (RL)          276  276  276  276  276  276  276  276  276  276
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator   NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
#> Linked pairs (augm. RL)    217  217  217  217  217  217  217  217  217  217
#> Linked pairs (RL)          276  276  276  276  276  276  276  276  276  276
#> fedmatch results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.09   0.09   0.09   0.09   0.09   0.09   0.09
#> Linked pairs (RL)         443.00 443.00 443.00 443.00 443.00 443.00 443.00
#>                             0.71   0.72   0.73   0.74 0.75 0.76 0.77 0.78 0.79
#> FDP model score estimator   0.09   0.09   0.09   0.09    0    0    0    0    0
#> Linked pairs (RL)         443.00 443.00 443.00 443.00  275  275  275  275  275
#>                           0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)          275  275  275  275  275  275  275  275  275  275
#>                           0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)          275  275  275  275  275  275  275  275  275  275

df_results[6,"fedmatch  "] = diag_fedmatch$FDP_measures$curve$RL_FDP_score[1]
df_results[7,"fedmatch  "] = diag_fedmatch$FDP_measures$curve$augmRL_FDP_score[1]
df_results[8,"fedmatch  "] = diag_fedmatch$FDP_measures$curve$augmRL_FDP_synth[1]

smd_idx = grep("smd", names(diag_fedmatch$discrepancy_measures))
iou_idx = grep("iou", names(diag_fedmatch$discrepancy_measures))
mmd_idx = grep("mmd", names(diag_fedmatch$discrepancy_measures))

df_results[9, "fedmatch  "] = max(diag_fedmatch$discrepancy_measures[1,smd_idx])
df_results[10,"fedmatch  "] = min(diag_fedmatch$discrepancy_measures[1,iou_idx])
df_results[11,"fedmatch  "] = max(diag_fedmatch$discrepancy_measures[1,mmd_idx])
```

``` r

options(digits = 2)

cat("Record Linkage ")
#> Record Linkage
df_results[1:5,]
#>                  Naive      FlexRL     BRL        multilink  fastLink  
#> TP                   258.00     270.00     235.00     243.00     229.00
#> FP                    94.00      46.00      32.00      41.00      46.00
#> FN                   142.00     130.00     165.00     157.00     171.00
#> sensitivity            0.65       0.68       0.59       0.61       0.57
#> FDP                    0.27       0.15       0.12       0.14       0.17
#>                  reclin2    diyar      fedmatch  
#> TP                   245.00     82.000     268.00
#> FP                    67.00      2.000     175.00
#> FN                   155.00    318.000     132.00
#> sensitivity            0.61      0.205       0.67
#> FDP                    0.21      0.024       0.40
```

``` r

options(digits = 2)

cat("\nFDP estimation ")
#> 
#> FDP estimation
df_results[6:8,]
#>                  Naive      FlexRL     BRL        multilink  fastLink  
#> hat FDP score            NA       0.10       0.14        NaN       0.23
#> hat FDP score*           NA       0.13       0.15        NaN       0.24
#> hat FDP synth*           NA       0.24       0.28       0.17       0.20
#>                  reclin2    diyar      fedmatch  
#> hat FDP score          0.30        NaN      0.095
#> hat FDP score*         0.31        NaN      0.094
#> hat FDP synth*         0.30          0        NaN
```

``` r

options(digits = 2)

cat("\nLinked divergence ")
#> 
#> Linked divergence
df_results[9:11,]
#>                  Naive      FlexRL     BRL        multilink  fastLink  
#> max SMD                  NA     0.0932     0.1515     0.0883      0.161
#> min support IoU          NA     0.8750     0.8750     0.8750      0.875
#> max MMD                  NA     0.0058     0.0082     0.0054      0.005
#>                  reclin2    diyar      fedmatch  
#> max SMD              0.1948      0.349     0.1747
#> min support IoU      0.8750      0.875     0.8750
#> max MMD              0.0066      0.082     0.0055
```
