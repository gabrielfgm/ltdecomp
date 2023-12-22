
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ltdecomp

<!-- badges: start -->
<!-- badges: end -->

The goal of ltdecomp is to build simple life-tables from mortality data
and then perform life-table decompositions to isolate the contribution
of specific causes of death to life-expectancy. The package is designed
for private use on historical mortality data, but should work for any
collection of information on populations and associated numbers of
deaths.

## Installation

You can install the development version of ltdecomp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabrielfgm/ltdecomp")
```

## Example

This is a basic example of constructing a lifetable:

``` r
library(ltdecomp)
## basic example code

# Age groups
Age <-  c("0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
          "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
          "70-74", "75-79", "80-84", "85-89", "90-94", "95+") 

# Death counts
Dx <- c(110934, 17693, 9168, 7378, 12192, 13356, 14223, 19210, 29136, 42961, 64310, 
        90630, 116759, 153504, 196682, 223796, 220065, 185305, 120414, 50299, 13887)

# mid year population
Px <- c(4126560, 16195304, 18659141, 16815965, 13287434, 10803165, 10870165,
        11951709, 12508316, 11567216, 10528878, 9696502, 8595947, 7111897, 6186763,
        4661136, 2977347, 1518206, 648581, 170653, 44551)

lifetable(age = Age, Dx = Dx, Px = Px)
#>      age           mx       ax          qx          lx          dx        Lx
#> 1      0 0.0268829243 0.350000 0.026421242 1.000000000 0.026421242 0.9828262
#> 2    1-4 0.0010924772 2.000000 0.004360381 0.973578758 0.004245175 3.8858247
#> 3    5-9 0.0004913409 2.500000 0.002453691 0.969333583 0.002378445 4.8407218
#> 4  10-14 0.0004387497 2.500000 0.002191345 0.966955138 0.002118932 4.8294784
#> 5  15-19 0.0009175586 2.500000 0.004577293 0.964836206 0.004416338 4.8131402
#> 6  20-24 0.0012363044 2.500000 0.006162475 0.960419868 0.005918563 4.7873029
#> 7  25-29 0.0013084438 2.500000 0.006520888 0.954501304 0.006224197 4.7569460
#> 8  30-34 0.0016073015 2.500000 0.008004344 0.948277108 0.007590336 4.7224097
#> 9  35-39 0.0023293303 2.500000 0.011579222 0.940686771 0.010892421 4.6762028
#> 10 40-44 0.0037140311 2.500000 0.018399316 0.929794350 0.017107580 4.6062028
#> 11 45-49 0.0061079633 2.500000 0.030080490 0.912686770 0.027454065 4.4947987
#> 12 50-54 0.0093466696 2.500000 0.045666279 0.885232704 0.040425283 4.3251003
#> 13 55-59 0.0135830293 2.500000 0.065684655 0.844807421 0.055490884 4.0853099
#> 14 60-64 0.0215841146 2.500000 0.102395294 0.789316537 0.080822299 3.7445269
#> 15 65-69 0.0317907765 2.500000 0.147250837 0.708494238 0.104326369 3.2816553
#> 16 70-74 0.0480131882 2.500000 0.214338281 0.604167869 0.129496302 2.6970986
#> 17 75-79 0.0739131180 2.500000 0.311926871 0.474671567 0.148062817 2.0032008
#> 18 80-84 0.1220552415 2.500000 0.467595119 0.326608750 0.152720658 1.2512421
#> 19 85-89 0.1856576125 2.500000 0.634014170 0.173888093 0.110247515 0.5938217
#> 20 90-94 0.2947443057 2.500000 0.848497219 0.063640578 0.053998853 0.1832058
#> 21   95+ 0.3117101749 3.208108 1.000000000 0.009641725 0.009641725 0.0309317
#>            Tx        ex
#> 1  69.5919472 69.591947
#> 2  68.6091210 70.471054
#> 3  64.7232963 66.770921
#> 4  59.8825745 61.929010
#> 5  55.0530962 57.059526
#> 6  50.2399560 52.310409
#> 7  45.4526531 47.619268
#> 8  40.6957070 42.915416
#> 9  35.9732973 38.241526
#> 10 31.2970945 33.660233
#> 11 26.6908917 29.244307
#> 12 22.1960930 25.073738
#> 13 17.8709927 21.153925
#> 14 13.7856828 17.465341
#> 15 10.0411559 14.172530
#> 16  6.7595006 11.188117
#> 17  4.0624020  8.558343
#> 18  2.0592012  6.304795
#> 19  0.8079591  4.646432
#> 20  0.2141375  3.364794
#> 21  0.0309317  3.208108
```