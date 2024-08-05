# spatInfer

The goal is to run a regression with a spatial basis made up of the first $p$ principal components of a $k \times k$ tensor with standard errors based on $c$ large clusters. The spatial basis controls for trends and other large scale structure in the variables, and the clusters deal with autocorrelation in residuals.

``` r
library(spatInfer)

basis_reg(y ~x1 + x2, dat,
           clusters=c,
           splines=k,
           num_pcs=p
          )
```

Naturally, before we run this we need to come up with sensible values of the parameters `c`, `k` and `p`. This requires two simple commands.

The first step is to run the `optimal_basis` command

``` r
optimal_basis(y ~ x1 + x2, dat)
```

This chooses the combination of $k$ and $p$ that optimally predicts the outcome `y` based on a Bayes Information Criterion.

Knowing $k$ and $p$, the second step is to run the `placebo` command to choose the optimal number of `c` of large clusters by comparing the p-values of spatial noise placebos against that of the treatment `x1`.

``` r
placebo(y ~ x1 + x2, dat,
            splines=k,
            num_pcs=p
        )
```

For each number of clusters `c`, `placebo` gives the proportion of simulations with p-values below 0.05. The optimal cluster number is where this proportion is in the range 0.05 to 0.07. The command also gives the placebo significance level of the treatment variable: how often a placebo p-value is below that of the real one.

Knowing spatial basis values $k$ and $p$ and cluster number $c$ we can now run `basis_reg`.

Besides `placebo`, `spatInfer` gives a second diagnostic, a synthetic outcome test, that tests the hypothesis that the outcome `y` is spatial noise, and thus independent of the treatment `x1`.

``` r
synth(y ~ x1 + x2, dat,
            splines=k,
            num_pcs=p
)
```
