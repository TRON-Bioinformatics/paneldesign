
<!-- README.md is generated from README.Rmd. Please edit that file -->

# paneldesign

<!-- badges: start -->

<!-- badges: end -->

The goal of this R package is to provide functionalities to design and
evaluate panels.

<https://tron.pages.gitlab.rlp.net/paneldesign>

## Installation

### Install Bioconductor dependencies

This package dependes on the Biocondocr packages
[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
and
[GenomeInfoDb](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html).
Install it by

``` r
install.packages("BiocManager") # if needed
BiocManager::install(c("GenomicRanges", "GenomeInfoDb"))
```

### Install this package from this GitLab

This R packge is not yet on [CRAN](https://CRAN.R-project.org) or
[Biodonducotr](https://www.bioconductor.org/). Therfore, you have to
install it form this GitLap repository.

However, this repository is a private GitLab reposoritory and therfore
you have to create an personal access token (PAT) first. This is
described
[here](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html).
As Scope use `api` or `read_api`.

``` r
install.packages("remotes") # if needed

remotes::install_gitlab("tron/paneldesign", host = "gitlab.rlp.net", auth_token = "YOUR_ACCESS_TOKEN")
```

## Example

This is a basic example to show you how to use the package:

``` r
library(paneldesign)
```

Given a data set of mutations in patients, the goal is to select
mutations that are present in a maximum number of patients. We use a toy
data set consisting of patients `p1`, `p2`, and `p3` of which each has 3
mutations.

``` r
mut_toy
#> # A tibble: 9 x 5
#>   patient_id mut_id chr   start   end
#>   <chr>      <chr>  <chr> <dbl> <dbl>
#> 1 p1         m01    1      1000  1000
#> 2 p1         m02    1      2000  2000
#> 3 p1         m03    2      3000  3000
#> 4 p2         m01    1      1000  1000
#> 5 p2         m04    2      4000  4000
#> 6 p2         m05    3      5000  5000
#> 7 p3         m02    1      2000  2000
#> 8 p3         m04    2      4000  4000
#> 9 p3         m06    2      6000  6000
```

### Select mutations

Here, we use a greedy algorithm to select a minimal set of mutations
that covers all patients.

``` r
select_greedy(mut_toy)
#> 100% covered by 2 sets.
#> # A tibble: 2 x 8
#>   mut_id     n n_samples  rank order n_cum coverage coverage_cum
#>   <chr>  <int>     <int> <dbl> <int> <int>    <dbl>        <dbl>
#> 1 m01        2         3     1     1     2     66.7         66.7
#> 2 m02        1         3     2     2     3     33.3        100
```

This results in two mutations `m01` and `m02`.

### Evaluate given panel of regions

Here we use the following toy genomic region set.

``` r
gr_toy
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> GRanges object with 2 ranges and 0 metadata columns:
#>      seqnames    ranges strand
#>         <Rle> <IRanges>  <Rle>
#>   r1        1      1000      *
#>   r3        2      3000      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

First, we associate each region to patients by overlap with the mutation
data

``` r
reg_to_patient <- panel_to_patient(gr_toy, mut_toy)

reg_to_patient
#> # A tibble: 6 x 9
#>   reg_id reg_chr reg_start reg_end mut_id patient_id chr   start   end
#>   <chr>  <chr>       <int>   <int> <chr>  <chr>      <chr> <dbl> <dbl>
#> 1 r1     1            1000    1000 m01    p1         1      1000  1000
#> 2 r1     1            1000    1000 m01    p2         1      1000  1000
#> 3 r1     1            1000    1000 m01    p1         1      1000  1000
#> 4 r1     1            1000    1000 m01    p2         1      1000  1000
#> 5 r3     2            3000    3000 m03    p1         2      3000  3000
#> 6 <NA>   <NA>           NA      NA <NA>   p3         <NA>     NA    NA
```

This results in the above-shown data set. The patient `p3` is not
*covered* because there are no mutations in `p3` that overlap with the
input regions. This results in ‘NA’ values for some region or
mutation-specific columns.

Next, we evaluate each region from the panel for the patient coverage.

``` r
eval_regions(reg_to_patient)
#> # A tibble: 2 x 15
#>   reg_id reg_chr reg_start reg_end n_region reg_size reg_size_cum n_patient
#>   <chr>  <chr>       <int>   <int>    <int>    <dbl>        <dbl>     <int>
#> 1 r1     1            1000    1000        1        1            1         2
#> 2 r3     2            3000    3000        2        1            2         1
#> # … with 7 more variables: n_patients_cum <dbl>, total_patients <int>,
#> #   percent_patients <dbl>, percent_patients_cum <dbl>, n_patient_gain <dbl>,
#> #   n_mut <int>, n_mut_cum <dbl>
```

This assumes some order of the input regions by decreasing priority and
in this way outputs additionally cumulative values such as
`percent_patients_cum` which the cumulative percent of input patients
that are covered by the current or any previous region. This is useful
to decide where to *cut* the panel.

Lastly, we can get some overall performance metrics of the entire panel:

``` r
eval <- eval_panel(reg_to_patient)
eval
#> # A tibble: 1 x 9
#>   n_reg size_total total_patients n_patient percent_patient n_mut
#>   <int>      <dbl>          <int>     <int>           <dbl> <int>
#> 1     2          2              3         2            66.7     2
#> # … with 3 more variables: mut_per_patient_median <int>,
#> #   mut_per_patient_mean <dbl>, mut_per_patient_df <list>
```

This data set contains a nested data set with information on how many
mutations each patient has represented on the panel.

``` r
eval$mut_per_patient_df[[1]]
#> # A tibble: 3 x 2
#>   patient_id mut_per_patient_n
#> * <chr>                  <int>
#> 1 p1                         2
#> 2 p2                         1
#> 3 p3                         0
```
