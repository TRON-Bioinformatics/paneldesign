
<!-- README.md is generated from README.Rmd. Please edit that file -->

# paneldesign

<!-- badges: start -->

[![R-CMD-check](https://github.com/TRON-Private/paneldesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TRON-Private/paneldesign/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This R package provides functions to design and evaluate seqsuencing
panels. A sequencing panel is a collection of genomic regions selected
for targeted sequencing, such as cancer hotspot mutations. Given a
cohort of individuals with associated genomic features (e.g., patients
with somatic mutations), the typical objective is to design a panel that
captures the largest number of individuals while using the smallest
possible subset of features.

## Installation

### Install Bioconductor dependencies

This package depended on the Biocondoctor packages
[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
and
[GenomeInfoDb](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html).
Install it by

``` r
install.packages("BiocManager") # if needed
BiocManager::install(c("GenomicRanges", "GenomeInfoDb"))
```

### Install this package from GitHub

This R package is not yet on [CRAN](https://CRAN.R-project.org) or
[Biodonducotr](https://www.bioconductor.org/). Therefore, you have to
install it form this GitHub repository.

``` r
# install.packages("remotes")
remotes::install_github("tron-bioinformatics/paneldesign")
```

## Usage example

This is a basic example to show you how to use the package:

``` r
library(GenomicRanges)
library(paneldesign)
```

Given a data set of mutations in patients, the goal is to select
mutations that are present in a maximum number of patients. We use a toy
data set consisting of patients `p1`, `p2`, and `p3` of which each has 3
mutations.

``` r
mut_toy
#> # A tibble: 9 × 6
#>   patient_id mut_id chr   start   end gene 
#>   <chr>      <chr>  <chr> <dbl> <dbl> <chr>
#> 1 p1         m01    1      1000  1000 g1   
#> 2 p1         m02    1      2000  2000 g1   
#> 3 p1         m03    2      3000  3000 g2   
#> 4 p2         m01    1      1000  1000 g1   
#> 5 p2         m04    2      4000  4000 g2   
#> 6 p2         m05    3      5000  5000 g3   
#> 7 p3         m02    1      2000  2000 g1   
#> 8 p3         m04    2      4000  4000 g2   
#> 9 p3         m06    2      6000  6000 g2
```

### Select mutations

Here, we use a greedy algorithm to select a minimal set of mutations
that covers all patients.

``` r
select_greedy(mut_toy)
#> 100% covered by 2 sets.
#> # A tibble: 2 × 8
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
#> # A tibble: 6 × 10
#>   reg_id reg_chr reg_start reg_end mut_id patient_id chr   start   end gene 
#>   <chr>  <chr>       <int>   <int> <chr>  <chr>      <chr> <dbl> <dbl> <chr>
#> 1 r1     1            1000    1000 m01    p1         1      1000  1000 g1   
#> 2 r1     1            1000    1000 m01    p2         1      1000  1000 g1   
#> 3 r1     1            1000    1000 m01    p1         1      1000  1000 g1   
#> 4 r1     1            1000    1000 m01    p2         1      1000  1000 g1   
#> 5 r3     2            3000    3000 m03    p1         2      3000  3000 g2   
#> 6 <NA>   <NA>           NA      NA <NA>   p3         <NA>     NA    NA <NA>
```

This results in the above-shown data set. The patient `p3` is not
*covered* because there are no mutations in `p3` that overlap with the
input regions. This results in ‘NA’ values for some region or
mutation-specific columns.

Next, we evaluate each region from the panel for the patient coverage.

``` r
eval_regions(reg_to_patient)
#> # A tibble: 2 × 15
#>   reg_id reg_chr reg_start reg_end n_region reg_size reg_size_cum n_patient
#>   <chr>  <chr>       <int>   <int>    <int>    <dbl>        <dbl>     <int>
#> 1 r1     1            1000    1000        1        1            1         2
#> 2 r3     2            3000    3000        2        1            2         1
#> # ℹ 7 more variables: n_patients_cum <dbl>, total_patients <int>,
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
#> # A tibble: 1 × 9
#>   n_reg size_total total_patients n_patient percent_patient n_mut
#>   <int>      <dbl>          <int>     <int>           <dbl> <int>
#> 1     2          2              3         2            66.7     2
#> # ℹ 3 more variables: mut_per_patient_median <int>, mut_per_patient_mean <dbl>,
#> #   mut_per_patient_df <list>
```

This data set contains a nested data set with information on how many
mutations each patient has represented on the panel.

``` r
eval$mut_per_patient_df[[1]]
#> # A tibble: 3 × 2
#>   patient_id mut_per_patient_n
#>   <chr>                  <int>
#> 1 p1                         2
#> 2 p2                         1
#> 3 p3                         0
```

### Structural variants

Structural variants (SVs) are represented as pairs of breakpoint
coordinates. And to be covered by a panel region either a single
(`sv_mode = "single"`) or both (`sv_mode = "both"`) breakpoints must be
included in a region of the panel.

An toy example SV dataset might look like this

``` r
sv_toy
#> # A tibble: 2 × 6
#>   patient_id mut_id bp1_chr bp1_pos bp2_chr bp2_pos
#>   <chr>      <chr>  <chr>     <dbl> <chr>     <dbl>
#> 1 p4         sv01   1          1000 1          1500
#> 2 p5         sv02   1          2000 2          2000
```

We can evaluate a dataset with mutations and SVs by passing it as
additional input to the `panel_to_patient()` function

``` r
reg_to_patient <- panel_to_patient(gr_toy, mut_toy, sv_toy, sv_mode = "single")
```

Now the `reg_to_patient` table includes additonal colums for the SV
breakpoint coordinates and includes here an association of SV `sv01`
from patient `p4` with region `r1`. In contrast patient `p5`is not
covered by the panel.

``` r
reg_to_patient
#> # A tibble: 9 × 14
#>   reg_id reg_chr reg_start reg_end mut_id patient_id chr   start   end gene 
#>   <chr>  <chr>       <int>   <int> <chr>  <chr>      <chr> <dbl> <dbl> <chr>
#> 1 r1     1            1000    1000 m01    p1         1      1000  1000 g1   
#> 2 r1     1            1000    1000 m01    p2         1      1000  1000 g1   
#> 3 r1     1            1000    1000 m01    p1         1      1000  1000 g1   
#> 4 r1     1            1000    1000 m01    p2         1      1000  1000 g1   
#> 5 r3     2            3000    3000 m03    p1         2      3000  3000 g2   
#> 6 r1     1            1000    1000 sv01   p4         <NA>     NA    NA <NA> 
#> 7 r3     2            3000    3000 <NA>   <NA>       <NA>     NA    NA <NA> 
#> 8 <NA>   <NA>           NA      NA <NA>   p3         <NA>     NA    NA <NA> 
#> 9 <NA>   <NA>           NA      NA <NA>   p5         <NA>     NA    NA <NA> 
#> # ℹ 4 more variables: bp1_chr <chr>, bp1_pos <dbl>, bp2_chr <chr>,
#> #   bp2_pos <dbl>
```
