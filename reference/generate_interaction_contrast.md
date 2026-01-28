# Fast factorial contrast generators

Returns a matrix **N_cells x N_contrasts** - *each row is a design
cell*, columns are independent contrasts (difference-coded for the
factors you ask for, grand-mean for the rest). Suitable for
`tcrossprod(dm, C)` or `lm.fit(design, y)` followed by `%*% coef` in the
usual way.

## Usage

``` r
generate_interaction_contrast(des, factors)

generate_main_effect_contrast(des, factor)
```

## Arguments

- des:

  data.frame with one column per factor (must be `factor`)

- factors:

  character vector: which factor(s) get **difference coding**. -
  `generate_main_effect_contrast()` takes a **single** factor name.  
  - `generate_interaction_contrast()` takes \>= 2 for an interaction (or
  1 to reproduce a main-effect matrix).

- factor:

  Single factor name for the main effect.

## Value

numeric matrix **nrow = prod levels(f) , ncol = prod (Li - 1)** for the
chosen factors.

## Examples

``` r
des <- expand.grid(Time = factor(1:4),
                   Cond = factor(c("face","scene")))

# Main effect of Time (4-1 = 3 contrasts)
M <- generate_main_effect_contrast(des, "Time")

# Full TimexCond interaction ( (4-1)*(2-1) = 3 contrasts )
I <- generate_interaction_contrast(des, c("Time","Cond"))
dim(I)   # 8 rows (cells) x 3 columns (contrasts)
#> [1] 8 3
```
