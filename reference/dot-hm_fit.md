# Iterative fit of smoothed doubled angles at the observation points

Objective per observation i (Hansen & Mount 1990, eq. 5; Müller et al.
2003, eq. 3): M_i = q_i delta_i y_i y_i^T + lambda \* sum_j w~\_ij
delta_j f_j f_j^T. The identity part of y y^T does not change the
dominant eigenvector, so only the doubled-angle components (cos 2a, sin
2a) are carried.

## Usage

``` r
.hm_fit(nb, n, th2, q, lambda, n_iter, tol, robust_scale = 0)
```

## Arguments

- robust_scale:

  minimum robustness scale in chord units (2 \* sin(deg / 2)), see
  circular_smooth().
