# Statistics

SOPRANO produces two sets of statistical constraints on the dN/dS values
computed by the pipeline.

***

## Analytic estimates

For each pipeline run, there is an analytic estimate of confidence intervals
(68 and 98 percentiles) and p-value. This uses Katz analytic method of
estimating the confidence intervals, which can then be used to estimate
a pvalue. (For example,
see [this reference](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Ratio_of_Two_Proportions.pdf).)

***

## Numerical estimates

SOPRANO computes dN/dS based on the inclusion/exclusion of genomic regions
defined by an input bed protein file.

If we instead randomize these regions before executing the pipeline,
the downstream results can be interpreted as a null hypothesis
simulation.

By randomizing and running the SOPRANO pipeline we are effectively drawing
a "sample" from the unknown distribution of null hypothesis dN/dS values.
By building an ensemble of many samples, it is therefore possible to build a
kernel density estimate (KDE) from the empirical distribution.

Once the kernel has been estimated, it is possible to assess the statistical
significance of the non-randomized dN/dS result. A p-value can be computed
by direct integration from the asymptotic tail of the kernel up until the
observed non-randomized dN/dS value(s).

### The Gaussian kernel

The kernel is derived by applying a Gaussian kernel to the distribution of
dN/dS from null hypothesis simulations. The kernel itself is computed with
`scikit-learn`, for which we apply a grid search to find an optimal bandwidth
for the window function.

By default, the grid bandwidth search is performed in 250 logarithmic
spaces between -3 and 10. This is fairly sensible given the
order-of-magnitudes one might expect for dN/dS in ON and OFF target regions.

These values can be adjusted at runtime via setting the environment variables:

- `SOPRANO_CV_LOG_MIN` - Base10 log value for minimum bandwidth.
- `SOPRANO_CV_LOG_MAX` - Base10 log value for maximum bandwidth.
- `SOPRANO_CV_LOG_STEPS` - Number of linear spaces in log space.

### Integration parameters

Once the kernel has been determined, one needs to determine the lower and
upper bounds for the integration from which, p-values are computed.

Based on the spread of dN/dS values from the null distribution,
a step distance is defined as some percentage of that spread.
To determine the integration bounds, the kernel is estimated iteratively by
"stepping" away from distribution maximum/minimum, until the absolute and
relative difference of the kernel evaluated at that step is below a critical
threshold. These values define the lower and upper bounds for the integration.
The integration itself is performed with the SciPy `quad` integrator.

By default, SOPRANO walks away from the null distribution in 1% steps;
until the absolute and relative difference between steps is < 1e-8.
This will continue for no more than 1000 steps.

These values can be adjusted at runtime via setting the environment variables:

- `SOPRANO_KDE_ABS_TOL` - The absolute tolerance for convergence.
- `SOPRANO_KDE_REL_TOL` - The relative tolerance for convergence.
- `SOPRANO_KDE_STEP_PCT` - Percentage of distributions spread to use as a step.
- `SOPRANO_KDE_MAX_ITER` - Maximum number of iterations before termination.
