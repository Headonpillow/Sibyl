utils::globalVariables(c(
  "Site", "data", "fit", "newdata", "preds", ".fitted", "Sample", "Species",
  "params", "term", "estimate", "ACE", "b", "sample_75", "density",
  "Threshold", "Avg_Distance", "Sample_ID", "Dim1", "Dim2", ".data",
  "sd", "Index", "Repeat_Amount", "after_stat", "xlim", "ylim",
  "b_est", "rmse", "max_species", "max_sample", "nrmse", "converged",
  "pass_local_range", "pass_global_range", "pass_nrmse", "fit_ok"
))

# Internal adapter: always return a plain samples-by-features numeric matrix.
# Converting an S4 otu_table directly with data.frame() can retain S4 columns.
sample_count_matrix <- function(input) {
  counts <- methods::as(phyloseq::otu_table(input), "matrix")
  if (phyloseq::taxa_are_rows(input)) counts <- t(counts)
  counts
}
