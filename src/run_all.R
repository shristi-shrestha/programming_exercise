# run_all.R
# ---------------------------------------------------------------------------
# Run the whole exercise end to end. Execute from the programming_exercise/
# folder:   Rscript src/run_all.R
# Each script is self-contained (sources 00-setup.R itself) and writes its
# figures / tables to results/.
# ---------------------------------------------------------------------------

scripts <- c(
  "src/01-explore.R",       # Q1: filtering counts + gene-list overlap
  "src/02-pca.R",           # Q2: combined PCA + variance explained
  "src/03-proportions.R",   # Q3: subtype proportions by race
  "src/04-model-white.R",   # Q4-5: train on White, apply to Black
  "src/05-model-combined.R" # Q6: combined model with 20% held-out
)

for (s in scripts) {
  message("\n===================================================")
  message("Running ", s)
  message("===================================================")
  source(s, echo = FALSE)
}
message("\nAll scripts finished. See results/ for figures and tables.")
