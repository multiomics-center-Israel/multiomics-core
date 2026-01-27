# Run Pipeline with Fresh Code

cat("Running tar_make() with freshly loaded code...\n")
targets::tar_destroy()
targets::tar_make()
