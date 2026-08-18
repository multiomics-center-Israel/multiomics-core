test_that("keep_kegg_pathways drops excluded classes and keeps the rest", {
  cache_dir <- withr::local_tempdir()
  # Synthetic stand-in for the BRITE hierarchy, so the test never hits KEGG.
  saveRDS(
    data.frame(
      pathway_id   = c("00010", "04260", "05416", "01523", "04390"),
      category     = c("Metabolism", "Organismal Systems", "Human Diseases",
                       "Human Diseases", "Environmental Information Processing"),
      subcategory  = c("Carbohydrate metabolism", "Circulatory system",
                       "Cardiovascular disease", "Drug resistance",
                       "Signal transduction"),
      pathway_name = c("Glycolysis", "Cardiac muscle contraction",
                       "Viral myocarditis", "Antifolate resistance",
                       "Hippo signaling pathway"),
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "kegg_pathway_categories.rds")
  )

  ids <- c("map00010", "map04260", "map05416", "map01523", "ko04390")
  keep <- keep_kegg_pathways(
    ids,
    exclude   = c("Human Diseases", "Organismal Systems"),
    cache_dir = cache_dir
  )

  expect_equal(keep, c(TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("keep_kegg_pathways strips organism prefixes before matching", {
  cache_dir <- withr::local_tempdir()
  saveRDS(
    data.frame(
      pathway_id = "04260", category = "Organismal Systems",
      subcategory = "Circulatory system", pathway_name = "Cardiac muscle contraction",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "kegg_pathway_categories.rds")
  )

  # The same map reached through three id spaces must classify identically.
  keep <- keep_kegg_pathways(
    c("map04260", "ko04260", "hsa04260"),
    exclude = "Organismal Systems", cache_dir = cache_dir
  )
  expect_equal(keep, c(FALSE, FALSE, FALSE))
})

test_that("keep_kegg_pathways matches on subcategory as well as category", {
  cache_dir <- withr::local_tempdir()
  saveRDS(
    data.frame(
      pathway_id = c("04260", "04620"), category = c("Organismal Systems", "Organismal Systems"),
      subcategory = c("Circulatory system", "Immune system"),
      pathway_name = c("Cardiac muscle contraction", "Toll and Imd signaling pathway"),
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "kegg_pathway_categories.rds")
  )

  # Excluding one subcategory must not take its siblings with it: an insect
  # project keeps immune maps while dropping vertebrate circulatory ones.
  keep <- keep_kegg_pathways(
    c("map04260", "map04620"),
    exclude = "Circulatory system", cache_dir = cache_dir
  )
  expect_equal(keep, c(FALSE, TRUE))
})

test_that("keep_kegg_pathways never silently drops unclassified pathways", {
  cache_dir <- withr::local_tempdir()
  saveRDS(
    data.frame(
      pathway_id = "00010", category = "Metabolism",
      subcategory = "Carbohydrate metabolism", pathway_name = "Glycolysis",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "kegg_pathway_categories.rds")
  )

  keep <- keep_kegg_pathways(
    c("map00010", "map99999"),
    exclude = "Human Diseases", cache_dir = cache_dir
  )
  expect_equal(keep, c(TRUE, TRUE))
})

test_that("keep_kegg_pathways is a no-op without an exclusion list", {
  ids <- c("map00010", "map05416")
  expect_equal(keep_kegg_pathways(ids, exclude = NULL), c(TRUE, TRUE))
  expect_equal(keep_kegg_pathways(ids, exclude = character(0)), c(TRUE, TRUE))
})

test_that("keep_kegg_pathways keeps everything when the classification is unavailable", {
  # Offline or first-run-without-cache must degrade to reporting more, not less.
  # The stub is installed by assignment rather than with_mocked_bindings(): these
  # functions are sourced into the global environment, not loaded from a package,
  # and the mocking helper requires a package namespace to rebind in.
  env <- environment(keep_kegg_pathways)
  original <- get("kegg_pathway_categories", envir = env)
  assign("kegg_pathway_categories", function(...) NULL, envir = env)
  on.exit(assign("kegg_pathway_categories", original, envir = env), add = TRUE)

  keep <- keep_kegg_pathways(c("map05416", "map04260"),
                             exclude = "Human Diseases",
                             cache_dir = withr::local_tempdir())
  expect_equal(keep, c(TRUE, TRUE))
})
