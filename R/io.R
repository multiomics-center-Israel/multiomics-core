read_config <- function(path = "config/config.yml") {
  yaml::read_yaml(path)
}

save_table <- function(df, path) {
  readr::write_csv(df, path)
}
