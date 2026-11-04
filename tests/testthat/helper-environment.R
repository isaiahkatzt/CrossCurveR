# testthat sets its message language locally; LC_ALL can override that on macOS.
# Restore the caller's locale override after the entire test run.
withr::local_envvar(c(LC_ALL = NA), .local_envir = testthat::teardown_env())
