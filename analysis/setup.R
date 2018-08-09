# startup.func <- function()
# {
#         # startup.func:
#         # Sources the GAM functions into R.
#         # Usage:
#         #   source("startup.func")
#         #   startup.func()
#         # This function only needs to be run once.
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/indsp.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/bootstrap.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/outer.boot.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/diffop.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/indmat.derivmat.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/index.ci.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/deriv.ci.func")
#         source("https://www.stat.auckland.ac.nz/~fewster/gams/R/sp.plot")
# }

# startup.func()


tt <- readLines("https://www.stat.auckland.ac.nz/~fewster/gams/R/cb")
cb <- read.table(textConnection(tt), skip = 1)
colnames(cb) <- c("site", "year", "count")
usethis::use_data(cb, overwrite = FALSE)
