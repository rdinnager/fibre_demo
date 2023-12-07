library(tidyverse)
library(phyf)
library(fibre)
library(glmnet)
library(survival)

clads_pf <- read_rds("data/clads_pf.rds")

#### check equivalence ####
clads_tree <- pf_as_phylo(clads_pf)
avonet_tree <- pf_as_phylo(avonet)

ape::all.equal.phylo(clads_tree, avonet_tree)

bird_pf <- clads_pf |>
  left_join(avonet |>
              select(label,
                     Beak.Length_Culmen:Mass))

bird_pf <- bird_pf |>
  mutate(len = pf_end_features(phlo),
         end_time = pf_flow_sum(phlo),
         start_time = end_time - len,
         cens = as.numeric(!is_tip))

tree_surv <- Surv(bird_pf$start_time / 100, bird_pf$end_time / 100,
                  bird_pf$cens)

tree_x <- pf_as_sparse(bird_pf)

fit <- glmnet(tree_x, tree_surv, family = "cox",
              trace.it = TRUE)

autoplot(bird_pf, div_rate, edge_traits = TRUE) +
  scale_color_viridis_c(trans = "log1p", option = "inferno")

set.seed(2)
nobs <- 1000; nvars <- 15
xvec <- rnorm(nobs * nvars)
xvec[sample.int(nobs * nvars, size = 0.4 * nobs * nvars)] <- 0
x <- matrix(xvec, nrow = nobs)  # non-sparse x
x_sparse <- Matrix::Matrix(xvec, nrow = nobs, sparse = TRUE)  # sparse x

# create start-stop data response
beta <- rnorm(5)*5
fx <- x_sparse[, 1:5] %*% beta / 3
ty <- rexp(nobs, drop(exp(fx)))
tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
starty <- runif(nobs)
yss <- Surv(starty, starty + ty, tcens)

# fit regularized Cox model with start-stop data
fit <- glmnet(x, yss, family = "cox")

glmnet_fit <- glmnet(x, yss, family = "cox", lambda = 0)
coxph_fit <- coxph(yss ~ x)
plot(coef(glmnet_fit), coef(coxph_fit))
abline(0, 1)
