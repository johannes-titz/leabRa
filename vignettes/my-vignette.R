## ---- message = F--------------------------------------------------------
library(leabRa)

## ------------------------------------------------------------------------
set.seed(07221904)

## ------------------------------------------------------------------------
dim_lays <- list(c(2, 5), c(2, 10), c(2, 5))

## ------------------------------------------------------------------------
connections <- matrix(c(0, 0, 0,
                        1, 0, 0.2,
                        0, 1, 0), nrow = 3, byrow = T)

## ------------------------------------------------------------------------
net <- network$new(dim_lays, connections)

## ------------------------------------------------------------------------
inputs_plus <- net$create_inputs(which_layers = c(1, 3),
                                 n_inputs = 15,
                                 prop_active = .3)

## ------------------------------------------------------------------------
inputs_minus <- lapply(inputs_plus, function(x) {x[3] <- list(NULL); return(x)})

## ------------------------------------------------------------------------
n_epochs <- 2
outs <- lapply(seq(n_epochs), function(x) net$learn_error_driven(inputs_minus,
                                                                 inputs_plus))#,
                                                                 #lrate = 0.5))
w0 <- net$get_weights()
net$set_weights(w0)

## ------------------------------------------------------------------------
mad <- net$mad_per_epoch(outs, inputs_plus, 3)

## ---- fig.height=4, fig.show='hold', fig.width=6-------------------------
plot(mad, axes = F, pch = 16, family = "serif", type = "b",
     xlab = "epoch [#]",
     ylab = "mean absolute deviation [activation]",
     ylim = c(min(mad), max(mad)))
axis(1, at = seq(length(mad)), tick = T, family = "serif")
axis(2, at = seq(0, 1, 0.05), labels = seq(0, 1, 0.05), tick = T,
     family = "serif", las = 2)

## ------------------------------------------------------------------------
w_init_fun = function(x) 0.3 + 0.4 * runif(x)

## ------------------------------------------------------------------------
net <- network$new(dim_lays, connections, w_init_fun = function(x) rnorm(x, mean = 0.6, sd = 0.1))

## ------------------------------------------------------------------------
set.seed(22071904)

## ------------------------------------------------------------------------
animals

## ------------------------------------------------------------------------
inputs <- plyr::alply(animals, 1)

## ------------------------------------------------------------------------
inputs <- lapply(inputs, function(x) list(x, NULL))

## ------------------------------------------------------------------------
dim_lays <- list(c(6, 1), c(3, 1))
connections <- matrix(c(0, 0,
                        1, 0), nrow = 2, byrow = T)

## ------------------------------------------------------------------------
run_sim <- function(dim_lays, connections, inputs){
  net <- network$new(dim_lays, connections)
  net$learn_self_organized(inputs, lrate = 0.8)
  net$test_inputs(inputs)
}

## ---- message=FALSE------------------------------------------------------
n_runs <- 10
outs <- lapply(seq(n_runs), function(x) run_sim(dim_lays, connections, inputs))

## ------------------------------------------------------------------------
outs_layer_two <- lapply(outs, function(x) lapply(x, function(y) y[[2]]))
outs_layer_two <- lapply(outs_layer_two, function(x) do.call(rbind, x))
outs_layer_two <- lapply(outs_layer_two, round, 2)
round(outs_layer_two[[1]], 2)

## ------------------------------------------------------------------------
set_won <- function(matrix){
  t(apply(matrix, 1, function(x) {x[-which.max(x)] <- 0; x}))
}

set_won2 <- function(matrix){
  t(apply(matrix, 1, function(x) {x[which.max(x)] <- 1; x}))
}

outs_layer_two <- lapply(outs_layer_two, set_won)
outs_layer_two <- lapply(outs_layer_two, set_won2)

## ------------------------------------------------------------------------
dists <- lapply(outs_layer_two, dist)
dists <- lapply(dists, function(x) {x[x != 0] <- 1; x})
dists[[1]]

## ------------------------------------------------------------------------
dists_mtrx <- lapply(dists, as.matrix)
mean_dists <- Reduce("+", dists_mtrx) / length(dists)
mean_dists

## ------------------------------------------------------------------------
colnames(mean_dists) <- rownames(animals)
rownames(mean_dists) <- rownames(animals)

## ---- fig.height=5, fig.width=5------------------------------------------
plot(hclust(as.dist(mean_dists)))

