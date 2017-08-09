# lapply(dir("R", full.names = T), source)
# to use the package, you have to run "Build - Load All" or "Ctr+Shift+L" in RStudio for now
set.seed(1)
# to later calculate error we use the cosine similarity
cosine <- function(x, y){
  sum(x * y) / (sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2)))
}
library(readr)
library(dplyr)
# This program constructs a 3-layer random associator using the 'network',
# 'layer', and 'unit' objects. It is a test of these
# clases, and a tutorial in how to construct networks with them.

# 1) First, specify the dimensions and connectivity of the layers At this
# point, layers are either fully connected or unconnected

# 1.1) Set the dimensions of the layers specifying 3 layers and their dimensions

# 1.1) Set the dimensions of the layers specifying 3 layers and their dimensions
cxn_tbl <- read_delim("cxn_table.csv", ";", col_types = list("i", "c", "c", "n", "l"))

lays_tbl <- read_delim("dim_table.csv", ";",
                       col_types = list("c", "i", "i", "n", "c"))

net <- network$new(lays_tbl, cxn_tbl)

create_rnd_input <- function(lays_tbl, n_inputs, prop_active_units){
  ext_inputs_tbl <- plyr::ddply(lays_tbl, 1, function(x)
    expand.grid(
      input = 1:n_inputs,
      layer = x$layer,
      type = x$type,
      unit_id = 1:(x$n_rows * x$n_cols)
    ))
  ext_inputs_tbl <- dplyr::mutate(ext_inputs_tbl,
                                  activation = ifelse(type == "hidden", NA,
                                                      sample(c(0.01, 0.96),
                                                             n(),
                                                             prop_active_units,
                                                             replace = T)))
  ext_inputs_tbl <- arrange(ext_inputs_tbl, c(input))
}

ext_inputs <- create_rnd_input(lays_tbl, 15, c(0.7, 0.3))
ext_inputs <- ext_input[ext_input == 1,]
## 3) Create the network using the constructor

self <- net
nxx1_df <- create_nxx1()
# go through all inputs
net$cycle(ext_input[ext_input == 1,], clamp_inp = 1)

## 4) Let's create some inputs
n_inputs <- 5  # number of input-output patterns to associate
patterns <- matrix(vector(mode = "list", length = n_inputs * 2), ncol = 2)
# is the i-th input pattern, and patterns[i,2] is the i-th output pattern.
# This will assume that layers 1 and 3 are input and output respectively.
# Patterns will be binary.

prop <- 0.3  # the proportion of active units in the patterns
for (i in seq(n_inputs)) {
    patterns[i, 1][[1]] <- sample(c(0, 1),
                                  n_units[1],
                                  prob = c(1 - prop, prop),
                                  replace = T) # either 0 or 1
    patterns[i, 2][[1]] <- sample(c(0, 1),
                                  n_units[3],
                                  prob = c(1 - prop, prop), replace = T)
    # in orig sim input act are 0.01 or 0.96, but if a pattern contains only
    # zeros (0.01) the output is also zero, so nothing will be learnt
    # therefore act is either 0.25 or 0.95
    patterns[i, 1][[1]] <- 0.25 + 0.7 * patterns[i, 1][[1]]
    patterns[i, 2][[1]] <- 0.25 + 0.7 * patterns[i, 2][[1]]
}

## 5) Train the network

# Specify parameters for training
n_epochs <- 10  # number of epochs. All input patterns are presented in one.
n_trials <- n_inputs  # number of trials. One input pattern per trial.
n_minus <- 50  # number of minus cycles per trial.
n_plus <- 25  # number of plus cycles per trial.
# learning rate schedule
lrate_sched <- seq(from = 0.8, to = 0.2, length.out = n_epochs)

# cosine error for each pattern
errors <- matrix(rep(0, n_epochs * n_trials), ncol = n_trials)

for (epoch in seq(n_epochs)) {
    order <- sample(seq(n_trials))  # order of presentation of inputs this epoch
    net$lrate <- lrate_sched[epoch]  # learning rate for this epoch
    for (trial in seq(n_trials)) {
        net$reset()  # randomize the acts for all units
        pat <- order[trial]  # input to be presented this trial
        #++++++ MINUS PHASE +++++++
        inputs <- list(unlist(patterns[pat, 1]), c(), c())

        for (minus in seq(n_minus)){
            # minus cycles: layer 1 is clamped
            net$cycle(inputs, clamp_inp = 1)
        }
        outs <- net$lays[[3]]$get_acts() # saving the output for testing

        #+++++++ PLUS PHASE +++++++
        inputs <- list(unlist(patterns[pat, 1]), c(), unlist(patterns[pat, 2]))
        for (plus in seq(n_plus)){
            # plus cycles: layers 1 and 3 are clamped
            net$cycle(inputs, clamp_inp = 1)
        }

        net$updt_pct_act() # update averages used for net input scaling

        #+++++++ LEARNING +++++++
        net <- net$chg_wt()  # updates the avg_l vars and applies XCAL learning

        #+++++++ ERRORS +++++++
        # Only the cosine error is used here
        errors[epoch, pat] <- 1 - cosine(unlist(patterns[pat, 2]), outs)
    }
    cat("\nepoch ", epoch, " finished")
}

## 6) A plot of the error by epoch figure
plot(rowMeans(errors))
