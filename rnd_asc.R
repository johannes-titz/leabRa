# lapply(dir("R", full.names = T), source)
# to use the package, you have to run "Build - Load All" or "Ctr+Shift+L" in RStudio for now
set.seed(1)
# to later calculate error we use the cosine similarity
cosine <- function(x, y){
  sum(x * y) / (sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2)))
}

# This program constructs a 3-layer random associator using the 'network',
# 'layer', and 'unit' objects. It is a test of these
# clases, and a tutorial in how to construct networks with them.

# 1) First, specify the dimensions and connectivity of the layers At this
# point, layers are either fully connected or unconnected

# 1.1) Set the dimensions of the layers specifying 3 layers and their dimensions
dim_lays <- list(c(1, 5), c(2, 5), c(1, 5))

# 1.2) Specify connectivity between layers
cxn <- matrix(c(0, 0, 0,
                1, 0, 0.2,
                0, 1, 0), nrow = 3, byrow = T)

# cxn(i,j) = c means that layer i receives cxn from j, and that
# they have a relative strength c. In this case, feedback cxn are 5
# times weaker. The relative weight scale of layer i comes from the non-zero
# entries of row i.  The network constructor will normalize this matrix so that
# if there are non-zero entries in a row, they add to 1.

net <- network$new(dim_lays, cxn)

create_all_inputs <- function(dim_lays, which_layers, n_inputs, prop = 0.3){
  lapply(1:n_inputs, function(x) create_one_input(dim_lays, which_layers, prop))
}

create_one_input <- function(dim_lays, which_layers, prop = 0.3){
  prob <- c(prop, 1 - prop)
  n_units <- lapply(dim_lays, function(x) x[1] * x[2])
  inputs <- mapply(create_one_input_for_one_layer,
                n_units,
                inv_which(which_layers, length(dim_lays)),
                MoreArgs = list(prob = prob))
  inputs
}

inv_which <- function(indices, totlength) is.element(seq_len(totlength), indices)

create_one_input_for_one_layer <- function(n, has_layer_input, prob){
  result <- ifelse(has_layer_input,
         return(sample(c(.96, 0.01), n, replace = T, prob = prob)),
         return(NULL))
  return(result)
}

patterns <- create_all_inputs(dim_lays, c(1, 3), 15)

# Train the network

# Specify parameters for training
n_epochs <- 10  # number of epochs. All input patterns are presented in one.
n_trials <- length(patterns)  # number of trials. One input pattern per trial.
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
        inputs <- list(patterns[[pat]][[1]], c(), c())

        for (minus in seq(n_minus)){
            # minus cycles: layer 1 is clamped
            net$cycle(inputs, clamp_inp = 1)
        }
        outs <- net$lays[[3]]$get_unit_acts() # saving the output for testing

        #+++++++ PLUS PHASE +++++++
        inputs <- patterns[[pat]]
        #list(unlist(patterns[[pat]][1], c(), unlist(patterns[pat, 2]))
        for (plus in seq(n_plus)){
            # plus cycles: layers 1 and 3 are clamped
            net$cycle(inputs, clamp_inp = 1)
        }

        #+++++++ LEARNING +++++++
        net <- net$chg_wt()  # updates the avg_l vars and applies XCAL learning

        #+++++++ ERRORS +++++++
        # Only the cosine error is used here
        errors[epoch, pat] <- 1 - cosine(patterns[[pat]][[3]], outs)
    }
    cat("\nepoch ", epoch, " finished")
}

## 6) A plot of the error by epoch figure
plot(rowMeans(errors, na.rm = T))
