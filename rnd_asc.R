# This program constructs a 3-layer random associator using the 'network',
# 'layer', and 'unit' objects. It is a tiny test of these classes.

# To reproduce the example we can use a seed. You can try to guess whose
# birthday is on July 22nd, 1904.
set.seed(07221904)

# Specify the dimensions and connectivity of the layers. At this point, layers
# are either fully connected or unconnected.

# Specifying 3 layers and their dimensions.
dim_lays <- list(c(2, 5), c(2, 10), c(2, 5))

# Specify connection strength between layers, if layer j sends projections to
# layer i, then cxn[i, j] = strength > 0 and 0 otherwise. Strength specifies the
# relative strength of that connection with respect to the other projections to
# layer i.
cxn <- matrix(c(0, 0, 0,
                1, 0, 0.2,
                0, 1, 0), nrow = 3, byrow = T)

# Create the network with default parameters; dim_lays and cxn is the minimum
# you need to specify a network, but if constructing other networks you should
# pay careful attention to g_i_gain, which controls overall inhibition in a
# layer. If this value is not set carefully, you will likely not get what you
# want.
net <- network$new(dim_lays, cxn)

# Create 15 random inputs with the convenience function create_inputs in the
# network class. Of course you could do this on your own.
inputs_plus <- net$create_inputs(which_layers = c(1, 3),
                                 n_inputs = 15,
                                 prop_active = .3)

# Remove the inputs for the output layer (3) for minus phase.
inputs_minus <- lapply(inputs_plus, function(x) {x[3] <- list(NULL); return(x)})
n_epochs <- 10

# Learn with default parameters, return value is the output activation after
# every trial, which we will use to calculate the error.
outs <- lapply(seq(n_epochs), function(x) net$learn_error_driven(inputs_minus,
                                                                 inputs_plus,
                                                                 lrate = 0.5))
# The network class can calculate the mean absolute deviation for each epoch.
# You can also use your own functions on these lists to calculate other types of
# errors.
mad <- net$mad_per_epoch(outs, inputs_plus, 3)
plot(mad)
