# to reproduce example
set.seed(22071904)

# to later calculate error we use the cosine similarity
cosine <- function(x, y){
  sum(x * y) / (sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2)))
}

# This program constructs a 3-layer random associator using the 'network',
# 'layer', and 'unit' objects. It is a test of these classes, and a tutorial in how to construct networks with them.

# 1) First, specify the dimensions and connectivity of the layers At this
# point, layers are either fully connected or unconnected

# 1.1) Set the dimensions of the layers specifying 3 layers and their dimensions
dim_lays <- list(c(1, 25), c(1, 100), c(1, 25))

# 1.2) Specify connectivity between layers
cxn <- matrix(c(0, 0, 0,
                1, 0, 0.2,
                0, 1, 0), nrow = 3, byrow = T)

# cxn(i,j) = c means that layer i receives cxn from j, and that
# they have a relative strength c. In this case, feedback cxn are 5
# times weaker. The relative weight scale of layer i comes from the non-zero
# entries of row i.  The network constructor will normalize this matrix so that
# if there are non-zero entries in a row, they add to 1.

net <- network$new(dim_lays, cxn, g_i_gain = c(2, 2, 2))

inputs_plus <- net$create_inputs(which_layers = c(1, 3),
                              n_inputs = 15, prop_active = .3)

# remove the inputs for the output layer for minus phase
inputs_minus <- lapply(inputs_plus, function(x) {x[3] <- list(NULL); return(x)})
n_epochs <- 10

outs <- lapply(seq(n_epochs), function(x) net$learn_error_driven(inputs_minus, inputs_plus, lrate = 0.8))

# calc cosine

extract_list_number = function(x, y){
  lapply(x, function(x) x[[y]])
}

lapply(outs, extract_list_number, 3)

overall_cosine <- function(outs, inputs_plus, layer){
  outs_layer <- extract_list_number(outs, layer)
  inputs_plus_layer <- extract_list_number(inputs_plus, layer)
  mean(unlist(Map(function(x, y) 1 - cosine(x, y), outs_layer, inputs_plus_layer)))
}

cos_error <- sapply(outs, overall_cosine, inputs_plus, 3)
plot(cos_error)
