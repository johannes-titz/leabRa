---
title: "Leabra Project Documentation"
output: github_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# todos
rechtschreibung in docu!

test building with cran options, go through karl broman guide

Sergio schreiben

https://cran.r-project.org/web/packages/URL_checks.html

test building on windows

check examples and time for running them (devtools?, see karl broman)

The ownership of copyright and intellectual property rights of all components of the package must be clear and unambiguous (including from the authors specification in the DESCRIPTION file). Where code is copied (or derived) from the work of others (including from R itself), care must be taken that any copyright/license statements are preserved and authorship is not misrepresented.
Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the authors of such code. Alternatively, the ‘Author’ field should list these authors as contributors.

Where copyrights are held by an entity other than the package authors, this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’ field, or using a ‘Copyright’ field (if necessary referring to an inst/COPYRIGHTS file).

Trademarks must be respected.

1. add type of layer for avg_l_lrn, it is 0.0004 on default, but 0 for output layers (they do not need self-organized learning).

5. Unit-testing.

2. Find a general logistic function that approximates the nxx1 or generate a smoother nxx1_df? Right now, I think the step function will be ok in general.
4. Checkout hogging-modification on leabra-ccn site when self-organizing is ready. I think it works already ok with g_i_gain set correctly.

spiking:

it is very strange that e_rev_e is set to 1, so that v_m (at least i think)
cannot grow larger than 1, but spike_thr is 1.2, so there will never be a spike
right; don't know

## errors

### weight bounding
I believe that in the matlob code there is a tiny mistake, on line 276 in the network class:

```{matlab}
idxp = logical(this.layers{rcv}.wt > 0);
``` 

This should be .dwt ant not .wt; depending on whether weights are increasing or decreasing we have to apply different calculations. The weight itself is always between 0 and 1. It is correct on this page: https://grey.colorado.edu/emergent/index.php/Leabra

### unused vars
* in unit class l_dn_dt is never used
* never used in unit  l_up_inc = 0.2  (increase in avg_l if avg_m has been "large")

## Problems
### links in matlab docu

https://grey.colorado.edu/ccnlab/index.php/Leabra_Hog_Prob_Fix#Adaptive_Contrast_Impl

### why is spk_thr 1.2; this is not a v value of -50mV...

### ext_input is activation and g_e
I do not really understand the part with clamped_input and ext_input

In the network class we can specify inputs and a binary flag whether inputs are
clamped or not.

We have different cases: 
the normal case is that we have input and it is clamped, then we set the input (activation) to this value without time integration. of course we can have layers that do not have any input, so for these after setting the clamped activations we can use the intern_input coming from these clamped layers to the no-input-layers and calculate their activation.

but what happens if we do not clamp inputs at all? Then, itnerestingly, cycle layer method is called with intern_inputs and extern_inputs
ce_wt \* intern_inputs is the g_e for the units, but now ce \*wt and extern_inputs are simply added together, so ext_inputs is a g_e too, although it is an activation. Do they have the same range? well in ccnbook they state that all parameters are approx between 0 and 1, still the parameter plays two roles. 

### why can avg_l_max be larger than 1? avg_l_max = 1.5 (max value of avg_l)
This is irrelevant in the current emergent version, because there is no max for avg_l anymore.

## errors in book
https://grey.colorado.edu/CompCogNeuro/index.php/CCNBook/Neuron adapt.dt_time is 144ms, but it should be 1/144ms

((avg_l.lrn_max - avg_l.lrn_min) / avg_l.gain - avg_l.min)) 

should probably be in braces: (avg_l.gain - avg_l.min)

verduzco, name correct on matlab site?

## Tested:
1. Vectorize xx1 (this does not bring any gains)

2. a new cxn class, this is just too much work right now
What will we need for this?

Note that the layer cycle function now needs the g_e_per_unit value as a parameter. Before it was calculated there, but now the layer class does not have the weights anymore, so it has to be done in the network class. There the g_e_per_unit is calcualted with the cxn-class and then given to the layer cycle class.

## Calculating contrast enhanced weights

The formula for contrast-enhanced weights in the learning chapter looks slightly different, but is actually equivalent if you look closely. In the code it is

```{r}
1 / (1 + (off * (1 - self$weights) / self$weights) ^ gain)
```

in the book the gain is negative and the fraction is turned around: 

```{r}
1 / (1 + (self$weights) / (off * (1 - self$weights)) ^ -gain)
```


# Changes in emergent since matlab-version (April 2015)

This is the comparison between the current version (as of 2017-09-07) and a version of Feb 2015 (the last one before April 2015, the matlab version release)
https://grey.colorado.edu/emergent/index.php?title=Leabra&type=revision&diff=13089&oldid=11046

I will only look at the critical parts:

I_net is reset when spiking, this is a minor thing and will not change anything I suppose, because i_net is calculated without the previous i_net value. I have adopted it anyway.

Weight updating seems to be like in matlab, although the docu is now more detailed. All other things seem to be optional.

The major changes will be described in the following sections:

## I_net_r 

I_net_r = net \* (e_rev.e - v_m_eq) + gc_l \* (e_rev.l - v_m_eq) + gc_i \* (e_rev.i - v_m_eq) + noise

rate-coded version of I_net, to provide adequate coupling with v_m_eq. 

v_m_eq +=  dt.integ \* dt.vm_dt \* (I_net_r - adapt)

In matlab version it is v_h (half step) instead of v_eq.

I now use this version and it seems to be almost identical to the matlab one.

## avg_l
avg_l is updated differently now, but it was already different in the matlab version anyway. Now, there is no max value for avg_l, only a gain factor, that is multiplied with avg_m to update avg_l. avg_l_prc is not needed anymore, it seems.

avg_l_prc is also used in layer, but we can simply use avg_l instead. 

So I just made the switch and it looks ok so far.

## avg_l_lrn
avg_l_lrn also depends on cosine between minus and plus, which i do not really like. We can use 0.0004 instead, but only for hidden layers; output layers do not need hebbian learning (see lebra docu in ccnbook). Add type of layer someday.

## momentum
I will not include momentum, because in my opinion it violates the "local" approach of leabra. Although there might be some biologically plausible implementations of it.

# Simulation stuff

### input
We will need a function that generates the input, checkout the pass-t repo

 * local inputs with one active unit, what about 100 unit layer 1 (10x10s, 100 stimuli) and layer 2?
 * check out norman \& Oreilly, how large were their networks?

### What to simulate?

check what happens to output competitive layer, when we present a stimulus
repeatedly. expectation: variance in activation increases over time
plot: variance in activation as a function of number of chg_wt()
one stimulis is ok here i suppose

then use learning weight to have freq-dominance
expectation: variance in activation increases over time, but non-linearly
plot: the same as above, weighted with learning rate...approx

does activation go down in layer overall? if not, increase gi_gain, if it is high enough, there should be a decrease!
plot: overall activation in network for stimulus A as a function of number of chg_wt()

attractor, how long does it take to go into attractor, we need a function for this
expectation: number of cycle reduces the more often a stimulus is presented

optional: attention, through avg_l as the threshold the xcal function should lead to less learning if a stimulus is presented for a continous time
plot: d_wt/d_t as a function of avg_l, on unit level? maybe on layer level? mean(d_wt/d_t) as a function of chg_wt()...the more often the shorter the change? it might oscillate at some point; would it be better to show something different so that the avg_l_receiver can take a break...?

gain in sigmoid function as a regulator of norepinephrin, stress...,
expectation: larger gain, produces more learning?
