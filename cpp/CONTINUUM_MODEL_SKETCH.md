# Minimal continuum reduction for load-stabilized assembly

## Particle-level constraints

The reduced model should preserve four observations rather than every particle
detail:

1. Motor load slows strong AM detachment approximately exponentially.
2. AA bond presence and force-bearing status are distinct variables.
3. Strong AM attachments establish force-bearing status, but their millisecond
   turnover makes instantaneous status too brittle. A fitted relaxation time is
   required; 5 ms is the provisional value from the two-actin assay.
4. Crosslinks affect AM kinetics through retained strain. There is no direct
   status-to-AM rate multiplier.

## Smallest closed state

Use fields or spatial averages

```text
m(x,t)       strong AM motor fraction
c(x,t)       total AA bond fraction
q(x,t)       force-bearing fraction of AA bonds, 0 <= q <= c
sigma(x,t)   active axial stress
S(x,t)       polarity/registry order parameter
```

If the nonmonotonic lambda response cannot be reproduced, split `c-q` into a
parallel/misregistered population `p` and an uncommitted population. That is the
first justified extension, because the particle runs implicate long-lived
status-1/parallel scaffolds in kinetic trapping.

## Candidate mean-field equations

```math
\dot m = k_{\rm on}^m a(t)(1-m)
       - k_{\rm off}^0 \exp[-|\sigma|/\sigma_*]m,
```

```math
\tau_q \dot q = c\,H(m,S,\sigma)-q,
```

```math
\dot c = k_{\rm on}^c(1-c)
       - \left[\frac{c-q}{\tau_1}
       + \frac{q}{\tau_1+\tau_2+\lambda |\sigma|}\right],
```

```math
\tau_\sigma \dot\sigma = \zeta a(t)m
       - \frac{\sigma}{1+\chi q},
```

```math
\dot S = \mu(c,q)\,K|\sigma|q(1-S)
       - D_S S.
```

`H` is the probability that two antiparallel actins are strongly attached to
distinct motors. Start with `H = h_0 m^2(1+alpha S)` and fit it from the event
tables. Set `tau_q = 5 ms` initially and retain 0, 1, and 20 ms as controls.
The stress equation implements mechanical feedback: force-bearing AA bonds slow
stress relaxation, and the retained stress then slows AM detachment.

## Capturing the intermediate-lambda maximum

The equations above give stabilization but may remain monotonic in `lambda`.
Two falsifiable closures can produce the observed optimum:

1. Turnover-limited search:

```math
\mu(c,q)=\mu_0/[1+\gamma q(\tau_2+\lambda|\sigma|)].
```

Large `lambda` preserves load-bearing structures but suppresses the rearrangements
needed to find the ordered state. This predicts history dependence: high lambda
preserves preassembled order but does not assemble efficiently from disorder.

2. Parallel scaffold trapping, if status 1 is stabilized:

```math
\dot p=k_p(1-S)(1-c-p)-p/\tau_p,
```

with `mu = mu_0/(1+gamma_p p)`. This predicts that status-2-only stabilization
removes the high-lambda suppression. Do not include `p` unless the toggle and
bond-type diagnostics require it.

## Parameter identification

- Fit `k_on^m`, `k_off^0`, and `sigma_*` from AM transition hazards.
- Fit `tau_q` from internal status-1 gaps and the memory sweep.
- Fit `tau_1`, `tau_2`, and `lambda` from conditional AA survival.
- Fit `chi` from stress relaxation with and without persistent AA crosslinks.
- Fit `K`, `D_S`, and the mobility closure from order trajectories and the
  medium-to-high-lambda history-switch test.

The key model-selection question is whether `(m,c,q,sigma,S)` reproduces the
lambda optimum and its history dependence. Only then decide whether the extra
parallel-trap field `p` is necessary.
