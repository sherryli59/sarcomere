# A–A / A–M overlap mechanism sandbox

See [REPORT.md](REPORT.md) for the conclusions and
[VARIABLE_LENGTH_TEST_PLAN.md](VARIABLE_LENGTH_TEST_PLAN.md) for the proposed
full-simulator diagnostic.  No production run was launched.

This directory is an isolated, reduced stochastic model used to discriminate
candidate mechanisms before changing the production C++ simulator.

It retains four mobile axial bodies (two antiparallel actin barbed ends and two
neighboring myosin centers), overlap-dependent active motor force, ACTN2
turnover and catch stabilization, translational Brownian diffusion, and
reciprocal mechanical forces.

Candidate interactions:

- turnover-capable elastic ACTN2;
- permanent elastic ACTN2 control;
- Johann-like overlap adhesion, `U = -gamma * A–A overlap`;
- Johann-like fixed-pool entropic force, `F = N kBT / overlap`;
- persistent Z-disc axial register;
- titin-like Z-disc–M-line springs;
- fixed-M-line boundary control;
- combined mechanisms.

Three levels are included:

1. `overlap_model.py`: a local four-body free model for force signs and creep;
2. `periodic_model.py`: the exact axial 2-myosin/4-actin force-balanced
   `two_sarcomere_unit` topology in a fixed box;
3. `variable_cell_model.py`: the same topology with one diagnostic
   variable-repeat coordinate, so global A–M selection can actually be tested.

These models are not replacements for the full simulator. Their job is to
reject mechanisms with the wrong force sign, unavoidable creep, or a remaining
free coordinate before those mechanisms are ported.

Run:

```bash
python -m unittest -v test_overlap_model.py test_periodic_models.py
python run_campaign.py
python run_periodic_campaign.py
python run_variable_cell_campaign.py
python make_periodic_animations.py
python make_variable_cell_animations.py
python make_phase_portraits.py
```

The campaigns write CSV tables, PNG summaries, and GIF animations under
`results/`.
