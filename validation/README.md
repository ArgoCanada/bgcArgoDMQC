# validation

Storage of code relating to validation of the `bgcArgo` module

## Comparison to direct SAGE-O2 output

## Comparison to completed Argo DMQC

## Comparison to DOXY audit distributed by Josh Plant

Comparisons to the DOXY audit were performed using `syn(wmo)` and matching
the cycle to the distributed audit gain calculated presumably using SAGE-O2
code, if not exactly the GUI distributed by SOCCOM, using analogous code.

The latest DOXY audit contained 5399 profiles and gains. Of these comparisons
made, 720 gains were within 0.01 of one another (effectively identical, as
the DOXY audit gains have a precision of 0.01), 1880 were greater than or
equal 0.01 but less than 0.05, 967 were between 0.05 and 0.2, and 496 were
greater than 0.2. There were also 804 gains that were NaN valued (NaN returned
by the python package, no NaN values are reported in the audit), and 532 were
infinite valued by both the python package and in the DOXY audit.

This figure summarizes all the 5000+ comparisons made:

![All comparisons](https://raw.githubusercontent.com/ArgoCanada/BGC-QC/master/figures/doxy_audit/DOXY_audit_comparison_waffle.png)

And this shows the success of comparisons as a percent of total comparisons:

![Percent of comparison](https://raw.githubusercontent.com/ArgoCanada/BGC-QC/master/figures/doxy_audit/DOXY_audit_comparison_breakdown.png)

Plan moving forward is to look into the NaN-valued and high (> 0.2) difference
gains to see if there is anything consistent about those data or metadata
that are causing the mismatch.
