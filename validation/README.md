# validation

Storage of code relating to validation of the `bgcArgo` module

## comparison to direct SAGE-O2 output

## comparison to completed Argo DMQC

## comparison to DOXY audit distributed by Josh Plant

Comparisons to the DOXY audit were performed using `sprof(wmo)` and matching
the cycle to the distributed audit gain calculated presumably using SAGE-O2
code, if not exactly the GUI distributed by SOCCOM, using analogous code.

The latest DOXY audit contained 5399 profiles and gains. Of these comparisons
made, 3086 (1560) gains were within 0.01 of one another (effectively identical, as
the DOXY audit gains have a precision of 0.01), 1237 (2333) were greater than or
equal 0.01 but less than 0.05, 360 (514) were between 0.05 and 0.2, and 81 (136) were
greater than 0.2. There were also 9 (225) gains that were NaN valued (NaN returned
by the python package, no NaN values are reported in the audit), and 626 (631) were
infinite valued by both the python package and in the DOXY audit (there are
634 infinite values in the DOXY audit, so still missing 3 matches).

And this shows the success of comparisons as a percent of total comparisons:

![Percent of comparison](https://raw.githubusercontent.com/ArgoCanada/BGC-QC/master/figures/doxy_audit/DOXY_audit_comparison_breakdown_20200920.png)

### change log

July 19, 2020: Plan moving forward is to look into the NaN-valued and high (> 0.2) difference
gains to see if there is anything consistent about those data or metadata
that are causing the mismatch.

As of July 24, 2020 the comparison was greatly improved (for more information
see [issue #10](https://github.com/ArgoCanada/BGC-QC/issues/10)). Figures have
been updated to reflect this, and plan still stands to investigate large and
NaN valued ones until we have full agreement.

As of July 30, 2020 all NaN values have been removed - bgcArgo no longer is
returning any NaN values during the analysis. Figures have been updated and
more information is available in
[issue #10](https://github.com/ArgoCanada/BGC-QC/issues/10).
