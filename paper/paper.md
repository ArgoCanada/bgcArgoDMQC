---
title: 'Performing Delayed-mode Quality Control on Argo Oxygen Measurements in Python'
tags:
  - Python
  - Argo
  - quality control
  - oceanography
  - dissolved oxygen
authors:
  - name: Christopher Gordon
    orcid: 0000-0002-1756-422X
    equal-contrib: true
    affiliation: 1
affiliations:
  - name: Bedford Institute of Oceanography, Department of Fisheries and Oceans Canada, Dartmouth, Canada
    index: 1
date: 4 August 2022
bibliography: paper.bib
---

# Summary

# Introduction

# Statement of need

The intention of this package is to provide the Argo community with an open
source package for performing delayed-mode quality control of BGC floats which
agrees with existing practices and softwares. It is not meant to superseed or
replace other softwares (e.g. 
[SAGE in matlab](https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING), 
@mauer:2021,  or [DM Filler in R](https://github.com/catsch/DM_FILLER)), but
rather to offer an option written in python for instances where programming
language preference or lack of access to propriatary software prevents use of
others. The primary function of the package is to perform DMQC on oxygen data
via gain adjustment and visual inspection, and then to export the adjusted
data to a delayed-mode file that complies with Argo file structure.

# Setup

# Using `bgcArgoDMQC`

# Conclusion

# Acknowledgements

# References
