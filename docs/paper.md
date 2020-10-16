---
title: 'bgcArgo: a python library for performing quality control on biogeochemcial (BGC) Argo floats'
tags:
  - Python
  - ocean
  - oceanography
  - observation
authors:
  - name: Christopher Gordon
    orcid: 0000-0002-1756-422X
    affiliation: 1
affiliations:
  - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada, Dartmouth, NS
    index: 1
date: 22 Oct 2020
bibliography: paper.bib
---

# bgcArgo: a python library for performing quality control on biogeochemcial (BGC) Argo floats

## Summary

Note - summary is like the abstract. The information within should be repeated
in the body of the paper.

- What is Argo?
- What is delayed mode? What is quality control?
- What is BGC data?

## Introduction

- Basic Argo motivation (ocean is difficult and costly to observe, etc.)
- The Argo network, and the biogeochemcial extension of the core mission
- High level quality control description
- Note current gap in community for an open source tool to QC BGC data

## Key Features

- Calculate gain as a simple ratio as in Johnson et al. (2015) and SAGE-O2
- Use in-air or surface data to calculate gain, comparing to NCEP or WOA data
- Calculate gain using carry-over factor (Bittig et al. 2018)
- Output required variables and information for delayed-mode Argo files
- Sufficient plotting features to do visual QC

## Conclusion

- Re-emphasize ability to calculate gain using various sources
- This package goes beyond what SAGE-O2 does already, adds features
- Provides an open source tool to the community where only a proprietary one previously existed

## Acknowledgements

Acknowledgements to developers of SAGE-O2 (especially Josh Plant, Tanya Mauer).
To Henry Bittig for valuable feedback in early stages. To to any expert alpha
users as well.

## References

Bibliography generated from papers.bib goes here.
