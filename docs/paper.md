---
title: 'bgcArgoDMQC: a python library for performing quality control on biogeochemcial (BGC) Argo floats'

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

date: 2 Feb 2021

bibliography: paper.bib

---

# bgcArgoDMQC: a python library for performing quality control on biogeochemcial (BGC) Argo floats

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

The ocean is difficult to observe, but has become more accessible over the past
20+ years thanks to autonomous platforms like floats and gliders. The Argo
program is perhaps the best implementation of autonomous technology for
observing the ocean, boasting over 3900 active floats (as of February 2, 2021,
via <ocean-ops.org>) and over 2 million profiles collected over its history.
Argo began as observations of exclusively temperature and salinity (the "core"
program). In **year** biogeochemical observations were introduced (the "BGC"
program), introducing measurements of chlorphyll fluorescence, optical
backscattering, dissolved oxygen, pH, and irradiance to the program.

For each variable measured by an Argo float, there is a corresponding series
of quality control (QC) methods. These methods are split into real-time (RT)
and delayed-mode (DM) methods. Real time quality control is automatically
conducted as data from the floats arrives, and is typically complete in 24
hours or less. Delayed mode quality control is conducted by a DMQC operator
usually after a profile is more than 1 year old. QC methods for RT and DM data
are different for each variable depending on issues with that sensor, but a
visual check is always included in DMQC.

Due to the global community nature of the Argo program, a wide variety of
software has been developed for performing QC on Argo data. This paper deals
specifically with delayed-mode quality control, and so will discuss software
designed for performing DMQC only. The dominant programming language in the
early aughts of the Argo program was MATLAB, which is a proprietary software.
The accepted methods for performing quality control on salinity is using the
Owens-Wong-Cabanes method (**OWC REFERENCES**) in MATLAB, but this has recently
been translated to python (see <https://github.com/euroargodev/argodmqc_owc>).
The software presented here is an effort to provide a complementary package
for performing DMQC on BGC-Argo variables.

The package, named `bgcArgoDMQC`, began as a python translation of MATLAB code
written by Mauer & Plant (<https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING>)
for correcting oxygen, particularly calculating oxygen gain following Johnson
et al. (**2014**). Oxygen quality control is currently the package's primary
function, and thus will be the primary subject of this paper, however it will
develop into a package that is capable of performing DMQC on all BGC variables.
This package seeks to fill a gap in the current array of available software
by offering the same high-quality, reliable QC methods in the MATLAB version
in an open-source language.

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
