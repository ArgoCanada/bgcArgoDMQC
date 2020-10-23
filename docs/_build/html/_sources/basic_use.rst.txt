Basic Funcationality
====================

`bgcArgo` is a python library of functions for quality controlling dissolved oxygen data. Heavily based on the `SOCCOM BGC Argo QC`_ methods program in matlab, it uses either NCEP_ or `World Ocean Atlas`_ data to calculate oxygen gains (`Johnson et al. 2015`_).

.. _NCEP: https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html
.. _`World Ocean Atlas`: https://www.nodc.noaa.gov/OC5/woa18/
.. _`Johnson et al. 2015`: https://doi.org/10.1175/JTECH-D-15-0101.1
.. _`SOCCOM BGC Argo QC`: https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING

Although functions in the `bgcArgo` module may be of use in other situations, the majority of the functionality is lies within two classes, `profiles` for typical profile files and `sprof` for synthetic profiles.

.. literalinclude:: basic_use.py

Version History
^^^^^^^^^^^^^^^

- 0.1: April 20, 2020 - Initial creation
- 0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class
- 0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient
- 0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe
- 0.2.3 - 0.2.6: September 3, 2020 - updates to pass all checks on conda-forge pull request, updated on PyPI as well
- 0.2.7 - 0.2.8: September 29, 2020 - re-spun for PyPI and PR to conda-feedstock