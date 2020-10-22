Basic Funcationality
====================

`bgcArgo` is a python library of functions for quality controlling dissolved oxygen data. Heavily based on the <SOCCOM BGC Argo QC https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING>_ methods program in matlab, it uses either NCEP_ or `World Ocean Atlas`_ data to calculate oxygen gains (`Johnson et al. 2015`_).

.. _NCEP: https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html
.. _`World Ocean Atlas`: https://www.nodc.noaa.gov/OC5/woa18/
.. _`Johnson et al. 2015`: https://doi.org/10.1175/JTECH-D-15-0101.1

Although functions in the `bgcArgo` module may be of use in other situations, the majority of the functionality is lies within two classes, `profiles` for typical profile files and `sprof` for synthetic profiles.

.. code-block:: python
    import bgcArgo as bgc

    # setup for your system - these directories need to already exist!
    argo_path = 'your/argo/data/path' # where to save Argo data
    ncep_path = 'your/ncep/data/path' # where to save NCEP reanalysis data
    woa_path  = 'your/woa18/data/path' # where to save WOA data

    # download the data - this can take some time depending on connection
    # Argo
    wmos = [4902480, 6902905]
    dacs = [bgc.get_dac(w) for w in wmos]
    dacpath = '/ifremer/argo/dac'
    fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dacs, wmos)]
    # NCEP
    bgc.io.get_ncep('pres', local_path=ncep_path)
    bgc.io.get_ncep('land', local_path=ncep_path)
    # WOA
    bgc.io.get_woa18('O2sat', local_path=woa_path)

    # tell the package where to look for data
    bgc.set_dirs(argo_path=argo_path, ncep_path=ncep_path, woa_path=woa_path)
    # load data from profiles for two floats
    flts = bgc.profiles(wmos)
    # calculate the dissolved oxygen gains
    gains = flts.calc_gains()
    print(gains)

    # load a synthetic profile
    syn = bgc.sprof(4902480)
    # plot a time vs. depth section for the top 500m
    g1 = syn.plot('cscatter', varname='DOXY', ylim=(0,500))
    # plot the first 10 profiles for temperature, practical salinity,
    # and adjusted oxygen
    g2 = syn.plot('profiles', varlist=['TEMP','PSAL', 'DOXY'], Ncycle=10, Nprof=10, ylim=(0,500))

Version History
^^^^^^^^^^^^^^^

- 0.1: April 20, 2020 - Initial creation
- 0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class
- 0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient
- 0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe
- 0.2.3 - 0.2.6: September 3, 2020 - updates to pass all checks on conda-forge pull request, updated on PyPI as well
- 0.2.7 - 0.2.8: September 29, 2020 - re-spun for PyPI and PR to conda-feedstock