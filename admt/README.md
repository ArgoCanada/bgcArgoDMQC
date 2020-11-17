# Presentation Outline

- What is `bgcArgoDMQC`, and why did we make it?
- Basic information
- Status of the package
- Key features with examples
- Validation against various sources
- Planning ahead/future work
- Outreach

## What is `bgcArgoDMQC`

`bgcArgoDMQC` is a python package for performing quality control on biogeochemical
variables recorded by the BGC-Argo program. At this stage of development,
dissolved oxygen (DOXY) is the only variable this package works on, but in the
long term we would like to see it expand to include QC methods for all 6 of the
BGC-Argo variables. The DOXY quality control methods are heavily based on the
SOCCOM matlab code for calculating oxygen gain.

In terms of *why* `bgcArgo`, we seek to provide a tool to the community that is
reliable, agrees with previously established methods, and is open source.

## Basic Information

### Requirements

- Must run on `python3.4` or higher, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- TEOS-10 package [gsw](https://teos-10.github.io/GSW-Python/), but will also work with the [seawater](https://pypi.org/project/seawater/) package, though it is deprecated in favor of gsw
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) is required

### Optional packages

- [seaborn](https://seaborn.pydata.org/) is recommended but not required, through there will be some reduced (non-essential) functionality
- [cmocean](https://matplotlib.org/cmocean/) is also recommended for nicer plots, but not required

### Installation

The recommended install is through the conda-forge channel, via the command:

```bash
conda install -c conda-forge bgcargodmqc
```

The package is also available through the python package index <https://pypi.org/project/bgcArgo/>, install with:

```bash
pip install bgcArgoDMQC
```

## Status

As described above, the package currently only includes code to perform quality
control on oxygen data. The DOXY QC has been under development for about 8
months now, but is still in what I would descibe as the 'alpha' phase of
development. The package is functional, documented, and has undergone
significant validation testing, but is still under active development.

## Key Features

- Load in Argo data from synthetic or profile files:

```python
import bgcArgoDMQC as bgc

syn  = bgc.sprof(4902480)
prof = bgc.profiles([4902480, 4902481])
```

- Calculate oxygen gains using NCEP or WOA data:

```python
gains = syn.calc_gains() # defaults to in-air data
woa_gains = syn.calc_gains(ref='WOA') # using %-saturation

>>> print(f'In-air gains: {gains}')
>>> print(f'WOA-based gains: {woa_gains}')
```

- Visualize gains or otherwise!:

```python
>>> syn.plot('gains')
>>> syn.plot('cscatter', varname='DOXY')
```

- Compare to independent data:

```python
df1 = pd.read_csv('my_winkler_data.csv')
>>> df1.head()

syn.add_independent_data(
    PRES=df1.pressure, # data arguments, match naming to Argo variables
    DOXY=df1.dissolved_oxygen,
    date='2020-10-04', # metadata arguments, optional, if no date matches to first profile
    LATITUDE=df1.lat, LONGITUDE=df1.lon, # again, optional
    label='Winkler' # label to classify the data - for if you have more than one source
)

df2 = pd.read_csv('my_CTD_data.csv')
syn.add_independent_data(
  PRES=df2.PRES, DOXY=df2.DOXY, TEMP=df2.TEMP, PSAL=df2.PSAL,
  date='2020-10-04',
  label='Shipboard CTD'
)

>>> syn.compare_independent_data()
>>> syn.compare_independent_data(data='Winkler')
```

- Other features:
  - Caclulate gain with carryover factor and temperature dependence (Bittig et al. 2018)
  - Perform time-response correction on optodes
  - Convert between oxygen units (based on SCOR WG 142)
  - Download Argo, NCEP, WOA data

## Validation

## Future Work

## Outreach and Feedback

At the moment, this package is almost exclusively concerned with QC of oxygen
measurements, but long term we would like to add in methods for all BGC-Argo
variables. If you're interested in collaborating as (1) and alpha user/tester,
(2) contributing to a python tool for another variable or (3) generally helping
with development of the package, please get in touch!! My email is
<chris.gordon@dfo.mpo.gc.ca>.
