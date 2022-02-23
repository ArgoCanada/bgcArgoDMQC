#!/usr/bin/python

from pathlib import Path
import bgcArgoDMQC as bgc

# this is a vignette of how I would like netcdf output to work
# many (most?) may not work or even exist yet, but I will try to 
# write the simplest script I can, and then create functions that
# allow that script to be carried out

# designate float ID, load the float
my_wmo = 4900497
syn = bgc.sprof(my_wmo)

# calculate gain using WMO
gain = syn.calc_gains(ref='WOA').mean()
### ACTION ITEM - make calc_gains() return a SERIES instead of a NUMPY ARRAY so if can accept trailing mean() action

# populate DOXY_ADJUSTED
# note - this returns doxy_adjusted, but it also does some work within the object, the output variable is not required
doxy_adjusted = syn.populate_doxy_adjusted(gain)
### ACTION ITEM - it may be prudent to allow this to take a value for drift (slope over time) or an array of gains

# populate adjusted error - by default do pO2=4, but allow other input. May use other methods later using kwargs
# same as above, returned argument is not necessary down the line, saved within sprof object
doxy_adjusted_error = syn.populate_doxy_adjusted_error(pO2=4)
### ACTION ITEM

# adjust flags away from 3, then check TEMP_ADJUSTED_QC and PSAL_ADJUSTED_QC
new_flags = syn.populate_doxy_adjusted_flags()
# adjust flag values that may be due to visual QC or any other non-computational reason
index = [0, 1, 2, 101]
values = 4*[4]
new_new_flags = syn.adjust_flags(index, values)
### ACTION ITEM

# create comments
comments = syn.generate_comments_equations(
    variable='DOXY',
    gain=gain,
    operator='Christopher Gordon',
    affiliation='Fisheries and Oceans Canada',
    orcid='0000-0002-1756-422X'
)
### ACTION ITEM - all arguments should also be able to be set to defaults using 'configure'

# change the data mode, input is list of parameters and mode to be made into
syn.update_data_mode(['DOXY'], ['D'])
### ACTION ITEM

# note that if all the above has been performed already, then the syn object knows about all this, so you can call it on its own
# alternatively, you could call it out of nowhere, but that should be done with caution and appropriate verification should be done afterwards
# this is the recommended use, if the above have all been ran
syn.Dmode_export()

# this method for running allows for custom things to be added
syn.Dmode_export(
    doxy_gain=gain,
    # the following arguments are all equivalent in function - to select which files get processed up to D-mode, using more than one causes error
    cycles=list(range(1, 10)),
    datelim=('2021-01-01', '2021-07-01'),
    Rfiles=list(Path('/path/to/my/argo/data/dac/wmo').glob('BR*.nc')),
    # varibles to be added
    add_variables=dict(DOXY_ADJUSTED=doxy_adjusted, DOXY_ADJUSTED_ERROR=doxy_adjusted_error),
    add_flags=dict(DOXY_ADJUSTED_ERROR=doxy_adjusted_error)
    # comments is a dict that can be expanded like kwargs with the default values
    **comments,
    # where to save files - saves to where R files are by default
    output_dir=Path('/where/to/output/Dfiles/')
)

# run the file checker
# calling Dewey's wrapper, redundancy to not run if java not installed?
### ACTION ITEM