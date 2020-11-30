#!/usr/bin/python

from pathlib import Path

import pandas as pd
import bgcArgoDMQC as bgc

# load in DOXY audit file - most recent one was on July 9 2020
audit_file = Path('/Users/gordonc/Documents/argo/doxy-audit/DOXY_audit_112020.TXT')
df = pd.read_csv(audit_file, sep='\t', header=25)

# download the synthetic, meta, and BRtraj files for each float in the audit
dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, bgc.get_dac(w), w) for w in df.WMO.unique()]
local_path = '/Users/gordonc/Documents/data/Argo'
bgc.io.get_argo(fltpath, local_path=local_path, mode='summary')
bgc.set_dirs(argo_path='/Users/gordonc/Documents/data/Argo', woa_path='/Users/gordonc/Documents/data/WOA18')

with open(Path('../data/doxy_audit_vs_bgcArgo_py_comparison_20200920.csv'),'w') as fid:
    fid.write('WMO,CYCLE,DAC,DATE,pyGAIN,sageGAIN')
    for wmo in df.WMO.unique():
        sub = df[df.WMO == wmo]
        syn = bgc.sprof(wmo)
        # syn.clean(bad_flags=4)
        syn.check_doxy_range()
        syn.calc_gains(ref='WOA')
        for i in range(sub.shape[0]):
            cycle = sub.cycle.iloc[i]
            ix = syn.CYCLE == cycle
            print('(pyGain, matlabGain) = ({:.2f}, {:.2f})'.format(syn.gains[ix][0], sub.G_raw.iloc[i]))
            fid.write('\n{:d},{:d},{:s},{:s},{:.2f},{:.2f}'.format(wmo,cycle,sub.DAC.iloc[i],sub['profile date'].iloc[i],syn.gains[ix][0],sub.G_raw.iloc[i]))
