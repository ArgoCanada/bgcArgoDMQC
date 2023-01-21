from pathlib import Path

global config_file
pkg_path = Path(__file__).parent.absolute()
config_file = pkg_path / '.config'

def check_config():
    '''
    Check that configuration file exists, and if not, create a blank one. 
    Having a blank file exist eases throwing any FileNotFoundError's, but
    if the config is blank then still no configuration is performed.
    '''
    # make sure its there, even as an empty file, eases handling of FileNotFoundError's
    if not config_file.exists():
        config_file.touch()

def read_config():
    '''
    Read in permanently stored configuration information for locations of
    Argo, NCEP, or WOA data.
    '''
    # dict object to load values into
    config_dict = dict()
    # loop through the config file - if there are no lines, will return an
    # empty dict, which is fine/right
    with open(config_file) as f:
        for line in f:
            llist = line.split('=')
            key = llist[0]
            val = llist[1].strip()
            config_dict[key] = val
    
    return config_dict

def reset_config():
    '''
    Permanently erase .config data, replace with empty file.
    '''
    config_file.unlink()
    config_file.touch()

def configure(**kwargs):
    '''
    Set up locations for Argo, NCEP, and/or WOA data on local machine.
    '''

    # existing and new configuration info
    existing_config = read_config()
    new_config = dict(**kwargs)

    paths = [k for k in new_config.keys() if "path" in k]
    for p in paths:
        new_config[p] = Path(new_config[p]).as_posix()

    # get rid of None values
    clean_new_config = {k: v for k, v in new_config.items() if v is not None}

    # config info already in .config file that should be kept
    keep_config = set(existing_config.keys()) - set(clean_new_config.keys())

    # config info that either exists in .config and should be overwritten, or is new to .config
    overwritten_config = set(existing_config.keys()).intersection(set(clean_new_config.keys()))
    add_config = set(clean_new_config.keys()) - set(existing_config.keys())
    write_config = overwritten_config | add_config

    final_config = dict()
    for k in keep_config:
        final_config[k] = existing_config[k]
    for k in write_config:
        final_config[k] = clean_new_config[k]

    # write the .config file anew
    with open(config_file, 'w') as f:
        i = 0
        for k, v in final_config.items():
            if i > 0:
                f.write('\n')
            f.write('{}={}'.format(k, v))
            i += 1