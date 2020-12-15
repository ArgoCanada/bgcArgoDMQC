from pathlib import Path

global config_file
pkg_path = Path(__file__).parent.absolute()
config_file = pkg_path / '.config'

def check_config():
    # make sure its there, even as an empty file, eases handling of FileNotFoundError's
    if not config_file.exists():
        config_file.touch()

def read_config():
    # dict object to load values into
    config_dict = dict()
    # loop through the config file - if there are no lines, will return an
    # empty dict, which is fine/right
    with open(config_file) as f:
        for line in f:
            llist = line.split('=')
            key = llist[0]
            val = llist[1]
            config_dict[key] = val
    
    return config_dict

def configure(argo_path=None, ncep_path=None, woa_path=None):
    # existing and new configuration info
    existing_config = read_config()
    new_config = dict(argo_path=argo_path, ncep_path=ncep_path, woa_path=woa_path)

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