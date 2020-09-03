import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

with open('requirements.txt') as fr:
    requirements = fr.read().splitlines()

setuptools.setup(
    name='bgcArgo',
    version='0.2.5',
    license='The MIT License (MIT)',
    author_email='chris.gordon@dfo-mpo.gc.ca',
    description='A python library for quality control of BGC-Argo data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ArgoCanada/BGC-QC',
    packages=setuptools.find_packages(),
    package_dir={'bgcArgo': 'bgcArgo'},
    install_requires=requirements,
    data_files=[('global_argo_index', ['bgcArgo/ref/ar_index_global_meta.txt.gz','bgcArgo/ref/ar_index_global_prof.txt.gz', 'bgcArgo/ref/argo_bio-profile_index.txt.gz', 'bgcArgo/ref/argo_synthetic-profile_index.txt.gz']),
                ('T_lL_tau_LUT',['bgcArgo/ref/T_lL_tau_3830_4330.dat'])],
    classifiers=[
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Development Status :: 3 - Alpha'
    ],
    python_requires='>=3.4',
)
