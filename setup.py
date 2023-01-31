import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

with open('requirements.txt') as fr:
    requirements = fr.read().splitlines()

setuptools.setup(
    name='bgcArgoDMQC',
    version='0.2.14',
    license='The MIT License (MIT)',
    author='Christopher Gordon',
    maintainer='cgrdn',
    author_email='chris.gordon@dfo-mpo.gc.ca',
    description='A python library for quality control of BGC-Argo data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ArgoCanada/bgcArgoDMQC',
    packages=setuptools.find_packages(),
    package_dir={'bgcArgoDMQC': 'bgcArgoDMQC'},
    install_requires=requirements,
    data_files=[('global_argo_index', ['bgcArgoDMQC/ref/ar_index_global_meta.txt.gz','bgcArgoDMQC/ref/ar_index_global_prof.txt.gz', 'bgcArgoDMQC/ref/argo_bio-profile_index.txt.gz', 'bgcArgoDMQC/ref/argo_synthetic-profile_index.txt.gz']),
                ('T_lL_tau_LUT',['bgcArgoDMQC/ref/T_lL_tau_3830_4330.dat']),
                ('resource_files', ['bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_BRtraj.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_meta.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_prof.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_Rtraj.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_Sprof.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/4901784_tech.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/profiles/BD4901784_000.nc', 'bgcArgoDMQC/resource/Argo/dac/meds/4901784/profiles/BD4901784_001.nc', 'bgcArgoDMQC/resource/NCEP/rhum/rhum.sig995.2019.nc', 'bgcArgoDMQC/resource/NCEP/rhum/rhum.sig995.2020.nc', 'bgcArgoDMQC/resource/WOA18/o2sat/ghost.txt'])],
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Development Status :: 3 - Alpha'
    ],
    python_requires='>=3.4',
)
