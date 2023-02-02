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
    include_package_data=True,
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
