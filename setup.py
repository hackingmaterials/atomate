#!/usr/bin/env python

import os

from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='atomate',
        version='1.0.1',
        description='atomate has implementations of FireWorks workflows for '
                    'Materials Science',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/hackingmaterials/atomate',
        author='Anubhav Jain, Kiran Mathew',
        author_email='anubhavster@gmail.com, kmathew@lbl.gov',
        license='modified BSD',
        packages=find_packages(),
        package_data={'atomate.vasp.workflows.base': ['library/*'],
                      'atomate.vasp.builders': ['*', 'examples/*']},
        zip_safe=False,
        install_requires=[
            'numpy',
            'scipy',
            'FireWorks>=1.4.0',
            'pymatgen-analysis-diffusion>=2021.4.29',
            'monty>=2.0.6',
            'paramiko',
            'pandas',
            'tqdm>=4.7.4',
            'networkx',
            'pymatgen>=2020.9.14',
            'custodian>=2019.8.24',
            'pydash>=4.1.0',
            'pyyaml>=5.1.2',
            'maggma>=0.26.0'
        ],
        extras_require={'rtransfer': ['paramiko>=2.4.2'],
                        'plotting': ['matplotlib>=1.5.2'],
                        'phonons': ['phonopy>=1.10.8'],
                        'complete': ['paramiko>=2.4.2',
                                     'matplotlib>=1.5.2',
                                     'phonopy>=1.10.8']},
        classifiers=["Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.6",
                     "Programming Language :: Python :: 3.7",
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Intended Audience :: Information Technology',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        test_suite='nose.collector',
        tests_require=['nose'],
        scripts=[os.path.join('scripts', f) for f in os.listdir(os.path.join(module_dir, 'scripts'))]
    )
