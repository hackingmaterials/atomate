#!/usr/bin/env python

import os

from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='atomate',
        version='0.8.3',
        description='atomate has implementations of FireWorks workflows for '
                    'Materials Science',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/hackingmaterials/atomate',
        author='Anubhav Jain, Kiran Mathew',
        author_email='anubhavster@gmail.com, kmathew@lbl.gov',
        license='modified BSD',
        packages=find_packages(),
        package_data={'atomate.vasp.workflows.base': ['library/*']},
        zip_safe=False,
        install_requires=['FireWorks>=1.4.0', 'pymatgen>=2018.6.11',
                          'custodian>=2018.6.11', 'monty>=1.0.2',
                          'tqdm>=4.7.4', 'six',
                          'pymatgen-diffusion>=2018.1.4',
                          'pydash>=4.1.0'],
        extras_require={'rtransfer': ['paramiko>=2.4.1'],
                        'plotting': ['matplotlib>=1.5.2'],
                        'phonons': ['phonopy>=1.10.8'],
                        'complete': ['paramiko>=2.4.1',
                                     'matplotlib>=1.5.2',
                                     'phonopy>=1.10.8']},
        classifiers=['Programming Language :: Python :: 2.7',
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.6",
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
