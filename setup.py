#!/usr/bin/env python

from setuptools import setup, find_packages
import os
import multiprocessing, logging  # AJ: for some reason this is needed to not have "python setup.py test" freak out

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='MatMethods',
        version='0.0.2',
        description='MatMethods has implementations of FireWorks workflows for Materials Science',
        long_description=open(os.path.join(module_dir, 'README.rst')).read(),
        url='https://github.com/hackingmaterials/MatMethods',
        author='Anubhav Jain',
        author_email='anubhavster@gmail.com',
        license='modified BSD',
        packages=find_packages(),
        package_data={},
        zip_safe=False,
        install_requires=['FireWorks>=1.1.4', 'pymatgen>=3.2.1', 'custodian>=0.8.8'],
        extras_require={},
        classifiers=['Programming Language :: Python :: 2.7',
                     'Development Status :: 4 - Beta',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Intended Audience :: Information Technology',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        test_suite='nose.collector',
        tests_require=['nose'],
        scripts=[]
        #scripts=[os.path.join('scripts', f) for f in os.listdir(os.path.join(module_dir, 'scripts'))]
    )
