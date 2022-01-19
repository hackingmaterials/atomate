import os
from os.path import abspath, dirname, join

from setuptools import find_packages, setup

module_dir = dirname(abspath(__file__))

if __name__ == "__main__":
    setup(
        name="atomate",
        version="1.0.2",
        description="atomate has implementations of FireWorks workflows for Materials Science",
        long_description=open(join(module_dir, "README.md")).read(),
        url="https://github.com/hackingmaterials/atomate",
        author="Anubhav Jain, Kiran Mathew",
        author_email="anubhavster@gmail.com, kmathew@lbl.gov",
        license="modified BSD",
        packages=find_packages(),
        package_data={
            "atomate.vasp.workflows.base": ["library/*"],
            "atomate.vasp.builders": ["*", "examples/*"],
        },
        zip_safe=False,
        install_requires=[
            "custodian>=2019.8.24",
            "FireWorks>=1.4.0",
            "maggma>=0.26.0",
            "monty>=2.0.6",
            "networkx",
            "numpy>=1.20.3,<1.22.0",
            "pandas",
            "paramiko",
            "pydash>=4.1.0",
            "pymatgen-analysis-diffusion>=2021.4.29",
            "pymatgen>=2020.9.14,<=2022.0.10",
            "pymongo<4.0.0",
            "pyyaml>=5.1.2",
            "ruamel.yaml",
            "scipy",
            "tqdm>=4.7.4",
        ],
        extras_require={
            "plotting": ["matplotlib>=1.5.2"],
            "phonons": ["phonopy>=1.10.8"],
            "qchem": ["openbabel"],
            "complete": [
                "matplotlib>=1.5.2",
                "phonopy>=1.10.8",
                "openbabel",
            ],
        },
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Intended Audience :: System Administrators",
            "Intended Audience :: Information Technology",
            "Operating System :: OS Independent",
            "Topic :: Other/Nonlisted Topic",
            "Topic :: Scientific/Engineering",
        ],
        scripts=[join("scripts", f) for f in os.listdir(join(module_dir, "scripts"))],
    )
