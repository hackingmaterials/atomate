import os
from os.path import abspath, dirname, join

from setuptools import find_packages, setup

module_dir = dirname(abspath(__file__))

if __name__ == "__main__":
    setup(
        name="atomate",
        version="1.0.3",
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
            "maggma>=0.39.0",
            "monty>=2.0.6",
            "networkx",
            "numpy",
            "pandas",
            "paramiko",
            "pydash>=4.1.0",
            "pymatgen-analysis-diffusion",
            "pymatgen",
            "pymongo",
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
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
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
