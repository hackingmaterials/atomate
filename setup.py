import os
from os.path import abspath, dirname, join

from setuptools import find_packages, setup

module_dir = dirname(abspath(__file__))

if __name__ == "__main__":
    setup(
        name="atomate",
        version="1.1.0",
        description="atomate has implementations of FireWorks workflows for Materials Science",
        long_description=open(join(module_dir, "README.md")).read(),
        long_description_content_type="text/markdown",
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
            "custodian>=2023.7.22",
            "FireWorks>=2.0.3",
            "maggma>=0.51.25",
            "monty>=2023.5.8",
            "networkx",
            "numpy",
            "pandas",
            "paramiko",
            "pydash>=7.0.6",
            "pymatgen-analysis-diffusion>=2023.8.15",
            "pymatgen-analysis-defects>=2023.7.24",
            "pymatgen>=2023.7.20",
            "pymongo",
            "pyyaml>=5.1.2",
            "ruamel.yaml",
            "scipy",
            "tqdm>=4.65.0",
        ],
        extras_require={
            "plotting": ["matplotlib>=1.5.2"],
            "phonons": ["phonopy>=1.10.8"],
            "qchem": ["openbabel-wheel"],
            "defects": ["pymatgen-analysis-defects"],
            "complete": [
                "matplotlib>=1.5.2",
                "phonopy>=1.10.8",
                "openbabel-wheel",
                "boto3>=1.28.15",
                "Flask>=2.3.2",
                "coverage>=7.2.7",
                "moto>=4.1.14",
                "pytest-cov>=4.1.0",
                "pytest>=7.4.0",
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
