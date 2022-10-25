from setuptools import setup, find_packages
import os

# first setuptools==58

install_requires = [
    "varcode==1.0.3",
    "PyVCF3==1.0.3",
    "pandas==1.4.3",
    "numpy==1.22.3",
    "biopython==1.79"
]

setup(
    include_package_data=False,
    name='clinparse',
    version='0.1',
    description='get missense pathogenic/benign variants in AminoAcid positions',
    author='chang li',
    license='MIT',
    packages=find_packages(),
    zip_safe=False,
    install_requires=install_requires,
)

os.system("pyensembl install --release 106 --species human")