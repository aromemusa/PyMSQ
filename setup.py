from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="PyMSQ",
    version="0.1.3",
    author="Abdulraheem Musa, Norbert Reinsch",
    author_email="musa@fbn-dummerstorf.de, reinsch@fbn-dummerstorf.de",
    description="A Python package for fast Mendelian sampling (co)variance and haplotype-based similarity in genomic selection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aromemusa/PyMSQ",

    project_urls={
        "Source": "https://github.com/aromemusa/PyMSQ",
        "Bug Tracker": "https://github.com/aromemusa/PyMSQ/issues",
        "Documentation": "https://github.com/aromemusa/PyMSQ/tree/main/docs",
        "Software paper": "https://doi.org/10.1186/s12859-026-06392-5",
        "Method paper": "https://doi.org/10.1111/jbg.12930",
        "Zenodo DOI": "https://doi.org/10.5281/zenodo.18760267",
    },

    packages=find_packages(include=["PyMSQ", "PyMSQ.*"]),
    license="MIT",
    keywords=[
        "Mendelian sampling",
        "variance",
        "covariance",
        "similarity",
        "selection",
        "haplotype diversity",
        "genomic selection",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy<1.25",
        "pandas",
        "scipy",
        "numba<0.58",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    include_package_data=True,
    package_data={"PyMSQ": ["data/*.txt"]},
)