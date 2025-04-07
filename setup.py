import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyMSQ",
    version="0.1.2",
    author="Abdulraheem Musa, Norbert Reinsch",
    author_email="musa@fbn-dummerstorf.de, reinsch@fbn-dummerstorf.de",
    description="A Python package for estimating Mendelian sampling-related quantities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aromemusa/PyMSQ",
    packages=setuptools.find_packages(),
    license="MIT",
    download_url="https://github.com/aromemusa/PyMSQ/dist/PyMSQ-0.1.2.tar.gz",
    keywords=["Mendelian sampling", "variance", "covariance", "similarity", "selection", "haplotype diversity"],
    python_requires=">=3.8",
    install_requires=[
        "numpy<1.25",    # pin to ensure compatibility with older numba versions
        "pandas",
        "scipy",
        "numba<0.58",    # ensure no conflict with numpy <1.25
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    include_package_data=True,
    package_data={"": ["data/*.txt"]},
)
