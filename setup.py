import setuptools

with open("README.md", "r", newline="", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="PyMSQ",
    packages=setuptools.find_packages(),
    version="0.1.1",
    license="MIT",
    description="A Python package for estimating Mendelian sampling-related quantities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Abdulraheem Musa, Norbert Reinsch",
    author_email="musa@fbn-dummerstorf.de, reinsch@fbn-dummerstorf.de",
    url="https://github.com/aromemusa/PyMSQ",
    download_url="https://github.com/aromemusa/PyMSQ/blob/main/dist/PyMSQ-0.1.tar.gz",
    keywords=["Mendelian sampling", "variance", "covariance", "similarity", "selection", "haplotype diversity"],
    install_requires=[
        "numpy<1.25",
        "pandas",
        "scipy",
        "numba<0.58",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Education",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    include_package_data=True,
    package_data={"": ["data/*.txt"]},
)
