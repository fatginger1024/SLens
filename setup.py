import setuptools


setuptools.setup(
    name="SLens",
    version="0.1.1",
    author="Qing Zhou",
    author_email="asmee.rodeus@gmail.com",
    description="A program for strong lensing simulation of 10^5 lenses.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fatginger1024/SLens",
    download_url="https://github.com/fmilthaler/FinQuant/archive/v{}.tar.gz".format(
        ver["release"]
    ),
    license="MIT",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 0 - Alpha",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords=[
        "astronomy",
        "cosmology",
        "strong lensing",
        "bayesian inference",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "tqdm",
    ],
)