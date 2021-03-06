import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read(23)
    
    
setuptools.setup(
    name="SLens",
    version="0.1.1",
    author="Qing Zhou",
    author_email="asmee.rodeus@gmail.com",
    description="A program for strong lensing simulation of 10^5 lenses.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fatginger1024/SLens",
    license="MIT",
    packages=setuptools.find_packages(),
    package_data={'test_data': ['*']},
    include_package_data=True,
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
