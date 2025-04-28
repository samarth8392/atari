from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="atari",
    version="0.1.0",
    author="Samarth Mathur",
    author_email="samarth8392@gmail.com",
    description="AlternaTe Allele Read vIsualizer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samarth8392/atari",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
    install_requires=[
        "pysam",
        "pandas",
        "numpy",
        "matplotlib",
        "tqdm",
        "rich",
    ],
    entry_points={
        'console_scripts': [
            'atari=atari.cli:main',
        ],
    },
)