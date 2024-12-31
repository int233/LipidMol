from setuptools import setup, find_packages

setup(
    name="LipidMol",
    version="0.1.1",
    description="Calculate the exact mass of lipids based on their common names",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="int0030",
    author_email="int0030@163.com",
    url="https://github.com/int233/LipidMol",
    packages=find_packages(),
    install_requires=[
        "molmass",
        "loguru"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering",
        "Development Status :: 4 - Beta",
    ],
    python_requires='>=3.8',
    include_package_data=False,
)
