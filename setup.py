from setuptools import setup, find_packages

setup(
    name="multigrid_py",
    version="3.0.2",
    description="Python port of the Fortran MARFA line-by-line molecular absorption code",
    url="https://github.com/os1832000-png/MultiGrid_Py",
    license="MIT",
    author="Osama M. M. Abdellatif",
    author_email="osama1832000@gmail.com",
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["multigrid", "multigrid.*"]),
    python_requires=">=3.6",
    install_requires=["numpy>=1.16", "scipy>=1.2"],
    extras_require={"plot": ["matplotlib>=3.0"], "dev": ["pytest"]},
    entry_points={"console_scripts": ["multigrid=multigrid.cli:main"]},
)
