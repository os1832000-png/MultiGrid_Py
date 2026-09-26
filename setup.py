from setuptools import setup, find_packages

setup(
    name="marfa",
    version="3.0.0",
    description="MARFA: Molecular atmospheric Absorption with Rapid and Flexible Analysis",
    url="https://github.com/yourname/MARFA_PY_SIMPLE",
    author="Osama",
    author_email="osama18@example.com",
    license="MIT",
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["marfa", "marfa.*"]),
    python_requires=">=3.6",
    install_requires=["numpy>=1.16", "scipy>=1.2"],
    extras_require={"plot": ["matplotlib>=3.0"], "dev": ["pytest"]},
    entry_points={"console_scripts": ["marfa=marfa.cli:main"]},
)
