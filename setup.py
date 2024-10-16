from setuptools import setup, find_packages

setup(
    name="svtoolbox",
    version="v1.1",
    packages=find_packages("src"),
    package_dir={"": "src"},
    test_suite="tests",
    entry_points={"console_scripts": ["svtoolbox = svtoolbox.client:run"]},
    python_requires=">=3.10",
    install_requires=["click", "pysam", "setuptools"],
    author="Michael Knudsen",
    author_email="micknudsen@gmail.com",
)
