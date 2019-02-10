# coding: utf-8

from setuptools import setup, find_packages

NAME = "bwb-scheduling-service"
VERSION = "1.0.0"

# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

REQUIRES = ["connexion"]

setup(
    name=NAME,
    version=VERSION,
    description="BWB Scheduler",
    author_email="",
    url="",
    keywords=["Swagger", "BWB Scheduler"],
    install_requires=REQUIRES,
    packages=find_packages(),
    package_data={'': ['swagger/swagger.yaml']},
    include_package_data=True,
    entry_points={
        'console_scripts': ['swagger_server=swagger_server.__main__:main']},
    long_description="""\
    Docker container scheduler for BWB
    """
)
