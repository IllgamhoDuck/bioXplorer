import os
from setuptools import setup, find_packages


requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
with open(requirements_path) as f:
    requirements = f.read().splitlines()

setup(
    name='BioXplorer',
    version='0.1.0',
    author='HyunByung Park, 42Robotics',
    author_email='hyunbyung87@gmail.com, gdanfx@gmail.com',
    packages=find_packages(),
    install_requires=requirements,
    description='A Python package for exploratory data analysis of cell data using cellxgene.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/IllgamhoDuck/BioXplorer',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='==3.10.*',  # Specifies that your package should only be used with Python 3.10.x
)