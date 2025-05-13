
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyTwobit',
    version='0.3.0',
    description='A fast reader for local or remote UCSC twobit sequence files.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    author='Jim Robinson',
    url='https://github.com/jrobinso/pyTwobit',
    packages=['pytwobit'],
    package_data={'pytwobit': ['tests/foo.2bit']},
    install_requires=[
        'requests',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
