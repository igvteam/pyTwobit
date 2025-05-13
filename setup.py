from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pytwobit',  # Normalized project name
    version='0.3.1',
    description='A fast reader for local or remote UCSC twobit sequence files.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    author='Jim Robinson',
    url='https://github.com/jrobinso/pytwobit',
    packages=['pytwobit'],  # Ensure this matches the normalized name
    package_data={'pytwobit': ['tests/foo.2bit']},  # Match normalized name
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