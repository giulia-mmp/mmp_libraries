import setuptools
#from mmp_theme import __version__

# Reads the content of your README.md into a variable to be used in the setup below
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='mmp_libraries',
    version='1.0.1',
    packages=['mmp_libraries'],  #, 'oct2pyMMP', 'shieldingMMP'],  # setuptools.find_packages(),
    license='MIT',
    description='Custom libraries and utilities for MMP',
    long_description=long_description,
    long_description_content_type="text/markdown",  
    author='Giulia Bortolato',
    author_email='giulia.bortolato@milanomultiphysics.com',
    url='https://github.com/giulia-mmp/mmp_libraries',
    install_requires=['requests'],
    
    download_url="https://github.com/giulia-mmp/mmp_libraries/archive/refs/tags/v1.0.1.tar.gz"
    
)
