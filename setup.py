from setuptools import setup, find_packages
import glob
import os
from casm_figures import __version__


# get console_scripts
def script_str(file):
    name = os.path.splitext(os.path.split(file)[1])[0]
    if name in ['casm_calc', 'casm_learn']:
        return name.replace('_', '-') + '=casm.scripts.' + name + ':main'
    else:
        return name.replace('_', '.') + '=casm.scripts.' + name + ':main'


console_scripts = [
    script_str(x) for x in glob.glob('casm/scripts/*.py')
    if x != 'casm/scripts/__init__.py'
]
print(console_scripts)

description = "Useful figure plotting for CASM"
long_description = "Useful figure plotting for CASM"
# with open(os.path.join('..', 'README.md'), encoding='utf-8') as f:
#     long_description = f.read()

setup(
    name='casm_figures',
    version=__version__,
    url='https://github.com/bpuchala/casm_figures',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Brian Puchala',
    author_email='bpuchala@umich.edu',
    license='MIT',
    packages=find_packages(),
    entry_points={'console_scripts': console_scripts},
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'scipy',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering'
    ],
    data_files=[('', [])])
