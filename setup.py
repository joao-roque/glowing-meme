import os

from setuptools import find_packages, setup
import re

REGEX_COMMENT = re.compile(r'[\s^]#(.*)')

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

dir_path = os.path.dirname(os.path.realpath(__file__))


def parse_requirements(filename):
    with open(filename, 'rt') as filehandle:
        return tuple(filter(None, (
            REGEX_COMMENT.sub('', line).strip()
            for line in filehandle
        )))

setup(
    name='glowingmeme',
    packages=find_packages(),
    include_package_data=True,
    author='Joao Almeida',
    author_email='joao_roque.eu@hotmail.com',
    classifiers=[
        'Environment :: Other Environment',
        'Intended Audience :: Other Audience',
        'License :: Other/Proprietary License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
    ],
    install_requires=parse_requirements('requirements.txt')
)
