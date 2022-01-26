from setuptools import setup, find_packages, Extension

classifiers = [
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: GNU General Public License V2 (GPLv2)',
    'Natural Language :: English',
    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: Implementation :: CPython',
]

ext_modules = [
    Extension(
        'mykmeanssp',
        ['kmeans.c'],
    ),
]

setup(
    name = 'mykmeanssp',
    version ='0.1.0',
    description="The k-means algorithm using given initial centroids",
    install_requires=['invoke'],
    packages=find_packages(),
    license='GPL-2',
    classifiers=classifiers,
    ext_modules=ext_modules
)