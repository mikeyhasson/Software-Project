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
        'capi_spk',
        ['spkmeansmodule.c', 'error_msg.c', 'file_translators.c',
                'format_printer.c', 'kmeans.c',
                'matrix.c', 'nsc.c', 'spkmeans.c',
                'vector.c'],
    ),
]

setup(
    name='capi_spk',
    version ='0.1.0',
    description="A Full Spectral Clustering Algorithm Interface",
    install_requires=['invoke'],
    packages=find_packages(),
    license='GPL-2',
    classifiers=classifiers,
    ext_modules=ext_modules
)