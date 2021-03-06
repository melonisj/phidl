#!/usr/bin/env python

from setuptools import setup

install_requires=[
   'gdspy>=1.2',
   'numpy',
   'matplotlib',
   'pyyaml',
   'scikit-image>=0.11',
   'webcolors',
]

setup(name='phidl',
      version='0.8.1',
      description='PHIDL',
      install_requires=install_requires,
      author='Adam McCaughan',
      author_email='amccaugh@gmail.com',
      packages=['phidl'],
      py_modules=['phidl.geometry', 'phidl.routing'],
      package_dir = {'phidl': 'phidl'},
     )