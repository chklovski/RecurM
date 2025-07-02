from setuptools import setup, find_packages

setup(
    name='RecurM',
    version='0.3.2',
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/chklovski/RecurM',
    license='',
    install_requires=(),
    author='Alex Chklovski',
    scripts=['bin/recurm'],
    author_email='chklovski@gmail.com',
    description='RecurM - recovering MGEs from metagenomic data'
)