from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dGinsertion',
      version='1.0a',
      description='deltaG of TM insertion Calculator',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
      ],
      keywords='biology sequence bioinformatics membrane protein',
      url='http://github.com/smsaladi/dGinsertion',
      author='Shyam Saladi',
      author_email='saladi@caltech.edu',
      packages=['dGinsertion'],
      install_requires=['pandas', 'numpy', 'scipy', 'biopython'],
      test_suite='pytest',
      tests_require=['pytest'],
      entry_points={
          'console_scripts': ['calcdG=dGinsertion.calcdG:main'],
      },
      include_package_data=True,
      zip_safe=False)
