from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='fragforce',
      version='0.1.1',
      description='Simulates the motion and surface forces of a fluid ellipsoid in shear flow.',
      url='http://github.com/MathBioCU/fragforce',
      author='Eric Kightley',
      author_email='kightley.1@gmail.com',
      license='MIT',
      packages=['fragforce'],
      install_requires=[
          'numpy',
          'scipy',
      ],
      include_package_data=True,
      zip_safe=False)
