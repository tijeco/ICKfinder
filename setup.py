from setuptools import setup


setup(name='cluck',
      version='0.1.1',
      description='The cluck joke in the world',
      url='http://github.com/tijeco/cluck',
      author='Flying Circus',
      author_email='flyingcircus@example.com',
      entry_points = {'console_scripts': ['cluck=src.main:main'],},
      license='MIT',
      packages=['src'],
      include_package_data=True,
      zip_safe=False)
