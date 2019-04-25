from setuptools import setup
import os
import fnmatch

def data_file_list():
    """
    Construct a list of all the data files, including the vey many SPICE kernals.
    """
    root_path = os.path.split(os.path.abspath(__file__))
    path = os.path.join(root_path, 'stereo_spice')
    data_files = []
    for root, dirnames, filenames in os.walk(path):
        for filename in filenames:
            full_path = os.path.join(root, filename)
            rel_path = os.path.relpath(full_path, path)
            if not rel_path.endswith('py'):
                data_files.append(rel_path)
    


def readme():
    """
    Function to parse the readme file into setup.
    """
    with open('README.md') as f:
        return f.read()

setup(name='stereo_spice',
      version='0.1',
      description='Compute state vectors of solar system bodies and the STEREO spacecraft, with functionality similar to that of SSWIDL. This is essentially a small wrapper around spiceypy and STEREO spice kernals.',
      long_description=readme(),
      long_description_content_type="text/markdown",
      license='MIT',
      url='https://github.com/LukeBarnard/stereo_spice.git',
      author='Luke Barnard',
      author_email='l.a.barnard@reading.ac.uk',
      packages=['stereo_spice'],
      package_data={'stereo_spice': data_file_list()},
      install_requires=['astropy', 'spiceypy', 'numpy'],
      include_package_data=True,
      zip_safe=False)