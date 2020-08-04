

# stereo_spice
--------

This module computes state vectors for the STEREO spacecraft and other solar system bodies, with a syntax and functionality very similar to (and largely copied from) that provided by SSWIDL. It essentially ports the old SSWIDL code to python, wrapping around spiceypy and the SPICE kernals, which do the actual computations. 

This is based on a class, StereoSpice, which on initialisation will load the relevant SPICE kernals into memory, and provide the functions required to work out the state vectors in a range of coordinate systems, and transform between them. Each instantiation of the StereoSpice class reloads these kernals into memory, which can cause a memory error if many seperate instances are made. So care should be taken not to create instances of StereoSpice inside loops ect. The kernals can be cleared from memory with StereoSpice.clear_kernals().

## Installation

``stereo_spice`` is written in Python 3.7.3, and depends upon ``numpy``, ``astropy`` and ``spiceypy``. 

The easiest way to install ``stereo_spice`` is as:
~~~
pip install https://github.com/LukeBarnard/stereo_spice/archive/master.zip
~~~

This command can also be used to install ``stereo_spice`` within a conda environment. For example, 
~~~
conda create -n myenv python=3.7.3 astropy=4.0 spiceypy=2.3.0 pip
conda activate myenv
pip install https://github.com/LukeBarnard/stereo_spice/archive/master.zip
~~~
Following these commands, ``stereo_spice`` will be available from within ``myenv``.

## Example use
The below script prints out the latitude and longitude of Earth, STEREO-A and STEREO-B, for a range of coordinate systems on 2010-01-01T12:00. Comparison with the [STEREO orbit tool](https://stereo-ssc.nascom.nasa.gov/where/) shows that the coordinates are correct. 
~~~
from astropy.time import Time
from stereo_spice.coordinates import StereoSpice
spice = StereoSpice()

time = Time('2010-01-01T12:00:00')
print("Date: {}".format(time.isot))
for system in ['HCI', 'CARRINGTON', 'HEEQ', 'HEE', 'HAE']:
    stb = spice.get_lonlat(time, 'stb', system, precess=False)
    ert = spice.get_lonlat(time, 'earth', system, precess=False)
    sta = spice.get_lonlat(time, 'sta', system, precess=False)
    print("{0} Longitude>> stb:{1:5.3f}  ert:{2:5.3f}  sta:{3:5.3f}".format(system, stb[1], ert[1], sta[1]))
    print("{0} Latitude>> stb:{1:5.3f}  ert:{2:5.3f}  sta:{3:5.3f}".format(system, stb[2], ert[2], sta[2]))
~~~
## Contact
Please contact [Luke Barnard](https://github.com/lukebarnard)
    
