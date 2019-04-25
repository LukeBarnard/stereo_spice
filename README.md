# stereo_spice
--------

This module computes state vectors for the STEREO spacecraft and other solar system bodies, with a syntax and functionality very similar to (and largely copied from) that provided by SSWIDL. It essentially ports the old SSWIDL code to python, wrapping around spiceypy and the SPICE kernals, which do the actual computations. 

This is based of a class, StereoSpice, which on initialisation will load the relevant SPICE kernals into memory, and provide the functions required to work out the state vectors in a range of coordinate systems, and transform between them. Each instantiation of the StereoSpice class reloads these kernals into memory, which can cause a memory error if many seperate instances are made. So care should be taken not to create instances of StereoSpice inside loops ect. The kernals can be cleared from memory with StereoSpice.clear_kernals().



Example:

~~~
from astropy.time import Time
from stereo_spice.coordinates import StereoSpice
import astropy.units as u
spice = StereoSpice()

time = Time('2010-01-01T12:00:00')
print "Date: {}".format(time.isot)
for system in ['HCI', 'CARRINGTON', 'HEEQ', 'HEE', 'HAE']:
    stb = spice.get_lonlat(time, 'stb', system, precess=False)
    ert = spice.get_lonlat(time, 'earth', system, precess=False)
    sta = spice.get_lonlat(time, 'sta', system, precess=False)
    print "{0} Longitude>> stb:{1:5.3f}  ert:{2:5.3f}  sta:{3:5.3f}".format(system, stb[1], ert[1], sta[1])
    print "{0} Latitude>> stb:{1:5.3f}  ert:{2:5.3f}  sta:{3:5.3f}".format(system, stb[2], ert[2], sta[2])
~~~
    