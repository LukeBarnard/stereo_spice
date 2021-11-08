The files in this directory tree are reduced resolution versions of the STEREO
Attitude History files.  They are provided as a service for those who do not
need the full resolution versions.

To reduce the size of the Attitude History files, the data are first passed
through a 5 minute smoothing filter.  The data are then subsampled to
sufficient cadence to keep the interpolated pointing within 2 arcseconds of the
smoothed pointing.  Note that the 5 minute smoothing operation has a
significant effect during maneuvers--these data are not appropriate for
detailed examination of the pointing during maneuvers.

To further restrict the total data size, files no longer in use will be removed
from this directory tree.  (The full resolution versions will be preserved.)

The original full resolution attitude history files are available as part of
the SolarSoft Database tree, in the directory $SSWDB/stereo/gen/spice/ah.  See

	     http://www.lmsal.com/solarsoft/sswdb_description.html

for more information on how to install the optional SolarSoft Database.
