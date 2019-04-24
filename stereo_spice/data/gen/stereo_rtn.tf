Dynamic Heliospheric Coordinate Frames developed for the NASA STEREO mission

The coordinate frames in this file all have ID values based on the pattern
18ccple, where

	18 = Prefix to put in the allowed 1400000 to 2000000 range
	cc = 03 for geocentric, 10 for heliocentric
	p  = Pole basis: 1=geographic, 2=geomagnetic, 3=ecliptic, 4=solar
	l  = Longitude basis: 1=Earth-Sun, 2=ecliptic, 4=STEREO-A, 5=STEREO-B
			      6=SOHO
	e  = Ecliptic basis: 0=J2000, 1=mean, 2=true

     Author:  William Thompson
	      NASA Goddard Space Flight Center
	      Code 612.1
	      Greenbelt, MD 20771

	      William.T.Thompson.1@gsfc.nasa.gov


History

    Version 1, 23-Feb-2005, WTT, initial release
    Version 2, 20-Jun-2005, WTT, corrected bug
    Version 3, 23-Feb-2005, WTT, added STEREO-A and -B Science Pointing Frames
    Version 4, 02-May-2006, WTT, added Heliocentric Ecliptic RTN frames and
				 Mission Plane frames
    Version 5, 10-Sep-2008, WTT, added SOHOHGRTN frame
    Version 6, 08-Jan-2015, WTT, added ANGLE_SEP_TOL keywords


STEREO Ahead - Heliocentric Radial Tangential Normal (STAHGRTN) Frame

     Definition of the STEREO-A HGRTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The solar rotation axis is the secondary vector: the Z axis is
	      the component of the solar north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STAHGRTN               =  1810440
        FRAME_1810440_NAME           = 'STAHGRTN'
        FRAME_1810440_CLASS          =  5
        FRAME_1810440_CLASS_ID       =  1810440
        FRAME_1810440_CENTER         =  10
        FRAME_1810440_RELATIVE       = 'J2000'
        FRAME_1810440_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810440_FAMILY         = 'TWO-VECTOR'
        FRAME_1810440_PRI_AXIS       = 'X'
        FRAME_1810440_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810440_PRI_OBSERVER   = 'SUN'
        FRAME_1810440_PRI_TARGET     = 'STEREO AHEAD'
        FRAME_1810440_PRI_ABCORR     = 'NONE'
        FRAME_1810440_PRI_FRAME      = 'IAU_SUN'
        FRAME_1810440_SEC_AXIS       = 'Z'
        FRAME_1810440_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1810440_SEC_FRAME      = 'IAU_SUN'
        FRAME_1810440_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1810440_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

STEREO Behind - Heliocentric Radial Tangential Normal (STBHGRTN) Frame

     Definition of the STEREO-B HGRTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The solar rotation axis is the secondary vector: the Z axis is
	      the component of the solar north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STBHGRTN               =  1810450
        FRAME_1810450_NAME           = 'STBHGRTN'
        FRAME_1810450_CLASS          =  5
        FRAME_1810450_CLASS_ID       =  1810450
        FRAME_1810450_CENTER         =  10
        FRAME_1810450_RELATIVE       = 'J2000'
        FRAME_1810450_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810450_FAMILY         = 'TWO-VECTOR'
        FRAME_1810450_PRI_AXIS       = 'X'
        FRAME_1810450_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810450_PRI_OBSERVER   = 'SUN'
        FRAME_1810450_PRI_TARGET     = 'STEREO BEHIND'
        FRAME_1810450_PRI_ABCORR     = 'NONE'
        FRAME_1810450_PRI_FRAME      = 'IAU_SUN'
        FRAME_1810450_SEC_AXIS       = 'Z'
        FRAME_1810450_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1810450_SEC_FRAME      = 'IAU_SUN'
        FRAME_1810450_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1810450_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

STEREO Ahead - Science Pointing Coordinate Frame

     Definition of the STEREO-A SCPNT Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the Sun relative to the spacecraft is the primary
              vector: the X axis points from the spacecraft to Sun center.
 
              The position of the Earth relative to the spacecraft is the
	      secondary vector: the -Z axis points toward the Earth.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STASCPNT               =  1834440
        FRAME_1834440_NAME           = 'STASCPNT'
        FRAME_1834440_CLASS          =  5
        FRAME_1834440_CLASS_ID       =  1834440
        FRAME_1834440_CENTER         =  10
        FRAME_1834440_RELATIVE       = 'J2000'
        FRAME_1834440_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1834440_FAMILY         = 'TWO-VECTOR'
        FRAME_1834440_PRI_AXIS       = 'X'
        FRAME_1834440_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1834440_PRI_OBSERVER   = 'STEREO AHEAD'
        FRAME_1834440_PRI_TARGET     = 'SUN'
        FRAME_1834440_PRI_ABCORR     = 'NONE'
        FRAME_1834440_PRI_FRAME      = 'J2000'
        FRAME_1834440_SEC_AXIS       = '-Z'
        FRAME_1834440_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1834440_SEC_OBSERVER   = 'STEREO AHEAD'
        FRAME_1834440_SEC_TARGET     = 'EARTH'
        FRAME_1834440_SEC_ABCORR     = 'NONE'
        FRAME_1834440_SEC_FRAME      = 'J2000'
	FRAME_1834440_ANGLE_SEP_TOL  = 0.0001

\begintext

STEREO Behind - Science Pointing Coordinate Frame

     Definition of the STEREO-B SCPNT Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the Sun relative to the spacecraft is the primary
              vector: the X axis points from the spacecraft to Sun center.
 
              The position of the Earth relative to the spacecraft is the
	      secondary vector: the -Z axis points toward the Earth.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STBSCPNT               =  1835450
        FRAME_1835450_NAME           = 'STBSCPNT'
        FRAME_1835450_CLASS          =  5
        FRAME_1835450_CLASS_ID       =  1835450
        FRAME_1835450_CENTER         =  10
        FRAME_1835450_RELATIVE       = 'J2000'
        FRAME_1835450_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1835450_FAMILY         = 'TWO-VECTOR'
        FRAME_1835450_PRI_AXIS       = 'X'
        FRAME_1835450_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1835450_PRI_OBSERVER   = 'STEREO BEHIND'
        FRAME_1835450_PRI_TARGET     = 'SUN'
        FRAME_1835450_PRI_ABCORR     = 'NONE'
        FRAME_1835450_PRI_FRAME      = 'J2000'
        FRAME_1835450_SEC_AXIS       = '-Z'
        FRAME_1835450_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1835450_SEC_OBSERVER   = 'STEREO BEHIND'
        FRAME_1835450_SEC_TARGET     = 'EARTH'
        FRAME_1835450_SEC_ABCORR     = 'NONE'
        FRAME_1835450_SEC_FRAME      = 'J2000'
	FRAME_1835450_ANGLE_SEP_TOL  = 0.0001

\begintext

STEREO Ahead - Heliocentric Ecliptic Radial Tangential Normal (STAHERTN) Frame

     Definition of the STEREO-A HERTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The ecliptic axis is the secondary vector: the Z axis is
	      the component of the ecliptic north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STAHERTN               =  1810341
        FRAME_1810341_NAME           = 'STAHERTN'
        FRAME_1810341_CLASS          =  5
        FRAME_1810341_CLASS_ID       =  1810341
        FRAME_1810341_CENTER         =  10
        FRAME_1810341_RELATIVE       = 'J2000'
        FRAME_1810341_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810341_FAMILY         = 'TWO-VECTOR'
        FRAME_1810341_PRI_AXIS       = 'X'
        FRAME_1810341_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810341_PRI_OBSERVER   = 'SUN'
        FRAME_1810341_PRI_TARGET     = 'STEREO AHEAD'
        FRAME_1810341_PRI_ABCORR     = 'NONE'
        FRAME_1810341_PRI_FRAME      = 'IAU_SUN'
        FRAME_1810341_SEC_AXIS       = 'Z'
        FRAME_1810341_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1810341_SEC_FRAME      = 'ECLIPDATE'
        FRAME_1810341_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1810341_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

STEREO Behind - Heliocentric Ecliptic Radial Tangential Normal (STBHERTN) Frame

     Definition of the STEREO-B HERTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The ecliptic axis is the secondary vector: the Z axis is
	      the component of the ecliptic north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STBHERTN               =  1810351
        FRAME_1810351_NAME           = 'STBHERTN'
        FRAME_1810351_CLASS          =  5
        FRAME_1810351_CLASS_ID       =  1810351
        FRAME_1810351_CENTER         =  10
        FRAME_1810351_RELATIVE       = 'J2000'
        FRAME_1810351_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810351_FAMILY         = 'TWO-VECTOR'
        FRAME_1810351_PRI_AXIS       = 'X'
        FRAME_1810351_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810351_PRI_OBSERVER   = 'SUN'
        FRAME_1810351_PRI_TARGET     = 'STEREO BEHIND'
        FRAME_1810351_PRI_ABCORR     = 'NONE'
        FRAME_1810351_PRI_FRAME      = 'IAU_SUN'
        FRAME_1810351_SEC_AXIS       = 'Z'
        FRAME_1810351_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1810351_SEC_FRAME      = 'ECLIPDATE'
        FRAME_1810351_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1810351_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

STEREO Ahead - STEREO Mission Plane (STAPLANE) Frame

     Definition of the STEREO-A Mission Plane Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The position of the other spacecraft relative to the first
	      defines the secondary vector: the -Y axis is the component of the
	      direction to the Behind spacecraft perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STAPLANE               =  1810940
        FRAME_1810940_NAME           = 'STAPLANE'
        FRAME_1810940_CLASS          =  5
        FRAME_1810940_CLASS_ID       =  1810940
        FRAME_1810940_CENTER         =  10
        FRAME_1810940_RELATIVE       = 'J2000'
        FRAME_1810940_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810940_FAMILY         = 'TWO-VECTOR'
        FRAME_1810940_PRI_AXIS       = 'X'
        FRAME_1810940_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810940_PRI_OBSERVER   = 'SUN'
        FRAME_1810940_PRI_TARGET     = 'STEREO AHEAD'
        FRAME_1810940_PRI_ABCORR     = 'NONE'
        FRAME_1810940_PRI_FRAME      = 'J2000'
        FRAME_1810940_SEC_AXIS       = '-Y'
        FRAME_1810940_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810940_SEC_OBSERVER   = 'STEREO AHEAD'
        FRAME_1810940_SEC_TARGET     = 'STEREO BEHIND'
        FRAME_1810940_SEC_ABCORR     = 'NONE'
        FRAME_1810940_SEC_FRAME      = 'J2000'

\begintext

STEREO Behind - STEREO Mission Plane (STBPLANE) Frame

     Definition of the STEREO-B Mission Plane Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The position of the other spacecraft relative to the first
	      defines the secondary vector: the +Y axis is the component of the
	      direction to the Ahead spacecraft perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_STBPLANE               =  1810950
        FRAME_1810950_NAME           = 'STBPLANE'
        FRAME_1810950_CLASS          =  5
        FRAME_1810950_CLASS_ID       =  1810950
        FRAME_1810950_CENTER         =  10
        FRAME_1810950_RELATIVE       = 'J2000'
        FRAME_1810950_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810950_FAMILY         = 'TWO-VECTOR'
        FRAME_1810950_PRI_AXIS       = 'X'
        FRAME_1810950_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810950_PRI_OBSERVER   = 'SUN'
        FRAME_1810950_PRI_TARGET     = 'STEREO BEHIND'
        FRAME_1810950_PRI_ABCORR     = 'NONE'
        FRAME_1810950_PRI_FRAME      = 'J2000'
        FRAME_1810950_SEC_AXIS       = 'Y'
        FRAME_1810950_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810950_SEC_OBSERVER   = 'STEREO BEHIND'
        FRAME_1810950_SEC_TARGET     = 'STEREO AHEAD'
        FRAME_1810950_SEC_ABCORR     = 'NONE'
        FRAME_1810950_SEC_FRAME      = 'J2000'

\begintext

SOHO - Heliocentric Radial Tangential Normal (SOHOHGRTN) Frame

     Definition of the SOHO HGRTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of the spacecraft relative to the Sun is the primary
              vector: the X axis points from the Sun center to the spacecraft.
 
              The solar rotation axis is the secondary vector: the Z axis is
	      the component of the solar north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.

\begindata

        FRAME_SOHOHGRTN              =  1810460
        FRAME_1810460_NAME           = 'SOHOHGRTN'
        FRAME_1810460_CLASS          =  5
        FRAME_1810460_CLASS_ID       =  1810460
        FRAME_1810460_CENTER         =  10
        FRAME_1810460_RELATIVE       = 'J2000'
        FRAME_1810460_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1810460_FAMILY         = 'TWO-VECTOR'
        FRAME_1810460_PRI_AXIS       = 'X'
        FRAME_1810460_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1810460_PRI_OBSERVER   = 'SUN'
        FRAME_1810460_PRI_TARGET     = 'SOHO'
        FRAME_1810460_PRI_ABCORR     = 'NONE'
        FRAME_1810460_PRI_FRAME      = 'IAU_SUN'
        FRAME_1810460_SEC_AXIS       = 'Z'
        FRAME_1810460_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1810460_SEC_FRAME      = 'IAU_SUN'
        FRAME_1810460_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1810460_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

