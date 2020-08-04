import os
import glob
import numpy as np
import spiceypy as spice
from astropy.time import Time


class StereoSpice:

    def __init__(self):
        """
        Load in the general spice kernals and stereo spice keranls
        """
        self.__load_general_kernals__()

        self.__load_stereo_kernals__()
        return
        

    def __get_spice_range__(self, filename):
        """
        Function to calculate the range of coverage of a spice file given by filename
        :param filename: String, full path to spice file.
        :return dates_min: Astropy time object giving start of period spanned by spice file
        :return dates_max: Astropy time object giving end of period spanned by spice file
        :return craft_ids: List of craft ids represented by the spice file
        """

        # Ephemeris files:
        if filename.endswith('bsp'):
            # Get craft id's
            craft_ids = spice.spkobj(filename)
            times = []
            for s in craft_ids:
                cover = spice.utils.support_types.SPICEDOUBLE_CELL(2000)
                spice.spkcov(filename, s, cover)
                times.append([c for c in cover])
        # Pointing files
        elif filename.endswith('bc'):
            # Get craft id's
            craft_ids = spice.ckobj(filename)
            times = []
            for s in craft_ids:
                cover = spice.utils.support_types.SPICEDOUBLE_CELL(2000)
                try:
                    print('spice.ckcov: compute segment')
                    spice.ckcov(filename, s, False, 'segment', 0.0, 'TDB', cover)
                except:
                    print('spice.ckcov: compute interval')
                    spice.ckcov(filename, s, False, 'interval', 0.0, 'TDB', cover)

                times.append([c for c in cover])
        else:
            print('Unrecognized file extension : ' + filename.split('.')[-1])
            dates_min = np.NaN
            dates_max = np.NaN
            craft_ids = np.NaN
            return dates_min, dates_max, craft_ids

        # Format the dates.
        min_time = min([min(t) for t in times])
        dates_min = Time(spice.et2utc(min_time, 'ISOC', 3))
        max_time = max([max(t) for t in times])
        dates_max = Time(spice.et2utc(max_time, 'ISOC', 3))
        return dates_min, dates_max, craft_ids
        

    def __get_kernal_files__(self):
        """
        Function to load in the paths to the relevant spice kernals and files with lists of spice kernals.
        :return kernals_dict: Dictionary containing the full paths to the kernal files needed.
        """

        root  = os.path.abspath(os.path.dirname(__file__))
        kernal_files = os.path.join(root,'config.dat')
        
        if os.path.exists(kernal_files):

            with open(kernal_files, "r") as f:
                kernals_dict = {}
                lines = f.read().splitlines()
                for l in lines:
                    name, rel_path = l.split(';')
                    abs_path = os.path.join(root, rel_path)
                    kernals_dict[name] = abs_path

            # Check that each file exists.
            for k, f in kernals_dict.items():
                
                if not os.path.exists(f):
                    print("Error: {} kernal file/file-list not found".format(k))
        else:
            print("Error: Kernal files list not found.")
            kernals_dict = {}

        return kernals_dict
        

    def __load_general_kernals__(self):
        """
        Function to load in all the general kernals.
        :return:
        """

        # Load in the necessary Spice kernels
        kernals_dict = self.__get_kernal_files__()
        
        # Load in the leap seconds kernal
        spice.furnsh(kernals_dict['leap_seconds'])
        self.__LeapSec__ = kernals_dict['leap_seconds']

        # Load in the solar system kernal
        spice.furnsh(kernals_dict['solar_system'])
        self.__SolarSystem__ = kernals_dict['solar_system']

        # Load in the planetary constants kernal
        spice.furnsh(kernals_dict['planet_constants'])
        self.__PlanetConstants__ = kernals_dict['planet_constants']

        # Load in the heliospheric frames kernal
        spice.furnsh(kernals_dict['heliospheric_frames'])
        self.__HeliosphericFrames__ = kernals_dict['heliospheric_frames']

        # Load in the stereo frames kernal
        spice.furnsh(kernals_dict['stereo_frames'])
        self.__StereoFrames__ = kernals_dict['stereo_frames']

        # Load in the stereo a clock kernal
        spice.furnsh(kernals_dict['stereo_frames'])
        self.__staclock__ = kernals_dict['sta_clock']

        # Load in the stereo b clock kernal
        spice.furnsh(kernals_dict['stereo_frames'])
        self.__stbclock__ = kernals_dict['stb_clock']

        return
        

    def __load_stereo_kernals__(self):
        """
        Function to load in the STEREO specific kernals
        """
        # Load in the necessary Spice kernels
        kernals_dict = self.__get_kernal_files__()
        
        # Load in the predicted ephemeris, and collect all ephem files.
        # >> Initialise record of maximum dates covered by the ephemeris
        self.__PredMaxdates__ = {'sta': Time(0.0, format='jd'),
                                'stb': Time(0.0, format='jd')}

        # >> Initialise the conic to cover the predicted ephemeris outside covered range
        self.__PredConic__ = {'sta': None, 'stb': None}
        self.__mu__ = 1.32712440018
        self.__AllPredictedEphem__ = []

        for craft in ['sta', 'stb']:
            # Load in the predicted ephemeris files
            key = craft + "_ephem_filelist"
            with open(kernals_dict[key]) as f:
                all_ephem_files = f.read().splitlines()
                for ephem_file in all_ephem_files:
                    full_ephem_path = os.path.join(kernals_dict['data_root'], ephem_file)
                    if os.path.exists(full_ephem_path):
                        spice.furnsh(full_ephem_path)
                        self.__AllPredictedEphem__.append(full_ephem_path)

                        # Extract the orbital parameters for conic extrapolation on predicted ephemeris
                        dates_min, dates_max, craft_ids = self.__get_spice_range__(full_ephem_path)
                        if dates_max > self.__PredMaxdates__[craft]:
                            self.__PredMaxdates__[craft] = dates_max
                            state = self.get_coord(dates_max, craft, system='HAE')
                            et = spice.utc2et(dates_max.isot)
                            elts = spice.oscelt(state, et, self.__mu__)
                            self.__PredConic__[craft] = elts
                    else:
                        print("Error: File not found - {}".format(full_ephem_path))

        # Load in the definitive ephemeris, and collect all ephem files.
        # >> Initialise record of maximum dates covered by the ephemeris
        self.__DefMaxdates__ = {'sta': Time(0.0, format='jd'),
                               'stb': Time(0.0, format='jd')}

        # >> Initialise the conic to cover the predicted ephemeris outside covered range
        self.__DefConic__ = {'sta': None, 'stb': None}
        self.__AllDefEphem__ = []

        for craft in ['sta', 'stb']:
            # Load in the predicted ephemeris files
            key = craft + "_def_ephem_filelist"
            with open(kernals_dict[key]) as f:
                all_ephem_files = f.read().splitlines()
                for ephem_file in all_ephem_files:
                    full_ephem_path = os.path.join(kernals_dict['data_root'], ephem_file)
                    if os.path.exists(full_ephem_path):
                        spice.furnsh(full_ephem_path)
                        self.__AllDefEphem__.append(full_ephem_path)

                        # Extract the orbital parameters for conic extrapolation on definitive ephemers
                        dates_min, dates_max, sid = self.__get_spice_range__(full_ephem_path)
                        if dates_max > self.__DefMaxdates__[craft]:
                            self.__DefMaxdates__[craft] = dates_max
                            state = self.get_coord(dates_max, craft, system='HAE')
                            et = spice.utc2et(dates_max.isot)
                            elts = spice.oscelt(state, et, self.__mu__)
                            self.__DefConic__[craft] = elts
                    else:
                        print("Error: File not found - {}".format(full_ephem_path))

        # TODO: Load in Attitude kernals if want to use CMAT stuff?

        return
    
    
    def convert_hpc_to_hpr(self, hpc_lon, hpc_lat, degrees=True):
        """
        Function to convert helioprojective cartesian coordinates (longitudes and latitudes) into helioprojective radial
        coordinates (elongations and position angles). Conversion done by Eqn. 19 in Thompson 2006.
        :param hpc_lon: Float or array of longitudes, in degrees.
        :param hpc_lat: Float or array of latitudes, in degrees.
        :param degrees: Boolean, if True (default), angles are parsed and returned in degrees.
        :return hpr_el: Float or array of elongations, in degrees.
        :return hpr_pa: Float or array of position angles, in degrees.
        """
        
        if degrees:
            hpc_lon = np.deg2rad(hpc_lon)
            hpc_lat = np.deg2rad(hpc_lat)
            
        # Elongation calc:
        # Get numerator and denomenator for atan2 calculation
        btm = np.cos(hpc_lat) * np.cos(hpc_lon)
        top = np.sqrt((np.cos(hpc_lat) ** 2) * (np.sin(hpc_lon) ** 2) + (np.sin(hpc_lat) ** 2))
        hpr_el = np.arctan2(top, btm)
        # Position angle calc:
        btm = np.sin(hpc_lat)
        top = -np.cos(hpc_lat) * np.sin(hpc_lon)
        hpr_pa = np.arctan2(top, btm)
        # Correct eastern longitudes so pa runs from 0>2pi, rather than 0>pi.
        if isinstance(hpr_pa, np.float):
            if hpc_lon >= 0:
                hpr_pa += 2 * np.pi
        else:
            hpr_pa[hpc_lon >= 0] += 2 * np.pi

        if degrees:
            # Put it back into degs
            hpr_el = np.rad2deg(hpr_el)
            hpr_pa = np.rad2deg(hpr_pa)
            
        return hpr_el, hpr_pa
        
        
    def convert_hpr_to_hpc(self, hpr_el, hpr_pa, degrees=True):
        """
        Function to convert helioprojective radial coordinates (elongations and position angles) into helioprojective
        cartesian coordinates (longitudes and latitudes) . Conversion done by Eqn. 20 in Thompson 2006.
        :param hpr_el: Array of elongations. Should have astropy unit of degrees.
        :param hpr_pa: Array of position angles. Should have astropy unit of degrees.
        :param degrees: Boolean, if True (default), angles are parsed and returned in degrees.
        :return hpc_lon: Array of longitudes with astropy unit of degrees.
        :return hpc_lat: Array of latitudes angles with astropy unit of degrees.
        """
        if degrees:
            hpr_el = np.deg2rad(hpr_el)
            hpr_pa = np.deg2rad(hpr_pa)
        
        # Longitude calc:
        # Get numerator and denomenator for atan2 calculation
        btm = np.cos(hpr_el)
        top = -np.sin(hpr_el) * np.sin(hpr_pa)
        hpc_lon = np.arctan2(top, btm)
        # Latitude calc:
        hpc_lat = np.arcsin(np.sin(hpr_el) * np.cos(hpr_pa))
        
        if degrees:
            hpc_lon = np.rad2deg(hpc_lon)
            hpc_lat = np.rad2deg(hpc_lat)
            
        return hpc_lon, hpc_lat

        
    def convert_hpc_to_rtn(self, hpc_lon, hpc_lat, degrees=True):
        """
        Function to convert helioprojective cartesian coordinates (longitudes and latitudes) into RTN longitude and latitudes.
        :param hpc_lon: Float or array of HPC longitudes, in degrees.
        :param hpc_lat: Float or array of HPC latitudes, in degrees.
        :param degrees: Boolean, if True (default), angles are parsed and returned in degrees.
        :return rtn_lon: Float or array of RTN longitudes, in degrees.
        :return rtn_lat: Float or array of RTN latitudes, in degrees.
        """
        rtn_lat = hpc_lat
        
        if degrees:
            hpc_lon = np.deg2rad(hpc_lon)
        
        rtn_lon = np.pi - hpc_lon
        
        if degrees:
            rtn_lon = np.rad2deg(rtn_lon)
            
        return rtn_lon, rtn_lat
    
    
    def convert_rtn_to_hpc(self, rtn_lon, rtn_lat, degrees=True):
        """
        Function to convert helioprojective cartesian coordinates (longitudes and latitudes) into RTN longitude and latitudes.
        :param hpc_lon: Float or array of HPC longitudes, in degrees.
        :param hpc_lat: Float or array of HPC latitudes, in degrees.
        :param degrees: Boolean, if True (default), angles are parsed and returned in degrees.
        :return rtn_lon: Float or array of RTN longitudes, in degrees.
        :return rtn_lat: Float or array of RTN latitudes, in degrees.
        """
        hpc_lat = rtn_lat
        
        if degrees:
            rtn_lon = np.deg2rad(rtn_lon)
            
        if isinstance(rtn_lon,float):
            hpc_lon = np.pi - rtn_lon
            if hpc_lon > np.pi:
                hpc_lon -= 2*np.pi

            if hpc_lon < -np.pi:
                hpc_lon += 2*np.pi
        else:
            hpc_lon = np.pi - rtn_lon
            id_over = hpc_lon > np.pi
            if any(id_over):
                hpc_lon[id_over] -= 2*np.pi

            id_under = hpc_lon < -np.pi
            if any(id_under):
                hpc[id_under] += 2*np.pi
        
        if degrees:
            hpc_lon = np.rad2deg(hpc_lon)
        
        return hpc_lon, hpc_lat

        
    def convert_coord(self, dates, coord_src, system_src, system_dst, observe_src=None, observe_dst=None, precess=False):
        """
        Function to convert coordinates betwen different reference frames.
        :param dates: Astropy time object of dates(s).
        :param coord_src: Array of coordinates to convert. Should be len(dates)*3 for positions only, or len(dates)*6 for full state
        :param system_src: String name of coordinate system of coord array.
        :param system_dst: String name of coordinate system to transform coord array to.
        :param observe_src: String name of observatory for origin of system from. Only needed for some systems.
        :param observe_dst: String name of observatory for origin of system to. Only needed for some systems.
        :param precess: Boolean. If True accounts for precession in coordinate system.
        :return state: Array, giving state at each dates. Either len(dates)x3 or len(dates)x6, depending on no_velocity
        :return ltime (optional): The light travel time between observatory and target
        """
        # If coordinates input as a list, then bung them into an array.
        if isinstance(coord_src, list):
            if all([isinstance(c, (float, int)) for c in coord_src]):
                coord_src = np.array(coord_src)
                coord_src = np.squeeze(coord_src)
            else:
                print("ERROR: coord_src should be a numpy array of coordinates or a list of floats/ints.")
            
        # If coord only has one dimension, set it so that time is zeroth.
        if coord_src.ndim == 1:
            n_coords = 1
            n_components = coord_src.size
        else:
            n_coords = coord_src.shape[0]
            n_components = coord_src.shape[1]
            
        # Check dates and coord sizes match.
        if dates.size != n_coords:
            print("Error: Number of dates does not correspond to number of coordinates.")
            
        # Check coord components are either 3 or 6, for position or state. Set no_velocity flag too.
        if n_components == 3:
            no_velocity = True
        elif n_components == 6:
            no_velocity = False
        else:    
            print("ERROR: Invalid dimension of position vector or state vector")
    
        # Get NAIF formated name for src and dst observer and frame.
        # Get SRC frame and osberver
        if observe_src is None:
            # Observer defined by the frame.
            frame_src, observe_src = self.get_system_frame_names(system_src, precess=precess)
        else:
            # Observer must be specified for this frame
            observe_src = self.get_naif_body_code(observe_src)
            # Get naif frame and observer
            frame_src, observe_src = self.get_system_frame_names(system_src, observatory=observe_src, precess=precess)
            
        # Repeat for DST frame and observer
        if observe_dst is None:
            # Observer defined by frame
            frame_dst, observe_dst = self.get_system_frame_names(system_dst, precess=precess)
        else:
            # Observer must be specified for this frame
            observe_dst = self.get_naif_body_code(observe_dst)
            # Get naif frame and observ
            frame_dst, observe_dst = self.get_system_frame_names(system_dst, observatory=observe_dst, precess=precess)
                
        # Convert to ephemeris time
        if dates.isscalar:
            et = spice.str2et(dates.isot)
        else:
            et = spice.str2et(dates.isot.tolist())
            
        # If observer changes, first do origin shift.
        if observe_src != observe_dst:
            
            # Get location of observe_dst relative to observe_src in frame_src
            corr = 'NONE'
            if no_velocity:
                observe_src_state, ltime = spice.spkpos(observe_dst, et, frame_src, corr, observe_src)
            else:
                # Velocity requested. This needs spkezr, which only accepts floats. So loop through et and call spkezr for
                #  each et. Preallocate state and ltime, in this case state is a len(et)x6 array.
                if dates.isscalar:
                    observe_src_state, ltime = spice.spkezr(observe_dst, et, frame_src, corr, observe_src)
                else:
                    observe_src_state = np.zeros(coord_src.shape, dtype=float)
                    ltime = np.zeros(dates.size, dtype=float)
                    for i in range(dates.size):
                        observe_src_state[i, :], ltime[i] = spice.spkezr(observe_dst, et[i], frame_src, corr, observe_src)
            
            # Now shift the origin         
            coord_src = coord_src - observe_src_state
        
        # Must loop through dates for the matrix multiplication. Preallocate space for output.
        # Now rotate from src frame to dst frame.
        if dates.isscalar:
            if no_velocity:
                # Get rotation matrix for position only
                transform = spice.pxform(frame_src, frame_dst, et)
                coord_dst = np.matmul(transform, coord_src)
            else:
                # Get rotation matrix for full state
                transform = spice.sxform(frame_src, frame_dst, et)
                coord_dst = np.matmul(transform, coord_src)
        else:
            # Must loop through dates for the matrix multiplication. Preallocate space for output.
            coord_dst = np.zeros(coord_src.shape, dtype=float)
            
            for i in range(dates.size):
                if no_velocity:
                    transform = spice.pxform(frame_src, frame_dst, et[i])
                else:
                    transform = spice.sxform(frame_src, frame_dst, et[i])
                
                coord_dst[i, :] = np.matmul(transform, coord_src[i,:])
        
        return coord_dst
        
    
    def convert_lonlat(self, dates, coord_src, system_src, system_dst, observe_src=None, observe_dst=None, degrees=True):
        """
        Function to convert latitudinal coordinates betwen different reference 
        frames.
        :param dates: Astropy time object of dates(s).
        :param coord_src: Array of latitudinal coordinates (rad, lon, lat) to convert. Should 
                          be either a numpy array of len(dates)*3, or a list of floats.
        :param system_src: String name of coordinate system of coord array.
        :param system_dst: String name of coordinate system to transform coord array to.
        :param observe_src: String name of observatory for origin of system from. Only needed for some systems.
        :param observe_dst: String name of observatory for origin of system to. Only needed for some systems.
        :param degrees: Boolean. If true indicates that units of coord_src, and returned coord_dst, are in degrees.
        :return coord_dst: Array, giving state at each dates. Shape of len(dates)x3.
        """
        
        # If coordinates input as a list, then bung them into an array.
        if isinstance(coord_src, list):
            if all([isinstance(c, (float, int)) for c in coord_src]):
                coord_src = np.array(coord_src)
                coord_src = np.squeeze(coord_src)
            else:
                print("ERROR: coord_src should be a numpy array of coordinates or a list of floats.")
            
        # If coord only has one dimension, set it so that time is zeroth.
        if coord_src.ndim == 1:
            n_coords = 1
            n_components = coord_src.size
        else:
            n_coords = coord_src.shape[0]
            n_components = coord_src.shape[1]
            
        # Check dates and coord sizes match.
        if dates.size != n_coords:
            print("Error: Number of dates does not correspond to number of coordinates.")
            
        # Check coords have 3 components.
        if n_components != 3:
            print("ERROR: Invalid dimension of position vector or state vector")
                
        if (system_src in ['HPC', 'hpc', 'HPR', 'hpr', 'RTN', 'rtn']) & (observe_src is None):
            print("ERROR: system_src given as {}, but no observe_src specified. Assuming Earth".format(system_src))
            observer_src = 'earth'
        
        if (system_dst in ['HPC', 'hpc', 'HPR', 'hpr', 'RTN', 'rtn']) & (observe_dst is None):
            print("ERROR: system_dst given as {}, but no observe_src specified. Assuming Earth".format(system_dst))
            observer_dst = 'earth'
        
        # Parse out the coordinates.
        if dates.isscalar:
            rad = coord_src[0]
            lon = coord_src[1]
            lat = coord_src[2]
        else:
            rad = np.squeeze(coord_src[:,0])
            lon = np.squeeze(coord_src[:,1])
            lat = np.squeeze(coord_src[:,2])
        
        # Put angles into radians for spice.
        if degrees:
            lon = np.deg2rad(lon)
            lat = np.deg2rad(lat)
            
        if system_src in ['HPC', 'hpc']:
            # Convert to RTN and updates system_src tag
            lon, lat = self.convert_hpc_to_rtn(lon, lat, degrees=False)
            system_src = 'RTN'
        elif system_src in ['HPR', 'hpr']:
            # Convert to HPC and then to RTN, updates system_src tag
            lon, lat = self.convert_hpr_to_hpc(lon, lat, degrees=False)
            lon, lat = self.convert_hpc_to_rtn(lon, lat, degrees=False)
            system_src = 'RTN'
        
        # Now convert to rectangular coords.
        if dates.isscalar:
            coord_src_rec = spice.latrec(rad, lon, lat)
        else:
            # Loop through state to do conversion (as spice.reclat doesn't handle arrays yet)
            coord_src_rec = np.zeros(coord_src.shape, dtype=float)
            for i in range(dates.size):
                coord_src_rec[i,:] = spice.latrec(rad[i], lon[i], lat[i])
    
        # If system_dst was HPC or HPR, do conversion in RTN and apply correction.
        if system_dst in ['HPC', 'hpc']:
            system_dst = 'RTN'
            calc_hpc = True
        else:
            calc_hpc = False
            
        if system_dst in ['HPR', 'hpr']:
            system_dst = 'RTN'
            calc_hpr = True
        else:
            calc_hpr = False
        
        coord_dst_rec = self.convert_coord(dates, coord_src_rec, system_src, \
                                           system_dst, observe_src=observe_src, \
                                           observe_dst=observe_dst)
                                           
        # Convert back to latitude coords.
        if dates.isscalar:
            rad_dst, lon_dst, lat_dst = spice.reclat(coord_dst_rec)
        else:
            # Loop through state to do conversion (as spice.reclat doesn't handle arrays yet)
            rad_dst = np.zeros(dates.size, dtype=float)
            lon_dst = np.zeros(dates.size, dtype=float)
            lat_dst = np.zeros(dates.size, dtype=float)
            for i in range(dates.size):
                rad_dst[i], lon_dst[i], lat_dst[i] = spice.reclat(coord_dst_rec[i,:])

        # Correct HPC coords if necessary
        if calc_hpc:
            lon_dst, lat_dst = self.convert_rtn_to_hpc(lon_dst, lat_dst, degrees=False)
            
        # Correct HPR coords if necessary
        if calc_hpr:
            lon_dst, lat_dst = self.convert_rtn_to_hpc(lon_dst, lat_dst, degrees=False)
            lon_dst, lat_dst = self.convert_hpc_to_hpr(lon_dst, lat_dst, degrees=False)

        # Correct Carrington longitudes if neccesary
        carrington_names = ['CARR', 'CARRINGTON', 'carr', 'carrington']
        if system_dst in carrington_names:
            if dates.size > 1:
                id_under = lon_dst < 0
                if any(id_under):
                    lon_dst[id_under] += 2*np.pi

            elif dates.size == 1:
                if lon_dst < 0:
                    lon_dst += 2*np.pi

        if degrees:
            lon_dst = np.rad2deg(lon_dst)
            lat_dst = np.rad2deg(lat_dst)

        # Bundle the output into one array.
        if dates.isscalar:
            coord_dst = np.array([rad_dst, lon_dst, lat_dst])
        else:
            rad_dst = np.expand_dims(rad_dst, axis=1)
            lon_dst = np.expand_dims(lon_dst, axis=1)
            lat_dst = np.expand_dims(lat_dst, axis=1)
            coord_dst = np.hstack((rad_dst, lon_dst, lat_dst))

        return coord_dst
        
    
    def convert_lat2rec(self, coord_lat, degrees=True):
        """
        Function to convert latitudinal coordinates to rectangular coordinates.
        :param coord_lat: Array of latitudinal coordinates (rad, lon, lat) to convert. Should 
                          be either a numpy array of len(dates)*3, or a list of floats.
        :param degrees: Boolean, True if units is in degrees.
        :return coord_rec: Array, giving position (km) for each latitudinal coordinate. Shape of coord_lat.
        """
        
        # If coordinates input as a list, then bung them into an array.
        if isinstance(coord_lat, list):
            if all([isinstance(c, (float, int)) for c in coord_lat]):
                coord_lat = np.array(coord_lat)
                coord_lat = np.squeeze(coord_lat)
            else:
                print("ERROR: coord_src should be a numpy array of coordinates or a list of floats.")
            
        # If coord only has one dimension, set it so that time is zeroth.
        if coord_lat.ndim == 1:
            n_coords = 1
            n_components = coord_lat.size
        else:
            n_coords = coord_lat.shape[0]
            n_components = coord_lat.shape[1]
                 
        # Check coords have 3 components.
        if n_components != 3:
            print("ERROR: Invalid dimension of position vector or state vector")
                
        # Get the coordinates.
        if n_coords == 1:
            rad = coord_lat[0]
            lon = coord_lat[1]
            lat = coord_lat[2]
        else:
            rad = np.squeeze(coord_lat[:,0])
            lon = np.squeeze(coord_lat[:,1])
            lat = np.squeeze(coord_lat[:,2])
        
        # Put angles into radians for spice.
        if degrees:
            lon = np.deg2rad(lon)
            lat = np.deg2rad(lat)
                   
        # Now convert to rectangular coords.
        if n_coords == 1:
            coord_rec = spice.latrec(rad, lon, lat)
        else:
            # Loop through state to do conversion (as spice.reclat doesn't handle arrays yet)
            coord_rec = np.zeros(coord_lat.shape, dtype=float)
            for i in range(n_coords):
                coord_rec[i,:] = spice.latrec(rad[i], lon[i], lat[i])
               
        return coord_rec
            
                
    def convert_rec2lat(self, coord_rec, degrees=True):
        """
        Function to convert rectangular coordinates to latitudinal coordinates.
        :param coord_rec: Array of rectangular coordinates (in km) to convert. Should 
                          be either a numpy array of len(dates)*3, or a list of floats.
        :param degrees: Boolean, True if units of latitudinal coords should be in degrees.
        :return coord_lat: Array, giving latitudinal position for each rectangular coord. Shape of coor_rec.
        """
        
        # If coordinates input as a list, then bung them into an array.
        if isinstance(coord_rec, list):
            if all([isinstance(c, (float, int)) for c in coord_rec]):
                coord_rec = np.array(coord_rec)
                coord_rec = np.squeeze(coord_rec)
            else:
                print("ERROR: coord_src should be a numpy array of coordinates or a list of floats.")
            
        # If coord only has one dimension, set it so that time is zeroth.
        if coord_rec.ndim == 1:
            n_coords = 1
            n_components = coord_rec.size
        else:
            n_coords = coord_rec.shape[0]
            n_components = coord_rec.shape[1]
                 
        # Check coords have 3 components.
        if n_components != 3:
            print("ERROR: Invalid dimension of position vector or state vector")
                                   
        # Now convert to latitudinal coords
        if n_coords == 1:
            coord_lat = spice.reclat(coord_rec)
        else:
            # Loop through state to do conversion (as spice.reclat doesn't handle arrays yet)
            coord_lat = np.zeros(coord_rec.shape, dtype=float)
            for i in range(n_coords):
                coord_lat[i,:] = spice.reclat(coord_rec[i,:])
        
        if degrees:
            if n_coords == 1:
                coord_lat[1:] = np.rad2deg(coord_lat[1:])
            else:
                coord_lat[:,1:] = np.rad2deg(coord_lat[:,1:])
        
        return coord_lat
        
    
    def get_naif_body_code(self, body):
        """
        Function to return the numeric NAIF body code of a spacecraft or solar system body.
        :param body: String name of solar system body or spacecraft
        :return naif_body: String of numeric NAIF body code, fit for use with spiceypy
        """

        # Get set of pseudonyms for stereo a and stereo b and Earth
        sta_names = {'A', 'sta', 'STA', 'Ahead', 'STEREO Ahead', 'STEREO-Ahead', 'STEREO_Ahead', '-234', -234}
        stb_names = {'B', 'stb', 'STB', 'Behind', 'STEREO Behind', 'STEREO-Behind', 'STEREO_Behind', '-235', -235}
        sun_names = {'S', 'Sun', 'sun', 'SUN', '10', 10}
        mercury_names = {'Mercury', 'mercury', 'MERCURY', '199', 199}
        venus_names = {'Venus', 'venus', 'VENUS', '299', 299}
        earth_names = {'E','Earth', 'earth', 'EARTH', 'ERT', 'ert', '399', 399}
        moon_names = {'Moon', 'moon', 'Moon', '301', 301}
        mars_names = {'Mars', 'mars', 'MARS', '499', 499}
        # Get NAIF code of user input observatory
        if body in sta_names:
            naif_body = '-234'
        elif body in stb_names:
            naif_body = '-235'
        elif body in sun_names:
            naif_body = '10'
        elif body in mercury_names:
            naif_body = '199'
        elif body in venus_names:
            naif_body = '299'
        elif body in earth_names:
            naif_body = '399'
        elif body in moon_names:
            naif_body = '301'
        elif body in mars_names:
            naif_body = '499'
        else:
            print(body)
            print('ERROR: Body name not recognised. Allowed bodies: STA, STB, Sun, Mercury, Venus, Earth, Moon, Mars.\n\
                  More codes available at Get more codes at http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html ')
            
            naif_body = 'ERROR'

        return naif_body
        

    def get_coord(self, dates, target, system, observatory=None, corr='NONE', precess=False, return_ltime=False,
                  no_velocity=False):
        """
        Function to calculate the coordinates of a spacecraft or solar system body
        :param dates: Astropy time object of dates(s).
        :param target: String name of target body to calculate the position of.
        :param system: String name of coordinate system to calculate position in.
        :param observatory: String name of observatory for origin of coordinate system.
        :param corr: Boolean. If True performs correction for planetary abberation.
        :param precess: Boolean. If True accounts for precession in ..
        :param return_ltime: Boolean. If True also returns light travel time between observatory and target.
        :param no_velocity: Boolean. If True, the velocity components are not returned in the state vector.
        :return state: Array, giving state at each dates. Either len(dates)x3 or len(dates)x6, depending on no_velocity
        :return ltime (optional): The light travel time between observatory and target
        """
        # TODO: Add in error checks to make sure target is specified.

        # The NAIF codes for stereo a and b are:
        sta_naif_code = '-234'
        stb_naif_code = '-235'

        target = self.get_naif_body_code(target)

        if observatory is not None:
            observatory = self.get_naif_body_code(observatory)

        frame, observatory = self.get_system_frame_names(system, observatory, precess)

        # Pull out the max valid ephemeris time for each craft
        if target == sta_naif_code:
            max_dates = self.__PredMaxdates__['sta']
        else:
            max_dates = self.__PredMaxdates__['stb']

        # Get ephemeris times for the input dates. If larger than max dates, set to max dates and set correct flag.
        correct_late_times = False
        if dates.isscalar:
            if dates <= max_dates:
                et = [spice.str2et(dates.isot)]
            else:
                et = [spice.str2et(max_dates.isot)]
                correct_late_times = True
        else:
            if all(dates <= max_dates):
                et = spice.str2et(dates.isot.tolist())
            else:
                et = spice.str2et(dates.isot.tolist())
                id_gt_max = dates > max_dates
                et[id_gt_max] = spice.str2et(max_dates.isot)
                correct_late_times = True

        # Now add in the calls to spice functions to get the state variables
        # If no velocity, use spkpos.
        if no_velocity:
            state, ltime = spice.spkpos(target, et, frame, corr, observatory)
            # Convert state from list of arrays to numpy array
            state = np.array(state)
        else:
            # Velocity requested. This needs spkezr, which only accepts floats. So loop through et and call spkezr for            
            state = np.zeros((dates.size, 6), dtype=float)
            ltime = np.zeros(dates.size, dtype=float)
            for i in range(dates.size):
                state[i, :], ltime[i] = spice.spkezr(target, et[i], frame, corr, observatory)

        # Add in correction to get states for times after max ephemeris time using
        # Conics. Loop over entries after maximum time
        if correct_late_times:
            if observatory == sta_naif_code:
                elts = self.__PredConic__['sta']
            elif observatory == stb_naif_code:
                elts = self.__PredConic__['stb']

        
            for idt in np.argwhere(id_gt_mxt is True):
                et = spice.str2et(dates.isot[idt])
                temp = spice.conics(elts, et)
                # Updates state with conics estimate.
                if no_velocity:
                    state[idt, :] = temp[0:3]
                else:
                    state[idt, :] = temp

                ltime[idt] = -1

            # Now convert updatesd coords from HAE to specified system
            temp = state[:, id_gt_mxt]
            temp = convert_stereo_coord(dates.isot[id_gt_mxt], temp, 'HAE', system)
            state[:, id_gt_mxt] = temp
        
        # If only one date, loose first dimension
        if dates.isscalar:
            state = np.squeeze(state)
            
        if return_ltime:
            return state, ltime
        else:
            return state
            

    def get_lonlat(self, dates, target, system, observatory=None, corr='NONE', precess=False, degrees=True):
        """
        Function to calculate the radius, longitude and lataitude of a target in coordinate system given by system,
        centered on an observatory. Observatory doesn't always need to be specified, but for some coordinate systems
        needs to be. Doesnt handle velocity components of state vector.
        :param dates: Astropy time object of dates(s) to get
        :param target: String name of target body.
        :param system: String name of coordinate system.
        :param observatory: String name of observatory.
        :param corr: String specifying whether to perform correction for planetary abberation.
        :param precess: Boolean. If true perform a calculation for precession.
        :param degrees: Boolean. If true return latitude and longitude in degrees
        :return coords: Float array of coords.
        """

        # Add in calculation for Helioprojective Cartesian coordinates.
        if system in ['HPC', 'hpc']:
            system = 'RTN'
            calc_hpc = True
        else:
            calc_hpc = False
            
        if system in ['HPR', 'hpr']:
            system = 'RTN'
            calc_hpr = True
        else:
            calc_hpr = False
            
        # Get the position of the target for this dates/system/observatory
        state = self.get_coord(dates, target, system=system, observatory=observatory, corr=corr, precess=precess,
                               no_velocity=True)
		
		
        # Use spice to convert to lon/lat
        state = np.array(state)
        
        if state.ndim == 1:
            rad, lon, lat = spice.reclat(state)
        else:
            # Loop through state to do conversion (as spice.reclat doesn't handle arrays yet)
            rad = np.zeros(len(dates), dtype=float)
            lon = np.zeros(len(dates), dtype=float)
            lat = np.zeros(len(dates), dtype=float)
            for i in range(len(dates)):
                rad[i], lon[i], lat[i] = spice.reclat(state[i])

        # Correct HPC coords if necessary
        if calc_hpc:
            lon, lat = self.convert_rtn_to_hpc(lon, lat, degrees=False)
            
        if calc_hpr:
            lon, lat = self.convert_rtn_to_hpc(lon, lat, degrees=False)
            lon, lat = self.convert_hpc_to_hpr(lon, lat, degrees=False)

        # Correct Carrington longitudes if neccesary
        carrington_names = ['CARR', 'CARRINGTON', 'carr', 'carrington']
        if system in carrington_names:
            if dates.size > 1:
                id_under = lon < 0
                if any(id_under):
                    lon[id_under] += 2*np.pi

            elif dates.size == 1:
                if lon < 0:
                    lon += 2*np.pi

        if degrees:
            lon = np.rad2deg(lon)
            lat = np.rad2deg(lat)

        # Make state vector with same structure as output by spkpos
        if dates.isscalar:
            coords = np.array([rad,lon,lat]).T
            # If only one date, loose first dimension
            coords = np.squeeze(coords)
        else:
            rad = np.expand_dims(rad, axis=1)
            lon = np.expand_dims(lon, axis=1)
            lat = np.expand_dims(lat, axis=1)
            coords = np.hstack((rad,lon,lat))
        return coords
        

    def get_system_frame_names(self, system, observatory=None, precess=False):
        """
        Function to return the spice friendly name of the frame relevant for a given coordinate system and observatory.
        For some systems, an observatory must also be specified.
        :param system: String name of the requested coordinate system.
        :param observatory: String name of observatory. Not needed for all systems (e.g HCI or GEO).
        :param precess: Boolean. If true select frame that accounts for precession.
        :return frame: String name of frame in spice friendly format.
        :return observatory: String of observatory in spice friendly format.
        """
        
        # Get NAIF codes needed to identify which observatory was parsed.
        sta_naif_code = self.get_naif_body_code('STA')
        stb_naif_code = self.get_naif_body_code('STB')
        sun_naif_code = self.get_naif_body_code('SUN')
        earth_naif_code = self.get_naif_body_code('EARTH')
                
        # If an observatory was parsed, make sure it is in naif format.
        if observatory is not None:
            observatory = self.get_naif_body_code(observatory)
         
        # Get lists of system names and their synonyms
        heeq_names = ['HEQ','HEEQ', 'heq', 'heeq']
        carrington_names = ['CARR', 'CARRINGTON', 'carr', 'carrington']
        hci_names = ['HCI', 'hci']
        hae_names = ['HAE', 'hae']
        hee_names = ['HEE', 'hee']
        hgrtn_names = ['HGRTN', 'hgrtn']
        rtn_names = ['RTN', 'rtn']
        sci_names = ['SCI', 'sci']
        hertn_names = ['HERTN', 'hertn']
        gei_names = ['GEI', 'gei']
        geo_names = ['GEO', 'geo']
        gse_names = ['GSE', 'gse']
        
        # If system is RTN then an observatory should have been specified.
        if (system == 'RTN') & (observatory is None):
            print("ERROR: RTN system specified, but no corresponding observatory. Assuming Earth.")
            observatory = earth_naif_code
            
        # Some systems have defined observatory. Check input observatory doesnt clash with this.
        # Helio systems:
        all_helio = heeq_names + carrington_names + hci_names + hae_names + hee_names + hgrtn_names
        if (system in all_helio) & (observatory not in [None, sun_naif_code]):
            print("ERROR: Observatory {0} specified for Sun based system ({1}).\n\
            Input observatory will be overridden".format(observatory, system))
            
        all_geo = gei_names + geo_names + gse_names
        if (system in all_geo) & (observatory not in [None, earth_naif_code]): 
            print("ERROR: Observatory {0} specified for Earth based system ({1}).\n\
            Input observatory will be overridden".format(observatory, system))

        # Find out which coordinate system / frame required
        if system in heeq_names:
            frame = 'HEEQ'
            observatory = 'Sun'

        elif system in carrington_names:
            frame = 'IAU_SUN'
            observatory = 'Sun'

        elif system in hci_names:
            frame = 'HCI'
            observatory = 'Sun'

        elif system in hae_names:
            observatory = 'Sun'
            if precess:
                frame = 'ECLIPdates'
            else:
                frame = 'ECLIPJ2000'

        elif system in hee_names:
            frame = 'HEE'
            observatory = 'Sun'
        
        elif system in rtn_names:
            if observatory == sta_naif_code:
                frame = 'STAHGRTN'
            elif observatory == stb_naif_code:
                frame = 'STBHGRTN'
            elif observatory == earth_naif_code:
                frame = 'GEORTN'
            elif observatory == sun_naif_code:
                frame = 'HGRTN'
        
        elif system in hgrtn_names:
                frame = 'HGRTN'
                observatory = 'Sun'

        elif system in sci_names:
            # Check observaotry input
            if observatory == sta_naif_code:
                frame = 'STASCPNT'
            elif observatory == stb_naif_code:
                frame = 'STBSCPNT'
            else:
                print('ERROR: INCORRECT SPACECARFT TARGET PARSED')

        elif system in hertn_names:
            # Check observatory input
            if observatory == sta_naif_code:
                frame = 'STAHERTN'
            elif observatory == stb_naif_code:
                frame = 'STBHERTN'
            else:
                print('ERROR: INCORRECT SPACECARFT TARGET PARSED')

        elif system in gei_names:
            frame = 'J200'
            observatory = 'Earth'

        elif system in geo_names:
            frame = 'J200'
            observatory = 'Earth'

        elif system in gse_names:
            frame = 'GSE'
            observatory = 'Earth'
        else:
            print("ERROR: System not recognised")

        return frame, observatory
        

    def clear_kernals(self):
        """
        A function to clear out the loaded spice kernals.
        :return:
        """
        spice.kclear()
        return
    