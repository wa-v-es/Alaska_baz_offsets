import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin, Catalog
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
import numpy as np
from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.earthworm import Client as EW_Client
from obspy.clients.fdsn.header import FDSNNoDataException
import datetime
import os
import sys
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
from obspy.io.sac.util import obspy_to_sac_header
import warnings

# Define some conversion factors
KM2M = 1000     # [m/km]
SEC2MIN = 1/60  # [min/s]

def get_aperture(stlat,stlong):
    stlat= [eval(i) for i in stlat]
    stlong= [eval(i) for i in stlong]
    max=[]
    for i in range(len(stlat)):
        lat,long=stlat[i],stlong[i]
        dist=[]
        for j in range(len(stlat)):

            d,v,b=gps2dist_azimuth(lat, long, stlat[j], stlong[j])
            dist.append(d)

        max.append(np.max(dist))
    return np.max(max),max

def get_one_st(st):
    stream=Stream()
    stream.append(st[0])
    return stream

def gather_stations_bulk(lon_0, lat_0, max_radius, starttime, endtime,
                          channel, network='*', station='*', location='*',
                          time_buffer=0, merge_fill_value=0, trim_fill_value=0,
                          remove_response=False, watc_url=None,
                          watc_username=None, watc_password=None, iris_only=True,
                          verbose=True):
    """
    Args:
        lon_0 (int or float): Longitude of search center [deg.]
        lat_0 (int or float): Latitude of search center [deg.]
        max_radius (int or float): Maximum radius to search for stations within
            [km]
        time_buffer (int or float): Extra amount of data to download after
            `endtime` [s]
        trim_fill_value (bool, int, float, or None): Controls trimming of the
            output :class:`~obspy.core.stream.Stream`, useful if precisely
            uniform start and end times are desired. If `False`, no trimming is
            performed. Otherwise, a trim is performed with the ``fill_value``
            provided to this parameter. For details, see the docstring of
            :meth:`obspy.core.stream.Stream.merge`
        remove_response (bool): Toggle response removal via
            :meth:`~obspy.core.trace.Trace.remove_sensitivity` or a simple
    """
    # log() does nothing if `verbose=False`
    log = print if verbose else lambda *args, **kwargs: None

    # Check for issues with fill value args
    if merge_fill_value is True or trim_fill_value is True:
        raise ValueError('Cannot provide True to fill value parameters.')

    log('-------------------')
    log('BULK GATHERING DATA')
    log('-------------------')

    log('Creating station list...')

    inventories = []  # Create empty list of inventories

    # Grab IRIS inventory
    iris_client = FDSN_Client('IRIS')
    try:
        iris_inv = iris_client.get_stations(starttime=starttime,
                                            endtime=endtime + time_buffer,
                                            network=network, station=station,
                                            location=location, channel=channel,
                                            level='channel')

        inventories.append(iris_inv)  # Add IRIS inventory to list

    except FDSNNoDataException:
        log('No stations found on IRIS FDSN.')

    # If the user supplied both a WATC password and WATC username, then search
    # through WATC database
    if watc_username and watc_password:

        # Grab WATC inventory
        log('Connecting to WATC FDSN...')
        watc_client = FDSN_Client(base_url=watc_url, user=watc_username,
                                  password=watc_password)
        log('Successfully connected.')
        try:
            watc_inv = watc_client.get_stations(starttime=starttime,
                                                endtime=endtime + time_buffer,
                                                network=network, station=station,
                                                location=location, channel=channel,
                                                level='channel')

            inventories.append(watc_inv)  # Add WATC inventory to list

        except FDSNNoDataException:
            log('No stations found on WATC FDSN.')

    requested_station_list = []  # Initialize list of stations to request

    # Big loop through all channels in all inventories!
    for inv in inventories:
        for nw in inv:
            for stn in nw:
                for cha in stn:
                    dist, _, _ = gps2dist_azimuth(lat_0, lon_0, cha.latitude,
                                                  cha.longitude)  # [m]
                    if dist <= max_radius * KM2M:
                        requested_station_list.append(stn.code)

    if not requested_station_list:
        raise ValueError('Station list is empty. Expand the station search '
                         'and try again.')

    # Put into the correct format for ObsPy (e.g., 'HOM,O22K,DLL')
    requested_stations = ','.join(np.unique(requested_station_list))

    return requested_stations, inventories
    log('Done')

###
def utc2str(t):
    str=str(t.year)+t.month+t.day+t.hour+t.minute



###
#######################################
#rm events_6.4mw_ilar.xml
catfile = 'events_6.5mw_ilar.xml'
catfile = 'events_6.7_2024mw_AK_allmag.xml'
# catfile = 'B_R_fig3.xml'

ar_lat= 64
ar_long=-152

client = Client("IRIS")

#
# station_name=[]
# for i in range(len(inventory.get_contents()['stations'])):
#     station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])
#
# station_name.remove('ILAR')#remove beam
# station_name.remove('IL31')#remove beam

#earthquakes
starttime= UTCDateTime('2020-01-01T00:00:01')
endtime= UTCDateTime('2024-07-10T00:00:01')

if not os.path.exists(catfile):
    catalog = client.get_events(starttime=starttime, endtime=endtime,minmagnitude=6.7,maxmagnitude=9.9,\
    latitude=ar_lat,longitude=ar_long, minradius=65, maxradius=120,mindepth=60)
    catalog.write(catfile, 'QUAKEML')

catalog=read_events(catfile)
print(catalog.__str__(print_all=True))

# sys.exit()

# stlist,inv = gather_stations_bulk(lon_0=-146.886597, lat_0=64.7714, max_radius=150,starttime=starttime, endtime=starttime+100,network='TA,AK',channel='H*Z,B*Z',iris_only=True)

# st_lat=[]; st_long=[]; st_elev=[]; station_name=[]; st_netw=[]
#
# for line in open('alaska_250.txt','r'):
#     line=line.split()
#     if line[0]=='#':
#         continue
#     station_name.append(line[0])
#     st_lat.append(line[1])
#     st_long.append(line[2])
#     st_elev.append(line[3])
#     st_netw.append(line[4])
#
# max,max_list=get_aperture(st_lat,st_long)
# file.close()

# print('Aperture of array = {} km'.format(int(max/1000)))


inventory_big = client.get_stations(network="AK,CN,AV,AT",minlatitude=55,maxlongitude=-137,
minlongitude=-168,starttime=starttime,endtime=endtime)#,level='response')
st_elev=[]
for nw in inventory_big:
    for stn in nw:
        st_elev.append(stn.elevation)
st_elev=np.array(st_elev)
print(f'Min elevation {np.min(st_elev)}.. Max: {np.max(st_elev)}')
sys.exit()
c_one=Catalog()
# c_one.append(catalog[4])
# c_one.append(catalog[3])
c_one.append(catalog[23])
# c_one.append(catalog[18])


#remove 4 events I've done before# 4,
# c_two.append(catalog[7])
# sys.exit()

# 2023-07-02T10:27:43.732000Z | -17.881, -174.948 | 6.9  mww
# st_list=['I23K','PS07','POKR','MDM','NEA2','MCK','WRH','CCB','HDA','PS08']
for event in c_one:
    # if event in c_one:
    #     print('done',event.preferred_origin().time)
    #     print('####')
    #     continue
    print('doing rest')
    ori= event.preferred_origin()
    mag=event.preferred_magnitude().mag
    if mag < 7.25:
        continue
    print('Earthquake depth:',ori['depth']/1000,'km\n')
    if -120<ori.longitude <-30:
        print('earthquake in SA\n')

    ##########
    ##
    print('#########  \n')
    print('Earthquake deets: \n')
    print(ori)

    ###
    st_all=Stream()
    # station_name=['IL01'] #[16:]
    # for st in inventory_big.select(network='CN').get_contents()['stations']:
    for st in inventory_big.get_contents()['stations']:
    # for st in st_list:


        st_name=st.split()[0]
        netw=st_name[:2]
        st_name=st_name[3:]

        # st_name=st
        netw='AK'
        # netw=st_netw[i]

        print('#########  \n')
        print('Doing', st_name)
        print('#########  \n')

        try:
            inventory = client.get_stations(network=netw,station=st_name,level='response',channel="B*Z,H*Z")#, minlatitude=55)#,level='stations')
        except:
            print('ehh: B or H channel not found for station: {}'.format(st_name))
            continue
        # netw=inventory.get_contents()['networks'][0]
        channel=inventory.get_contents()['channels'][0]
        # try:
        station_coord=inventory.get_coordinates(channel)
        # except:
        # station_coord=inventory.get_coordinates('{}.{}..H*Z'.format(netw,st_name))


        model = TauPyModel(model="ak135")
        dist=calc_dist(ori.latitude,ori.longitude,station_coord['latitude'],station_coord['longitude'],6400,0)
        dist_azi_bazi=calc_dist_azi(ori.latitude,ori.longitude,station_coord['latitude'],station_coord['longitude'],6400,0)

        arrivals = model.get_travel_times(source_depth_in_km=ori.depth/1000,distance_in_degree=dist,phase_list=["PP"])
        arr_PP=arrivals[0]
        try:
            arr_pP=model.get_travel_times(source_depth_in_km=ori.depth/1000,distance_in_degree=dist,phase_list=["pP"])[0]
        except:
            print('no pP :/')
        arris_P = model.get_travel_times(source_depth_in_km=ori.depth/1000,distance_in_degree=dist,phase_list=["P"])
        try:
            arr_P=arris_P[0]
        except:
            print('no P')
        # arr_PP=arrivals[1]
        pad=30 # pad in sec
        ##get waveforms P -pad to PP + pad
        ### get waveform PP - 300 to PP + 30. P doesn't work for dist >95 deep quakes
        starttime=ori.time+arr_PP.time-400
        # starttime=ori.time+arr_P.time+30
        endtime=ori.time+arr_PP.time+200
        try:
            st=client.get_waveforms(netw, st_name, "*", "**Z", starttime,endtime,attach_response=True)
            print('got st')
            if len(st) > 1:
                print('Found more than one traces for station {}\n'.format(st_name))
                print('saving first')
                st=get_one_st(st)

            st[0].stats.sac = obspy_to_sac_header(st[0].stats)
            # prepare header for sac
            #tr.stats.sac.o = (origin.time - tr.stats.starttime)
            st.filter('bandpass', freqmin=0.005, freqmax=9.5)
            sachdr = {
                'b':arr_PP.time-400,
                'evla': ori.latitude,
                'evlo': ori.longitude,
                'evdp': ori.depth,
                'stla': station_coord['latitude'],
                'stlo': station_coord['longitude'],
                'stel': inventory.select(station=st_name)[0][0].elevation,
                'stdp': 0.0,
                'mag': mag,
               'gcarc':dist_azi_bazi[0],
               'dist':dist_azi_bazi[0]*111,
               'az':dist_azi_bazi[1],
               'baz':dist_azi_bazi[2]
                }
            st[0].stats.sac=sachdr
            st[0].stats.sac = obspy_to_sac_header(st[0].stats)
            pre_filt = [0.01, 0.02, 9, 10]
            inv=inventory.select(station=st_name)
            st.detrend('linear')
            st.detrend('demean').taper(type='cosine',max_percentage=0.05)
            # st.taper(type='cosine',max_percentage=.1)#increase taper..
            time_string=UTCDateTime.strftime(ori.time,'%y%m%d_%H%M%S')

            try:
                # st[0].remove_response(inventory=inv,pre_filt=pre_filt,water_level=None, plot='response_{}_{}.png'.format(time_string,st_name))
                st[0].remove_response(inventory=inv,pre_filt=pre_filt,water_level=None)
            except:
                print('response removal didnt work')

            if np.max(st[0].data) < 1e-4 and np.max(st[0].data) !=0:
                st_all.extend(st)
            ###

        except:
            print('data not found for station {}'.format(st_name))


    st_all.resample(20.0)
    print(len(st_all), 'traces for earthquake:', ori.time,'\n')
    print('########\n')
    ####

    # plt.show()
    # sys.exit() # exit here to see

    time_string=UTCDateTime.strftime(ori.time,'%y%m%d_%H%M%S')
    for tr in st_all:

        if -120 <ori.longitude <-30: # if earthquake in South America, add SA in folder name
            os.makedirs(os.path.join('sac_files/{}_SA/'.format(time_string), ''),exist_ok=True)
            data = os.path.join('sac_files/{}_SA/'.format(time_string), '')
            sst='{}_{}.sac'.format(time_string,tr.id)
            tr.write(data+sst,format='sac')
        else:
            os.makedirs(os.path.join('sac_files/{}_PA/'.format(time_string), ''),exist_ok=True)
            data = os.path.join('sac_files/{}_PA/'.format(time_string), '')
            sst='{}_{}.sac'.format(time_string,tr.id)
            tr.write(data+sst,format='sac')



    folder_path='sac_files/eq_{}_seismo'.format(time_string[:6])
    os.makedirs(folder_path, exist_ok=True)
    N_old=0
    for i in range(int((len(st_all) + 14) / 15)):  # Adjusted range to include all chunks
        N_new = (i + 1) * 15

        if N_new >= len(st_all):  # last bit
            st_temp = st_all[N_old:]
        else:
            st_temp = st_all[N_old:N_new]  # current chunk

        N_old = N_new

        fig = plt.figure(figsize=(11,8))
        plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        st_temp.plot(color='cadetblue',bgcolor='ivory', fig=fig, orientation="horizontal")
        # try:
        #     plt.axvline(ori.time+arr_P.time,.1,.9,ls='--',c='darkred')
        # except:
        #     print(' ')
        # plt.axvline(ori.time+arr_pP.time,.1,.9,ls='--',c='darkred')
        # plt.axvline(ori.time+arr_PP.time,.1,.9,ls='--',c='darkred')
        fig_name=os.path.join(folder_path,'eq_{}_{}.png'.format(time_string,i))
        plt.savefig(fig_name, dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        # plt.show()
        # arr_P=[]
        # print('Eq done\n')
    ###
    # sys.exit()

# cp *sac ../VespapackBasic/sac_try/
# sys.exit()
#following is a script to remove some traces within a folder using station name as key
# had to do this coz when I saved traces, I didn't implement tr.data.max < 1e-4.

# st_ref=Stream()
#
# for tr in st_all:
#     if np.max(tr.data) < 1e-4:
#         st_ref.append(tr)
#     else:
#         file_path = 'sac_files/230702_102743_PA/*{}*'.format(tr.stats.station)
#         for filename in Path(".").glob(file_path):
#             os.remove(filename)
#             print('removed station: {}'.format(tr.stats.station))
#
# print(len(st_ref))

## script to save list of stations as text file
#
# with open('AK_4arrays.txt', 'w') as file:
#     for st in inventory_big.get_contents()['stations']:
#         st_name=st.split()[0]
#         netw=st_name[:2]
#         st_name=st_name[3:]
#         try:
#             inventory = client.get_stations(network=netw,station=st_name,level='response',channel="B*Z,H*Z")#, minlatitude=55)#,level='stations')
#         except:
#             print('ehh: B or H channel not found for station: {}'.format(st_name))
#             continue
#         channel=inventory.get_contents()['channels'][0]
#         station_coord=inventory.get_coordinates(channel)
#         file.write(f"{station_coord['latitude']:.3f} {station_coord['longitude']:.3f} {inventory.select(station=st_name)[0][0].elevation} {st_name} {netw}\n")
# file.close()
