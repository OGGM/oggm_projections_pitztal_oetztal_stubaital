import oggm.cfg as cfg
from oggm import utils, workflow, tasks, global_tasks, graphics, entity_task
from oggm.core import gis
from oggm.core.flowline import FileModel
import logging
import geopandas as gpd
import os
import sys
import matplotlib.pyplot as plt
from oggm.shop import gcm_climate
from time import gmtime, strftime
import xarray as xr
import numpy as np
import pandas as pd
import salem
import shapely.geometry as shpg

# help function to get area with min h
log = logging.getLogger(__name__)

@entity_task(log)
def process_spartacus_data(gdir, output_filesuffix=None, buffer=500, fp_clim=None, fp_dem=None,
                           aggregation_method='mean', add_plot=False):
    """Processes and writes climate data from the SPARTACUS dataset.
    
    """
    
    dataset = 'SPARTACUS'
    tvar = 'Tm'
    pvar = 'RR'

    # read outline and add buffer
    outline = gdir.read_shapefile('outlines')
    outline_buffered = salem.transform_geopandas(
        gpd.GeoDataFrame(geometry=outline.to_crs("EPSG:32634").buffer(buffer)),
        to_crs=salem.wgs84)

    # open and flatten climate datasets and add dem
    with xr.open_dataset(fp_clim) as ds:
        ds_spartacus = ds.load().stack(points=('x', 'y'))

    with xr.open_dataset(fp_dem) as ds:
        ds_dem = ds.load().stack(points=('x', 'y'))

    ds_spartacus['dem'] = ds_dem.loc[{'points': ds_spartacus.points}].dem    

    # get closest grid point to center point for reference
    lon = gdir.cenlon + 360 if gdir.cenlon < 0 else gdir.cenlon
    lat = gdir.cenlat

    c = (ds_spartacus.lon - lon)**2 + (ds_spartacus.lat - lat)**2
    ds_center_point = ds_spartacus.isel(points=np.argmin(c.data))

    ref_lon = float(ds_center_point.lon.values)
    ref_lat = float(ds_center_point.lat.values)

    # get data points inside from outline
    in_outline = [outline_buffered.geometry.contains(point)[0] for
                  point in gpd.points_from_xy(ds_spartacus.lon, ds_spartacus.lat,
                                              crs='EPSG:3416')]
    if sum(in_outline) == 0:
        # if no point is inside the outline we just use the nearest grid point
        ds_in_outline = ds_center_point
    else:
        # select all points inside of outline
        ds_in_outline = ds_spartacus.isel(points=np.where(in_outline)[0])

    

    # condense all climate data points to single timeseries
    # reference height is the mean height of all
    hgt = ds_in_outline.mean(dim='points', keep_attrs=True).dem.values
    
    if aggregation_method == 'mean':
        # precipitation alculate the mean
        prcp = ds_in_outline.mean(dim='points', keep_attrs=True)[pvar].values

        # for aggregating temperature take laps rate into account
        temp_grad = cfg.PARAMS['temp_default_gradient']
        temp = (ds_in_outline[tvar] + temp_grad * (hgt - ds_in_outline.dem)).mean(dim='points', keep_attrs=True).values
    elif aggregation_method == 'median':
        raise NotImplementedError(f'median aggregation method not implemented so far')
    else:
        raise NotImplemetedError(f'aggregation method not understood: {aggregation_method}')

    time = ds_in_outline.time.data

    # OK, ready to write
    gdir.write_monthly_climate_file(time, prcp, temp, hgt, ref_lon, ref_lat,
                                    filesuffix=output_filesuffix,
                                    source=dataset)
    
    if add_plot:
        plot_margin = 0.02

        all_points_plot = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(ds_spartacus.lon, ds_spartacus.lat,
                                        crs='EPSG:3416'))
        selected_points_plot = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(ds_in_outline.lon, ds_in_outline.lat,
                                        crs='EPSG:3416'))
        outline_plot = salem.transform_geopandas(outline, to_crs=salem.wgs84)
        center_point_plot = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy([ds_center_point.lon], [ds_center_point.lat],
                                        crs='EPSG:3416'))
        
        fig, ax = plt.subplots(figsize=(10,10))
        outline_buffered.plot(ax=ax, color='C2')
        outline_plot.plot(ax=ax, color='C0')
        all_points_plot.plot(ax=ax, color='black', label='climate input grid')
        selected_points_plot.plot(ax=ax, color='C3', label='selected input points')
        center_point_plot.plot(ax=ax, color='C1', label='reference grid point')
        ax.set_xlim([outline_buffered.min_x.values[0] - plot_margin,
                     outline_buffered.max_x.values[0] + plot_margin])
        ax.set_ylim([outline_buffered.min_y.values[0] - plot_margin,
                     outline_buffered.max_y.values[0] + plot_margin])
        ax.legend()
        plt.show()
        
        temp_points = ds_in_outline[tvar] + temp_grad * (hgt - ds_in_outline.dem)
        
        fig, axs = plt.subplots(2, 1, figsize=(10, 10))
        
        for point in temp_points.points:
            axs[0].plot(time, temp_points.loc[{'points': point}])
        axs[0].plot(time, temp, label='mean')
        axs[0].legend()
        for point in ds_in_outline.points:
            axs[1].plot(time, ds_in_outline.loc[{'points':point}][pvar])
        axs[1].plot(time, prcp, label='mean')


@entity_task(log)
def calculate_area_min_h(gdir, filesuffix, min_thickness=1):
    """TODO
    """
    fp_diag = gdir.get_filepath('model_diagnostics',
                                filesuffix=filesuffix)
    with xr.open_dataset(fp_diag) as ds:
        ds_diag = ds.load()
    
    fp = gdir.get_filepath('model_geometry',
                           filesuffix=filesuffix)
    fmod = FileModel(fp)
    
    area_min_h = []
    for year in ds_diag.time.values:
        fmod.run_until(year)
        fl = fmod.fls[0]
        area_min_h.append(
            np.sum(
                np.where(fl.thick > min_thickness,
                         fl.widths_m * fl.dx_meter,
                         0)
            )
        )

    ds_diag[f'area_m2_min_h'] = (('time'), area_min_h)
    #ds_diag[f'area_min_{min_thickness}_m'] = (('time'), area_min_h)
    
    ds_diag.to_netcdf(fp_diag)

# set directory of INPUTDATA
INPUT_DIR = os.environ.get('INPUTDIR')

# to import the user provided function for spartacus data
sys.path.append(INPUT_DIR)
from process_spartacus import process_spartacus_data

# Define data files (should be stored in INPUTDIR)
DEM_file = 'tir1718_clip_extratiles_5m_crs.tif'
Outline_file = 'mergedGI5_withYEAR_3.shp'
Outline_2006_file = 'mergedGI3_withYEAR_3.shp'
Outline_Hochjoch_2017_file = 'Hochjoch_clipGI_5.shp'
Outline_Hochjoch_2006_file = 'Hochjoch_clipGI_3.shp'
dv_Hochjoch_file = 'Hochjoch_clipped.csv'
Volume_file = 'GI_3_all_h.shp'
dv_97_06_file = 'volumechange1997_2006.csv'
dv_06_17_file = 'volumechange2006_2017.csv'
Clim_file = 'SPARTACUS_Oetztal_Stubai_196101_202212.nc'
Clim_dem_file = 'SPARTACUS-MONTHLY_DEM_MASK.nc'

output_file_suffix = '_v1'

cfg.initialize(logging_level='WARNING')

cfg.PARAMS['trapezoid_lambdas'] = 4
cfg.PARAMS['dynamic_spinup_min_ice_thick'] = 2

prcp_fac = 2.4
temp_bias = 0

# Local working directory (where OGGM will work)
if temp_bias < 0:
    temp_bias_suffix = 'm'
else:
    temp_bias_suffix = ''
WORKING_DIR = os.path.join(os.environ["OGGM_WORKDIR"],
                           'working_dir'
)

utils.mkdir(WORKING_DIR)
cfg.PATHS['working_dir'] = WORKING_DIR

# OGGM params
cfg.PARAMS['border'] = 250

# activate multiprocessing
cfg.PARAMS['use_multiprocessing'] = True

# Make it robust
cfg.PARAMS['continue_on_error'] = True

# This is important! We tell OGGM to recompute the glacier area for us
cfg.PARAMS['use_rgi_area'] = False

# We need to keep multipolygons, otherwise the volume estimate is not correct
cfg.PARAMS['keep_multipolygon_outlines'] = True
#cfg.PARAMS['elevation_band_flowline_binsize'] = 10
cfg.PARAMS['map_proj'] = 'utm'

cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['store_fl_diagnostics'] = True
cfg.PARAMS['store_fl_diagnostic_variables'] = [
    'area', 'thickness', 'volume',
    #'volume_bsl', 'volume_bwl', 'calving_bucket',
    #'ice_velocity', 'dhdt', 'climatic_mb', 'flux_divergence'
]

cfg.PARAMS['use_intersects'] = False

# Define DEM source
cfg.PATHS['dem_file'] = os.path.join(INPUT_DIR, DEM_file)

# If we want to use CRU, default is W5E5 with a variable prcp factor
# cfg.PARAMS['baseline_climate'] = 'CRU'  # HISTALP only up to 2014
cfg.PARAMS['baseline_climate'] = 'CUSTOM'
cfg.PARAMS['prcp_fac'] = prcp_fac
cfg.PARAMS['use_winter_prcp_fac'] = False

temp_bias = temp_bias

# this geodetic_mb_period is overwritten for each glacier
cfg.PARAMS['geodetic_mb_period'] = '2006-01-01_2017-01-01'
cfg.PARAMS['hydro_month_nh']=1
cfg.PARAMS['hydro_month_sh']=1


# ---------------  build initial dataframe including all relevant observations  ---------------------
# load outlines and prepare for OGGM
df_outlines = gpd.read_file(os.path.join(INPUT_DIR, Outline_file))

# replace hochjoch with clipped outline
df_outline_hochjoch = gpd.read_file(os.path.join(INPUT_DIR,
                                                 Outline_Hochjoch_2017_file))
hochjoch_index = df_outlines.nr == df_outline_hochjoch.nr.values[0]
df_outlines.loc[hochjoch_index, 'area'] = df_outline_hochjoch.area.values[0]
df_outlines.loc[hochjoch_index, 'geometry'] = df_outline_hochjoch.geometry.values[0]

df_outlines['area_2017'] = df_outlines.geometry.area

# some glaciers are from 2017 and some from 2018
bgndate = [f'{year}9999' for year in df_outlines.YEAR]
geodetic_mb_period = [f'2006-01-01_{year}-01-01' for year in df_outlines.YEAR]

# load additional data used for the model calibration
df_outlines_2006 = gpd.read_file(os.path.join(INPUT_DIR, Outline_2006_file))
df_dv_97_06 = pd.read_csv(os.path.join(INPUT_DIR, dv_97_06_file))
df_dv_06_17 = pd.read_csv(os.path.join(INPUT_DIR, dv_06_17_file))
df_volume = gpd.read_file(os.path.join(INPUT_DIR, Volume_file))
df_outline_hochjoch_2006 = gpd.read_file(os.path.join(INPUT_DIR,
                                                      Outline_Hochjoch_2006_file))
df_hochjoch_data = pd.read_csv(os.path.join(INPUT_DIR, dv_Hochjoch_file),
                               index_col=0)

dv_97_06 = []
dv_06_17 = []
dv_97_17 = []
area_2006 = []
volume_2006 = []
volume_2017 = []

for i, el in df_outlines.iterrows():
    glacier_id = el.nr
    # special treatment for hochjochferner
    if glacier_id == 2121:
        area_2006.append(df_outline_hochjoch_2006.geometry.area.values[0])
        dv_97_06.append(df_hochjoch_data.loc[2006].volchange)
        dv_06_17.append(df_hochjoch_data.loc[2017].volchange)
        dv_97_17.append(dv_97_06[-1] + dv_06_17[-1])
        volume_2006.append(df_hochjoch_data.loc[2006].volume)
    else:
        area_2006.append(df_outlines_2006[df_outlines_2006['nr'] == glacier_id].geometry.area.values[0])
        dv_97_06.append(df_dv_97_06[df_dv_97_06['nr'] == glacier_id]['dV'].values[0])
        dv_06_17.append(df_dv_06_17[df_dv_06_17['nr'] == glacier_id]['dV'].values[0])
        dv_97_17.append(dv_97_06[-1] + dv_06_17[-1])
    
        # some things are wrong in volume file
        if glacier_id == 2072:
            volume_tmp = df_volume[df_volume['Gletschern'] == 'Langtaler Ferner']
        elif glacier_id == 2125:
            # here I must adapt the code maybe for hintereisferner total
            volume_tmp = df_volume[df_volume['ID'] == 2125000]
        else:
            volume_tmp = df_volume[df_volume['ID'] == glacier_id]
        try:
            # sum is only needed for HEF Toteis and without Toteis
            volume_2006.append((volume_tmp['AREA'] * volume_tmp['Mean_h']).sum())
        except IndexError:
            print(glacier_id)
            #raise
    
    # calculate the volume 2017
    volume_2017.append(volume_2006[-1] + dv_06_17[-1])
    
df_outlines['area_2006'] = area_2006
df_outlines['dv_97_06'] = dv_97_06
df_outlines['dv_06_17'] = dv_06_17
df_outlines['dv_97_17'] = dv_97_17
df_outlines['mb_period'] = geodetic_mb_period
df_outlines['volume_06'] = volume_2006
df_outlines['volume_17'] = volume_2017

# test all values are correct
for i, wl in df_outlines.iterrows():
    glacier_id = el.nr

    # special treatment for Hochjochferner
    if glacier_id != 2121:
        # check area 2006
        ref_value = df_outlines_2006[df_outlines_2006['nr'] == glacier_id].geometry.area.values[0]
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['area_2006'].values[0]
        assert np.isclose(ref_value, actual_value)
        
        # check dv 97 06
        ref_value = df_dv_97_06[df_dv_97_06['nr'] == glacier_id]['dV'].values[0]
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['dv_97_06'].values[0]
        assert np.isclose(ref_value, actual_value)
        
        # check dv 06 17
        ref_value_dv_06_17 = df_dv_06_17[df_dv_06_17['nr'] == glacier_id]['dV'].values[0]
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['dv_06_17'].values[0]
        assert np.isclose(ref_value_dv_06_17, actual_value)
        
        # check volume
        # HEF toteis special treetment
        glacier_id_vol = glacier_id if glacier_id != 2125 else 2125000
        volume_tmp = df_volume[df_volume['ID'] == glacier_id_vol]
        ref_value_volume_06 = (volume_tmp['AREA'] * volume_tmp['Mean_h']).sum()
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['volume_06'].values[0]
        assert np.isclose(ref_value_volume_06, actual_value)
        
        # check volume 2017
        ref_value_volume_17 = ref_value_volume_06 + ref_value_dv_06_17
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['volume_17'].values[0]
        assert np.isclose(ref_value_volume_17, actual_value)
    else:
        # check area 2006
        ref_value = df_outline_hochjoch_2006.geometry.area.values[0]
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['area_2006'].values[0]
        assert np.isclose(ref_value, actual_value)
        
        # check dv 97 06
        ref_value = df_hochjoch_data.loc[2006].volchange
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['dv_97_06'].values[0]
        assert np.isclose(ref_value, actual_value)
        
        # check dv 06 17
        ref_value_dv_06_17 = df_hochjoch_data.loc[2017].volchange
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['dv_06_17'].values[0]
        assert np.isclose(ref_value_dv_06_17, actual_value)
        
        # check volume
        ref_value_volume_06 = df_hochjoch_data.loc[2006].volume
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['volume_06'].values[0]
        assert np.isclose(ref_value_volume_06, actual_value)
        
        # check volume 2017
        ref_value_volume_17 = ref_value_volume_06 + ref_value_dv_06_17
        actual_value = df_outlines[df_outlines['nr'] == glacier_id]['volume_17'].values[0]
        assert np.isclose(ref_value_volume_17, actual_value)

rgidf_simple = utils.cook_rgidf(df_outlines, o1_region='11', bgndate=bgndate,
                                ids=df_outlines.nr.values, 
                                assign_column_values={'nr': 'ID',
                                                      'area_2017': 'area_2017',
                                                      'Gletschern': 'Name',
                                                      'YEAR': 'year',
                                                      'area_2006': 'area_2006',
                                                      'dv_97_06': 'dv_97_06',
                                                      'dv_06_17': 'dv_06_17',
                                                      'dv_97_17': 'dv_97_17',
                                                      'mb_period': 'mb_period',
                                                      'volume_06': 'volume_06',
                                                      'volume_17': 'volume_17'})
rgidf_simple = rgidf_simple.sort_values('area_2017', ascending=False)

# store the dataframe for analysis
rgidf_simple.to_csv(os.path.join(WORKING_DIR,
                                 'provided_data.csv'))

# start with initialisation
gdirs = workflow.init_glacier_directories(rgidf_simple, reset=True, force=True)


# ---------------------------------------- gis prepro tasks ------------------------------------------------
workflow.execute_entity_task(tasks.define_glacier_region, gdirs)
workflow.execute_entity_task(gis.rasterio_glacier_mask, gdirs)
workflow.execute_entity_task(tasks.simple_glacier_masks, gdirs)
workflow.execute_entity_task(tasks.elevation_band_flowline, gdirs)
workflow.execute_entity_task(tasks.fixed_dx_elevation_band_flowline, gdirs)
workflow.execute_entity_task(tasks.compute_downstream_line, gdirs)
workflow.execute_entity_task(tasks.compute_downstream_bedshape, gdirs)
opath = os.path.join(cfg.PATHS['working_dir'], 'glacier_statistics_gis_prepro' + output_file_suffix + '.csv')
ds = utils.compile_glacier_statistics(gdirs, path=opath)


# ----------------------------------------- climate tasks --------------------------------------------------
fp_clim = os.path.join(INPUT_DIR, Clim_file)
fp_dem = os.path.join(INPUT_DIR, Clim_dem_file)

workflow.execute_entity_task(process_spartacus_data, gdirs,
                             fp_clim=fp_clim,
                             fp_dem=fp_dem)


# ------------------------------------ static mb calibration -----------------------------------------------
# help function to convert dV into dmdtda
def get_dmdtda(gdir, rgidf_simple):
    """convert volume change into kg m-2 yr-1, for the period 2006 - 2017/18"""
    extra_data = rgidf_simple.loc[rgidf_simple['RGIId'] == gdir.rgi_id]
    ice_density = 850  # kg m-3
    yr0 = 2006
    yr1 = int(extra_data.year.values[0])
    dV = float(rgidf_simple.loc[rgidf_simple['RGIId'] == gdir.rgi_id]['dv_06_17'].values)
    ref_mb = dV / gdir.rgi_area_m2 / (yr1 - yr0) * ice_density
    return ref_mb

workflow.execute_entity_task(
    tasks.mb_calibration_from_scalar_mb,
    [(gdir, {'ref_mb': get_dmdtda(gdir, rgidf_simple),
             'ref_period': rgidf_simple.loc[rgidf_simple['RGIId'] == gdir.rgi_id]['mb_period'].values[0],
             'temp_bias': temp_bias})
     for gdir in gdirs])
workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs)
opath = os.path.join(cfg.PATHS['working_dir'], 'glacier_statistics_climate_tasks' + output_file_suffix + '.csv')
ds = utils.compile_glacier_statistics(gdirs, path=opath)


# --------------------------- inversion and static initialisation -------------------------------------------
for gdir in gdirs:
    if gdir.rgi_id == 'RGI60-11.02125':
        cfg.PARAMS['min_slope'] = 10
    else:
        cfg.PARAMS['min_slope'] = 1.5
    workflow.calibrate_inversion_from_consensus(
        [gdir],
        volume_m3_reference=rgidf_simple[rgidf_simple.RGIId == gdir.rgi_id].volume_17.values[0],
        apply_fs_on_mismatch=True,
        error_on_mismatch=False)
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)
opath = os.path.join(cfg.PATHS['working_dir'], 'glacier_statistics_inversion' + output_file_suffix + '.csv')
ds = utils.compile_glacier_statistics(gdirs, path=opath)

# ----------------------------- fixed geometry initialisation -----------------------------------------------
start_year = 1979
workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,
                             fixed_geometry_spinup_yr=start_year,
                             output_filesuffix=f'_fixed_geometry_spinup_{start_year}')
opath = os.path.join(cfg.PATHS['working_dir'], 'fixed_historical_run_output' + output_file_suffix + '.nc')
utils.compile_run_output(gdirs, path=opath,
                         input_filesuffix=f'_fixed_geometry_spinup_{start_year}')

# --------------------------------- dynamic initialisation --------------------------------------------------
workflow.execute_entity_task(
    tasks.run_dynamic_melt_f_calibration,
    [(gdir,
      {'ys': start_year,
       'ref_dmdtda': get_dmdtda(gdir, rgidf_simple),
       'err_ref_dmdtda': 1.,  # kg m-2 yr-1, just assumed
       'ref_period': rgidf_simple[rgidf_simple.RGIId == gdir.rgi_id].mb_period.values[0],
       'ignore_errors': True,
       'melt_f_max_step_length_minimum': 0.1,
       'maxiter': 100,
       'target_yr': 2006,
       'output_filesuffix': f'_spinup_historical_{start_year}',
       'kwargs_run_function': {
           'target_value':
           rgidf_simple[rgidf_simple.RGIId == gdir.rgi_id].area_2006.values[0] * 1e-6,  # km2
           'spinup_start_yr_max': 1996,
           'first_guess_t_bias': -1,
           'maxiter': 100},
       'kwargs_fallback_function': {
           'target_value':
           rgidf_simple[rgidf_simple.RGIId == gdir.rgi_id].area_2006.values[0] * 1e-6,  # km2
           'spinup_start_yr_max': 1996,
           'first_guess_t_bias': -1,
           'maxiter': 100}})
      for gdir in gdirs])
opath = os.path.join(cfg.PATHS['working_dir'], f'glacier_statistics_dyn_mu' + output_file_suffix + '.csv')
ds = utils.compile_glacier_statistics(gdirs, path=opath)
opath = os.path.join(cfg.PATHS['working_dir'], 'spinup_historical_run_output' + output_file_suffix + '.nc')
utils.compile_run_output(gdirs, path=opath,
                         input_filesuffix=f'_spinup_historical_{start_year}')


# ---------------------------------- projection runs -------------------------------------------------------

# define CMIP6 GCMs and scenarios
all_GCM_ssps = [
    'ACCESS-CM2_ssp126',
    'ACCESS-CM2_ssp585',

    'ACCESS-ESM1-5_ssp126',
    'ACCESS-ESM1-5_ssp585',

    'BCC-CSM2-MR_ssp126',
    'BCC-CSM2-MR_ssp245',
    'BCC-CSM2-MR_ssp370',
    'BCC-CSM2-MR_ssp585',

    'CAMS-CSM1-0_ssp119',
    'CAMS-CSM1-0_ssp126',
    'CAMS-CSM1-0_ssp245',
    'CAMS-CSM1-0_ssp370',
    'CAMS-CSM1-0_ssp585',

    'CESM2_ssp126',
    'CESM2_ssp245',
    'CESM2_ssp370',
    'CESM2_ssp585',

    'CESM2-WACCM_ssp126',
    'CESM2-WACCM_ssp245',
    'CESM2-WACCM_ssp370',
    'CESM2-WACCM_ssp534-over',
    'CESM2-WACCM_ssp585',

    'CMCC-CM2-SR5_ssp245',
    'CMCC-CM2-SR5_ssp585',

    'CanESM5_ssp126',
    'CanESM5_ssp534-over',
    'CanESM5_ssp585',

    'EC-Earth3_ssp126',
    'EC-Earth3_ssp245',
    'EC-Earth3_ssp370',
    'EC-Earth3_ssp585',

    'EC-Earth3-Veg_ssp119',
    'EC-Earth3-Veg_ssp126',
    'EC-Earth3-Veg_ssp245',
    'EC-Earth3-Veg_ssp370',
    'EC-Earth3-Veg_ssp585',

    'FGOALS-f3-L_ssp126',
    'FGOALS-f3-L_ssp245',
    'FGOALS-f3-L_ssp370',
    'FGOALS-f3-L_ssp585',

    'GFDL-ESM4_ssp119',
    'GFDL-ESM4_ssp126',
    'GFDL-ESM4_ssp245',
    'GFDL-ESM4_ssp370',
    'GFDL-ESM4_ssp585',

    'INM-CM4-8_ssp126',
    'INM-CM4-8_ssp245',
    'INM-CM4-8_ssp370',
    'INM-CM4-8_ssp585',

    'INM-CM5-0_ssp126',
    'INM-CM5-0_ssp245',
    'INM-CM5-0_ssp370',
    'INM-CM5-0_ssp585',

    'IPSL-CM6A-LR_ssp126',
    'IPSL-CM6A-LR_ssp534-over',
    'IPSL-CM6A-LR_ssp585',
    
    'MPI-ESM1-2-HR_ssp126',
    'MPI-ESM1-2-HR_ssp245',
    'MPI-ESM1-2-HR_ssp370',
    'MPI-ESM1-2-HR_ssp585',

    'MRI-ESM2-0_ssp119',
    'MRI-ESM2-0_ssp126',
    'MRI-ESM2-0_ssp245',
    'MRI-ESM2-0_ssp370',
    'MRI-ESM2-0_ssp434',
    'MRI-ESM2-0_ssp460',
    'MRI-ESM2-0_ssp534-over',
    'MRI-ESM2-0_ssp585',
    
    'NorESM2-MM_ssp126',
    'NorESM2-MM_ssp245',
    'NorESM2-MM_ssp370',
    'NorESM2-MM_ssp585',

    'TaiESM1_ssp585',
    ]
#all_GCM = ['BCC-CSM2-MR',
#           'CESM2',
#           'CESM2-WACCM',
#           'EC-Earth3',
#           'EC-Earth3-Veg',
#           'FGOALS-f3-L',
#           'GFDL-ESM4',
#           'INM-CM4-8',
#           'INM-CM5-0',
#           'MPI-ESM1-2-HR',
#           'MRI-ESM2-0',
#           'NorESM2-MM'
#          ]
#
#all_ssp = ['ssp126', 'ssp245', 'ssp370', 'ssp585']

# download locations for precipitation and temperature
bp = 'https://cluster.klima.uni-bremen.de/~oggm/cmip6/GCM/{}/{}_{}_r1i1p1f1_pr.nc'
bt = 'https://cluster.klima.uni-bremen.de/~oggm/cmip6/GCM/{}/{}_{}_r1i1p1f1_tas.nc'

# 'download' and process GCM data, should be not downloaded (linked to local downloads in slurm file)
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')

    # 'Download' the files
    ft = utils.file_downloader(bt.format(GCM, GCM, ssp))
    fp = utils.file_downloader(bp.format(GCM, GCM, ssp))
    try:
        # bias correct them
        workflow.execute_entity_task(gcm_climate.process_cmip_data, gdirs, 
                                     filesuffix='_{}_{}'.format(GCM, ssp),  # recognize the climate file for later
                                     fpath_temp=ft,  # temperature projections
                                     fpath_precip=fp,  # precip projections
                                     #year_range=('2000', '2019'),
                                     #apply_bias_correction=True,
                                     )
    except ValueError:
        print('No ' + GCM + ' run with scenario ' + ssp + ' available')

# --------------------  actual projection runs using FIXED GEOMETRY SPINUP -----------------------------
#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')

    rid = '_{}_{}'.format(GCM, ssp)
    try:
        workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, 
                                     climate_filename='gcm_data',  # use gcm_data, not climate_historical
                                     climate_input_filesuffix=rid,  # use the chosen scenario
                                     init_model_filesuffix=f'_fixed_geometry_spinup_{start_year}',
                                     output_filesuffix=f'{rid}_fixed_spinup',  # recognize the run for later
                                    )
    except FileNotFoundError:
        print('No ' + GCM + ' run with scenario ' + ssp + ' available')

# add area using min h
#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')

    rid = '_{}_{}_fixed_spinup'.format(GCM, ssp)

    try:
        workflow.execute_entity_task(calculate_area_min_h, gdirs,
                                     filesuffix=rid,
                                     min_thickness=cfg.PARAMS['dynamic_spinup_min_ice_thick'])
    except FileNotFoundError:
        print('No ' + GCM + ' run with scenario ' + ssp + ' available')

# merge all runs into one nc file
ds_all = []
creation_date = strftime("%Y-%m-%d %H:%M:%S", gmtime())

#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')
    try:
        rid = '_{}_{}_fixed_spinup'.format(GCM, ssp)

        ds_tmp = utils.compile_run_output(gdirs, input_filesuffix=rid)
        ds_tmp.coords['GCM'] = GCM
        ds_tmp.coords['GCM'].attrs['description'] = 'used Global Circulation Model'
        ds_tmp = ds_tmp.expand_dims("GCM")
        ds_tmp.coords['SSP'] = ssp
        ds_tmp.coords['SSP'].attrs['description'] = 'used Shared Socioeconomic Pathway'
        ds_tmp = ds_tmp.expand_dims("SSP")
        ds_all.append(ds_tmp)
        ds_tmp.attrs['creation_date'] = creation_date
    except RuntimeError as err:
        if str(err) == 'Found no valid glaciers!':
            print(f'No data for GCM {GCM} with SSP {ssp} found!')
        else:
            raise RuntimeError(err)

ds_merged = xr.combine_by_coords(ds_all, fill_value=np.nan)
opath = os.path.join(cfg.PATHS['working_dir'], 'projections_run_output_fixed_spinup' + output_file_suffix + '.nc')
ds_merged.to_netcdf(opath)

# ------------------  actual projection runs using DYNAMIC SPINUP  -----------------------------------
#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')

    rid = '_{}_{}'.format(GCM, ssp)
    try:
        workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, 
                                     climate_filename='gcm_data',  # use gcm_data, not climate_historical
                                     climate_input_filesuffix=rid,  # use the chosen scenario
                                     init_model_filesuffix=f'_spinup_historical_{start_year}',
                                     output_filesuffix=rid,  # recognize the run for later
                                    )
    except FileNotFoundError:
        print('No ' + GCM + ' run with scenario ' + ssp + ' available')

# add area using min h
#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')

    rid = '_{}_{}'.format(GCM, ssp)

    try:
        workflow.execute_entity_task(calculate_area_min_h, gdirs,
                                     filesuffix=rid,
                                     min_thickness=cfg.PARAMS['dynamic_spinup_min_ice_thick'])
    except FileNotFoundError:
        print('No ' + GCM + ' run with scenario ' + ssp + ' available')

# merge all runs into one nc file
ds_all = []
creation_date = strftime("%Y-%m-%d %H:%M:%S", gmtime())

#for GCM in all_GCM:
#    for ssp in all_ssp:
for GCM_ssp in all_GCM_ssps:
    GCM, ssp = GCM_ssp.split('_')
    try:
        rid = '_{}_{}'.format(GCM, ssp)

        ds_tmp = utils.compile_run_output(gdirs, input_filesuffix=rid)
        ds_tmp.coords['GCM'] = GCM
        ds_tmp.coords['GCM'].attrs['description'] = 'used Global Circulation Model'
        ds_tmp = ds_tmp.expand_dims("GCM")
        ds_tmp.coords['SSP'] = ssp
        ds_tmp.coords['SSP'].attrs['description'] = 'used Shared Socioeconomic Pathway'
        ds_tmp = ds_tmp.expand_dims("SSP")
        ds_all.append(ds_tmp)
        ds_tmp.attrs['creation_date'] = creation_date
    except RuntimeError as err:
        if str(err) == 'Found no valid glaciers!':
            print(f'No data for GCM {GCM} with SSP {ssp} found!')
        else:
            raise RuntimeError(err)

ds_merged = xr.combine_by_coords(ds_all, fill_value=np.nan)
opath = os.path.join(cfg.PATHS['working_dir'], 'projections_run_output' + output_file_suffix + '.nc')
ds_merged.to_netcdf(opath)
