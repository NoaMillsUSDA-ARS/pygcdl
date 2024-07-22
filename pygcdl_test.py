import pygcdl

# Create pygcdl object
pygcdl_obj = pygcdl.PyGeoCDL(url_base="http://127.0.0.1:8000")

# Test list_datasets call
def test_list_datasets():
    r = pygcdl_obj.list_datasets()
    expected_response = {'DaymetV4': 'Daymet Version 4', 'GTOPO30': 'Global 30 Arc-Second Elevation', 'MODIS_NDVI': 'MODIS NDVI Data, Smoothed and Gap-filled, for the Conterminous US: 2000-2015', 'NASS_CDL': 'NASS Cropland Data Layer', 'NLCD': 'National Land Cover Database', 'PRISM': 'PRISM', 'RAPV3': 'Rangeland Analysis Platform Version 3', 'SMAP-HB1km': 'SMAP HydroBlocks - 1 km', 'Soilgrids250mV2': 'SoilGrids — global gridded soil information', 'VIP': 'Vegetation Index and Phenology (VIP) Vegetation Indices Daily Global 0.05Deg CMG V004'}
    assert r == expected_response
