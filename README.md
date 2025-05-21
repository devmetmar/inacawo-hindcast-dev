# CAWO Hindcast investigation

## Input data
- ROMS
    - [Mercator Ocean 1/12°](https://doi.org/10.48670/moi-00021) | [HYCOM 1/12°](https://tds.hycom.org/thredds/catalog.html)
    - [ERA5 0.25°](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)
    - TPXO (/scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/cawo_g2_tide_18581117_TPXO9altasv5.nc)
    - varinfo.dat.swell (/scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/varinfo.dat.swell)
    - Grid: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/cawo_grid_3km_min3_max5000_wrf_bmkgmask_sponge5_v2.6.7.nc
    - Namelist: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/roms.in.template_valid3b.06hrs.ROMS_GRD_VERS.swell
- SWAN
    - [Mercator Wave 0.2°](https://doi.org/10.48670/moi-00022)
    - [ERA5 0.25°](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)
    - Grid: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan_coord_v255.grd
    - Bathymetry: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan_bathy_v255.bot
    - Namelist: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan.in.template_gfswave_HHMMSS.maxiter3.cold.10day
- MCT Coupler 
    - SCRIP: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/scrip_mar2023.nc
    - Namelist: /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/coupling_cawo.in.template_valid3b