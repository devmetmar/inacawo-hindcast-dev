#!/bin/bash

cd /home/cawofcst_dev/devmetmar/cawo-anal-202412

echo "Linking binary and coupler"
ln -sf /home/cawofcst_dev/devmetmar/cawo_3way_build/coawstM .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/scrip_mar2023.nc .

echo "Linking static files"
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/varinfo.dat.swell .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/cawo_grid_3km_min3_max5000_wrf_bmkgmask_sponge5_v2.6.7.nc .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/cawo_g2_tide_18581117_TPXO9altasv5.nc .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan_coord_v255.grd .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan_bathy_v255.bot .

echo "Linking templates"
cd /home/cawofcst_dev/devmetmar/cawo-anal-202412/templates
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/coupling_cawo.in.template_valid3b .
ln -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/swan.in.template_gfswave_HHMMSS.maxiter3.cold.10day .
ls -sf /scratch/cawofcst_dev/data/cawo_models_run/anal_cycles/cawo_3way_template/roms.in.template_valid3b.06hrs.ROMS_GRD_VERS.swell .



