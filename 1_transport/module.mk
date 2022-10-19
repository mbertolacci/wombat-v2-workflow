$(shell mkdir -p 1_transport/intermediates)

CLIMATOLOGY_PARTIAL = 1_transport/src/partials/climatology.R
export CLIMATOLOGY_PARTIAL

1_TRANSPORT_TARGETS += 1_transport/intermediates/runs/RUNS 1_transport/intermediates/runs-r10-r15-rNZ/RUNS

# Flux inventories

INVENTORY_OUTPUT_YEARS = 2014 2015 2016 2017 2018 2019 2020 2021

## SiB4 inventory

SIB4_MONTHLY = 1_transport/intermediates/sib4-monthly.nc
SIB4_INPUT_YEARS = 2014 2015 2016 2017 2018 2019 2020
SIB4_HOURLY = $(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),1_transport/intermediates/sib4-hourly-$(SIB4_YEAR).nc)

SIB4_CLIMATOLOGY_ASSIM = 1_transport/intermediates/sib4-climatology-assim.nc
SIB4_CLIMATOLOGY_RESP_TOT = 1_transport/intermediates/sib4-climatology-resp-tot.nc
SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY = $(foreach SIB4_YEAR,$(INVENTORY_OUTPUT_YEARS),1_transport/intermediates/sib4-hourly-climatology-inventory-assim-$(SIB4_YEAR).nc)
SIB4_CLIMATOLOGY_INVENTORY_RESP_TOT_HOURLY = $(foreach SIB4_YEAR,$(INVENTORY_OUTPUT_YEARS),1_transport/intermediates/sib4-hourly-climatology-inventory-resp-tot-$(SIB4_YEAR).nc)
SIB4_RESIDUAL_ASSIM_HOURLY = $(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),1_transport/intermediates/sib4-hourly-residual-assim-$(SIB4_YEAR).nc)
SIB4_RESIDUAL_RESP_TOT_HOURLY = $(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),1_transport/intermediates/sib4-hourly-residual-resp-tot-$(SIB4_YEAR).nc)

$(SIB4_MONTHLY):
	cdo -v -z zip_6 \
		-setattribute,assim@units=kg/m2/s \
		-setattribute,resp_tot@units=kg/m2/s \
		-mulc,44.01 \
		-mulc,1e-9 \
		-monmean \
		-vertsum \
		-select,name=assim,resp_tot \
		data/sib4-hourly/sib4-hourly-* \
		$@

$(SIB4_HOURLY) &:
	cdo -w -z zip_6 \
		-splityear \
		-setattribute,assim@units=kg/m2/s \
		-setattribute,resp_tot@units=kg/m2/s \
		-shifttime,-30minutes \
		-settunits,minutes \
		-mulc,44.01 \
		-mulc,1e-9 \
		-vertsum \
		-select,name=assim,resp_tot \
		data/sib4-hourly/sib4-hourly-{2014,2015,2016,2017,2018,2019,2020}-* \
		1_transport/intermediates/sib4-hourly-

$(SIB4_CLIMATOLOGY_ASSIM): \
	1_transport/src/estimate-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_MONTHLY)
	Rscript $< \
		--input $(SIB4_MONTHLY) \
		--field-name assim \
		--time-resolution monthly \
		--time-alignment middle \
		--origin "2000-01-01 00:00:00" \
		--harmonics 3 \
		--output $@

$(SIB4_CLIMATOLOGY_RESP_TOT): \
	1_transport/src/estimate-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_MONTHLY)
	Rscript $< \
		--input $(SIB4_MONTHLY) \
		--field-name resp_tot \
		--time-resolution monthly \
		--time-alignment middle \
		--origin "2000-01-01 00:00:00" \
		--harmonics 3 \
		--output $@

$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY) &: \
	1_transport/src/expand-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_CLIMATOLOGY_ASSIM) \
	$(SIB4_HOURLY)
	Rscript $< \
		--year $(INVENTORY_OUTPUT_YEARS) \
		--split-year \
		--climatology $(SIB4_CLIMATOLOGY_ASSIM) \
		--field-name assim \
		--time-resolution hourly \
		--output $(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY)

$(SIB4_CLIMATOLOGY_INVENTORY_RESP_TOT_HOURLY) &: \
	1_transport/src/expand-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_CLIMATOLOGY_RESP_TOT) \
	$(SIB4_HOURLY)
	Rscript $< \
		--year $(INVENTORY_OUTPUT_YEARS) \
		--split-year \
		--climatology $(SIB4_CLIMATOLOGY_RESP_TOT) \
		--field-name resp_tot \
		--time-resolution hourly \
		--output $(SIB4_CLIMATOLOGY_INVENTORY_RESP_TOT_HOURLY)

$(SIB4_RESIDUAL_ASSIM_HOURLY) &: \
	1_transport/src/decompose-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_CLIMATOLOGY_ASSIM) \
	$(SIB4_HOURLY)
	Rscript $< \
		--input $(SIB4_HOURLY) \
		--climatology $(SIB4_CLIMATOLOGY_ASSIM) \
		--field-name assim \
		--time-resolution hourly \
		--time-alignment start \
		--output $(SIB4_RESIDUAL_ASSIM_HOURLY)

$(SIB4_RESIDUAL_RESP_TOT_HOURLY) &: \
	1_transport/src/decompose-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(SIB4_CLIMATOLOGY_RESP_TOT) \
	$(SIB4_HOURLY)
	Rscript $< \
		--input $(SIB4_HOURLY) \
		--climatology $(SIB4_CLIMATOLOGY_RESP_TOT) \
		--field-name resp_tot \
		--time-resolution hourly \
		--time-alignment start \
		--output $(SIB4_RESIDUAL_RESP_TOT_HOURLY)

## Oceans inventories

LANDSCHUTZER_INVENTORY_INPUT = data/spco2_MPI-SOM_FFN_v2020.nc
LANDSCHUTZER_INVENTORY = 1_transport/intermediates/spco2_MPI-SOM_FFN_v2020_raw.nc
LANDSCHUTZER_CLIMATOLOGY = 1_transport/intermediates/spco2_MPI-SOM_FFN_v2020_climatology.nc
LANDSCHUTZER_CLIMATOLOGY_INVENTORY = 1_transport/intermediates/spco2_MPI-SOM_FFN_v2020_climatology_inventory.nc
LANDSCHUTZER_INVENTORY_RESIDUAL = 1_transport/intermediates/spco2_MPI-SOM_FFN_v2020_residual.nc

$(LANDSCHUTZER_INVENTORY): \
	1_transport/src/process-ocean.sh \
	$(LANDSCHUTZER_INVENTORY_INPUT)
	bash $< $(LANDSCHUTZER_INVENTORY_INPUT) $@

$(LANDSCHUTZER_CLIMATOLOGY): \
	1_transport/src/estimate-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(LANDSCHUTZER_INVENTORY)
	Rscript $< \
		--input $(LANDSCHUTZER_INVENTORY) \
		--field-name fgco2_raw \
		--time-resolution monthly \
		--time-alignment start \
		--origin "2000-01-01 00:00:00" \
		--harmonics 2 \
		--output $@

$(LANDSCHUTZER_CLIMATOLOGY_INVENTORY): \
	1_transport/src/expand-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(LANDSCHUTZER_CLIMATOLOGY) \
	$(LANDSCHUTZER_INVENTORY)
	Rscript $< \
		--year $(INVENTORY_OUTPUT_YEARS) \
		--climatology $(LANDSCHUTZER_CLIMATOLOGY) \
		--field-name fgco2_raw \
		--time-resolution monthly \
		--output $@

$(LANDSCHUTZER_INVENTORY_RESIDUAL): \
	1_transport/src/decompose-climatology.R \
	$(CLIMATOLOGY_PARTIAL) \
	$(LANDSCHUTZER_CLIMATOLOGY) \
	$(LANDSCHUTZER_INVENTORY)
	Rscript $< \
		--input $(LANDSCHUTZER_INVENTORY) \
		--climatology $(LANDSCHUTZER_CLIMATOLOGY) \
		--field-name fgco2_raw \
		--time-resolution monthly \
		--time-alignment start \
		--output $@

# Runs

## Basis functions

TRANSCOM_MASK_ORIGINAL_1X1 = data/TRANSCOM_mask_original_1x1.nc
GEOS_CHEM_RESTART_FILE = data/GEOSChem.Restart.20140901_0000z.nc4
GEOS_2X25_GRID = data/geos.2x25.grid
FOSSIL_INVENTORY_PATH = data/fossil-mipv10
FIRES_INVENTORY_PATH = data/GFED4.1s_Aleya_hemco_080321
GEOS_CHEM_TEMPLATE_DIRECTORY = 1_transport/src/GEOS_Chem

TRANSCOM_MASK_GEOS_2X25 = 1_transport/intermediates/TRANSCOM_mask_GEOS_Chem_2x2.5.nc
REGION_MASK = 1_transport/intermediates/region-mask.rds
REGION_MASK_BASE = 1_transport/intermediates/region-mask-base.rds

1_transport/intermediates/runs/RUNS: \
	1_transport/src/basis-runs.R \
	$(GEOS_CHEM_TEMPLATE_DIRECTORY)/* \
	$(LANDSCHUTZER_CLIMATOLOGY_INVENTORY) \
	$(LANDSCHUTZER_INVENTORY_RESIDUAL) \
	$(SIB4_RESIDUAL_ASSIM_HOURLY) \
	$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY) \
	$(SIB4_RESIDUAL_RESP_TOT_HOURLY) \
	$(SIB4_CLIMATOLOGY_INVENTORY_RESP_TOT_HOURLY) \
	$(REGION_MASK) \
	$(GEOS_CHEM_RESTART_FILE)
	Rscript $< \
		--template-directory $(GEOS_CHEM_TEMPLATE_DIRECTORY) \
		--region-mask $(REGION_MASK_BASE) \
		--fossil-inventory-path $(FOSSIL_INVENTORY_PATH) \
		--fires-inventory-path $(FIRES_INVENTORY_PATH) \
		--ocean-climatology-inventory-path $(LANDSCHUTZER_CLIMATOLOGY_INVENTORY) \
		--ocean-residual-path $(LANDSCHUTZER_INVENTORY_RESIDUAL) \
		--sib4-assim-climatology-inventory-path '1_transport/intermediates/sib4-hourly-climatology-inventory-assim-$$YYYY.nc' \
		--sib4-assim-residual-path '1_transport/intermediates/sib4-hourly-residual-assim-$$YYYY.nc' \
		--sib4-resp-tot-climatology-inventory-path '1_transport/intermediates/sib4-hourly-climatology-inventory-resp-tot-$$YYYY.nc' \
		--sib4-resp-tot-residual-path '1_transport/intermediates/sib4-hourly-residual-resp-tot-$$YYYY.nc' \
		--geos-chem-restart $(GEOS_CHEM_RESTART_FILE) \
		--geos-chem-extdata $(GEOS_CHEM_EXTDATA) \
		--geos-chem-source $(GEOS_CHEM_SOURCE) \
		--output-directory 1_transport/intermediates/runs

1_transport/intermediates/runs-r10-r15-rNZ/RUNS: \
	1_transport/src/basis-runs.R \
	$(GEOS_CHEM_TEMPLATE_DIRECTORY)/* \
	$(LANDSCHUTZER_CLIMATOLOGY_INVENTORY) \
	$(LANDSCHUTZER_INVENTORY_RESIDUAL) \
	$(SIB4_RESIDUAL_ASSIM_HOURLY) \
	$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY) \
	$(SIB4_RESIDUAL_RESP_TOT_HOURLY) \
	$(SIB4_CLIMATOLOGY_INVENTORY_RESP_TOT_HOURLY) \
	$(REGION_MASK) \
	$(GEOS_CHEM_RESTART_FILE)
	Rscript $< \
		--template-directory $(GEOS_CHEM_TEMPLATE_DIRECTORY) \
		--region-mask $(REGION_MASK) \
		--region Region10 Region15 RegionNZ \
		--fossil-inventory-path $(FOSSIL_INVENTORY_PATH) \
		--fires-inventory-path $(FIRES_INVENTORY_PATH) \
		--ocean-climatology-inventory-path $(LANDSCHUTZER_CLIMATOLOGY_INVENTORY) \
		--ocean-residual-path $(LANDSCHUTZER_INVENTORY_RESIDUAL) \
		--sib4-assim-climatology-inventory-path '1_transport/intermediates/sib4-hourly-climatology-inventory-assim-$$YYYY.nc' \
		--sib4-assim-residual-path '1_transport/intermediates/sib4-hourly-residual-assim-$$YYYY.nc' \
		--sib4-resp-tot-climatology-inventory-path '1_transport/intermediates/sib4-hourly-climatology-inventory-resp-tot-$$YYYY.nc' \
		--sib4-resp-tot-residual-path '1_transport/intermediates/sib4-hourly-residual-resp-tot-$$YYYY.nc' \
		--geos-chem-restart $(GEOS_CHEM_RESTART_FILE) \
		--geos-chem-extdata $(GEOS_CHEM_EXTDATA) \
		--geos-chem-source $(GEOS_CHEM_SOURCE) \
		--output-directory 1_transport/intermediates/runs-r10-r15-rNZ

validate_runs: \
	1_transport/src/validate-runs.R
	Rscript $< \
		--runs 1_transport/intermediates/runs


$(REGION_MASK): \
	1_transport/src/region-mask.R \
	$(TRANSCOM_MASK_GEOS_2X25)
	Rscript $< \
		--transcom-mask $(TRANSCOM_MASK_GEOS_2X25) \
		--output $@

$(REGION_MASK_BASE): \
	1_transport/src/region-mask-base.R \
	$(TRANSCOM_MASK_GEOS_2X25)
	Rscript $< \
		--transcom-mask $(TRANSCOM_MASK_GEOS_2X25) \
		--output $@

$(TRANSCOM_MASK_GEOS_2X25): \
	1_transport/src/transcom-mask.sh \
	$(TRANSCOM_MASK_ORIGINAL_1X1) \
	$(GEOS_2X25_GRID)
	bash $< \
		$(TRANSCOM_MASK_ORIGINAL_1X1) \
		$(GEOS_2X25_GRID) \
		$@
