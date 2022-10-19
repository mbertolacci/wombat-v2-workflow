4_RESULTS_SRC_DIR = 4_results/src
4_RESULTS_INTERMEDIATES_DIR = 4_results/intermediates
4_RESULTS_FIGURES_DIR = 4_results/figures
4_RESULTS_PRODUCTS_DIR = 4_results/products

$(shell mkdir -p $(4_RESULTS_INTERMEDIATES_DIR))
$(shell mkdir -p $(4_RESULTS_FIGURES_DIR))
$(shell mkdir -p $(4_RESULTS_PRODUCTS_DIR))

# Intermediates
PERTURBATIONS_AUGMENTED = $(4_RESULTS_INTERMEDIATES_DIR)/perturbations-augmented.fst
LANDSCHUTZER_CLIMATOLOGY_2X25 = $(4_RESULTS_INTERMEDIATES_DIR)/spco2_MPI-SOM_FFN_v2020_climatology-2x25.nc
SIB4_CLIMATOLOGY_ASSIM_2X25 = $(4_RESULTS_INTERMEDIATES_DIR)/sib4-climatology-assim-2x25.nc
SIB4_CLIMATOLOGY_RESP_TOT_2X25 = $(4_RESULTS_INTERMEDIATES_DIR)/sib4-climatology-resp-tot-2x25.nc
CLIMATOLOGY_BY_REGION = $(4_RESULTS_INTERMEDIATES_DIR)/climatology-by-region.rds
TREND_GRID = $(4_RESULTS_INTERMEDIATES_DIR)/trend-grid.fst
REGION_SF = $(4_RESULTS_INTERMEDIATES_DIR)/region-sf.rds
REGION_GRID = $(4_RESULTS_INTERMEDIATES_DIR)/region-grid.rds
SIX_YEAR_AVERAGE = $(4_RESULTS_INTERMEDIATES_DIR)/six-year-average.fst

4_RESULTS_TARGETS += \
	$(4_RESULTS_PRODUCTS_DIR)/WOMBAT_v2_CO2_gridded_flux_samples.nc4 \
	$(4_RESULTS_PRODUCTS_DIR)/WOMBAT_v2_CO2_gridded_climatology_samples.nc4 \
	$(4_RESULTS_FIGURES_DIR)/average-and-trend-map.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-global-zonal.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-land-transcoms.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-global-zonal-table.txt \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-land-transcoms-table.txt \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-nee.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-gpp.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-resp.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-nee-2.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-gpp-2.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-resp-2.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-nee-3.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-gpp-3.pdf \
	$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-resp-3.pdf \
	$(4_RESULTS_FIGURES_DIR)/trend-table.txt \
	$(4_RESULTS_FIGURES_DIR)/flux-components-global.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-n-boreal.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-n-temperate.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-n-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-s-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-components-s-extratropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-global.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-n-boreal.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-n-temperate.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-n-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-s-tropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/flux-net-s-extratropical.pdf \
	$(4_RESULTS_FIGURES_DIR)/el-nino-2015-residual-map.pdf \
	$(4_RESULTS_FIGURES_DIR)/global-bottom-up.pdf \
	$(4_RESULTS_FIGURES_DIR)/landschutzer-fluxes.pdf \
	$(4_RESULTS_FIGURES_DIR)/region-map.pdf \
	$(4_RESULTS_FIGURES_DIR)/sib4-fluxes.pdf \
	$(4_RESULTS_FIGURES_DIR)/observation-map.pdf \
	$(4_RESULTS_FIGURES_DIR)/effective-sample-size.txt \
	$(4_RESULTS_FIGURES_DIR)/hyperparameter-table.tex \
	$(4_RESULTS_FIGURES_DIR)/ell-table.tex \
	$(4_RESULTS_FIGURES_DIR)/traceplots.pdf

## Products

$(4_RESULTS_PRODUCTS_DIR)/WOMBAT_v2_CO2_gridded_flux_samples.nc4: \
	$(4_RESULTS_SRC_DIR)/gridded-flux-samples.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples $(SAMPLES_LNLGIS) \
		--output $@

$(4_RESULTS_PRODUCTS_DIR)/WOMBAT_v2_CO2_gridded_climatology_samples.nc4: \
	$(4_RESULTS_SRC_DIR)/gridded-climatology-samples.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--region-grid $(REGION_GRID) \
		--samples $(SAMPLES_LNLGIS) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
		--landschutzer-climatology $(LANDSCHUTZER_CLIMATOLOGY_2X25) \
		--output $@

## Figures

$(4_RESULTS_FIGURES_DIR)/average-and-trend-map.pdf: \
	$(4_RESULTS_SRC_DIR)/average-and-trend-map.R \
	$(TREND_GRID) \
	$(SIX_YEAR_AVERAGE) \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--climatology-trend $(TREND_GRID) \
		--six-year-average $(SIX_YEAR_AVERAGE) \
		--region-sf $(REGION_SF) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/harmonic-shift-%-table.txt: \
	$(4_RESULTS_SRC_DIR)/harmonic-shift-table.R \
	$(CLIMATOLOGY_BY_REGION) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(BASE_PARTIAL)
	Rscript $< \
		--regions $* \
		--climatology-by-region $(CLIMATOLOGY_BY_REGION) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/harmonic-shift-%.pdf: \
	$(4_RESULTS_SRC_DIR)/harmonic-shift.R \
	$(CLIMATOLOGY_BY_REGION) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--regions $* \
		--region-sf $(REGION_SF) \
		--climatology-by-region $(CLIMATOLOGY_BY_REGION) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/harmonic-shift-map-%.pdf: \
	$(4_RESULTS_SRC_DIR)/harmonic-shift-map.R \
	$(CONTROL_EMISSIONS) \
	$(REGION_SF) \
	$(REGION_GRID) \
	$(SAMPLES_LNLGIS) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--region-sf $(REGION_SF) \
		--region-grid $(REGION_GRID) \
		--samples $(SAMPLES_LNLGIS) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
		--inventory-harmonic $* \
		--output $@

$(4_RESULTS_FIGURES_DIR)/trend-table.txt: \
	$(4_RESULTS_SRC_DIR)/trend-table.R \
	$(CLIMATOLOGY_BY_REGION) \
	$(BASE_PARTIAL)
	Rscript $< \
		--climatology-by-region $(CLIMATOLOGY_BY_REGION) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/flux-components-%.pdf: \
	$(4_RESULTS_SRC_DIR)/flux-components.R \
	$(AREA_1X1) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--area-1x1 $(AREA_1X1) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples $(SAMPLES_LNLGIS) \
		--region $* \
		--output $@

$(4_RESULTS_FIGURES_DIR)/flux-components-%.pdf: \
	$(4_RESULTS_SRC_DIR)/flux-components.R \
	$(AREA_1X1) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--area-1x1 $(AREA_1X1) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples $(SAMPLES_LNLGIS) \
		--region $* \
		--output $@

$(4_RESULTS_FIGURES_DIR)/flux-net-%.pdf: \
	$(4_RESULTS_SRC_DIR)/flux-net.R \
	$(AREA_1X1) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--area-1x1 $(AREA_1X1) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples $(SAMPLES_LNLGIS) \
		--region $* \
		--output $@

$(4_RESULTS_FIGURES_DIR)/el-nino-2015-residual-map.pdf: \
	$(4_RESULTS_SRC_DIR)/el-nino-2015-residual-map.R \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS_AUGMENTED) \
	$(REGION_SF) \
	$(SAMPLES_LNLGIS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--region-sf $(REGION_SF) \
		--samples $(SAMPLES_LNLGIS) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/global-bottom-up.pdf: \
	$(4_RESULTS_SRC_DIR)/global-bottom-up.R \
	$(CONTROL_EMISSIONS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/landschutzer-fluxes.pdf: \
	$(4_RESULTS_SRC_DIR)/landschutzer-fluxes.R \
	$(AREA_1X1) \
	$(LANDSCHUTZER_CLIMATOLOGY) \
	$(LANDSCHUTZER_INVENTORY) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--area-1x1 $(AREA_1X1) \
		--landschutzer-climatology $(LANDSCHUTZER_CLIMATOLOGY) \
		--landschutzer $(LANDSCHUTZER_INVENTORY) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/region-map.pdf: \
	$(4_RESULTS_SRC_DIR)/region-map.R \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--region-sf $(REGION_SF) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/sib4-fluxes.pdf: \
	$(4_RESULTS_SRC_DIR)/sib4-fluxes.R \
	$(AREA_1X1) \
	$(SIB4_CLIMATOLOGY_ASSIM) \
	$(SIB4_CLIMATOLOGY_RESP_TOT) \
	$(SIB4_MONTHLY) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--area-1x1 $(AREA_1X1) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT) \
		--sib4-monthly $(SIB4_MONTHLY) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/observation-map.pdf: \
	$(4_RESULTS_SRC_DIR)/observation-map.R \
	$(OBSERVATIONS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/traceplots.pdf: \
	$(4_RESULTS_SRC_DIR)/traceplots.R \
	$(SAMPLES_LNLGIS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--samples $(SAMPLES_LNLGIS) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/effective-sample-size.txt: \
	$(4_RESULTS_SRC_DIR)/effective-sample-size.R \
	$(SAMPLES_LNLGIS) \
	$(BASE_PARTIAL)
	Rscript $< \
		--samples $(SAMPLES_LNLGIS) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/ell-table.tex: \
	$(4_RESULTS_SRC_DIR)/ell-table.R \
	$(HYPERPARAMETER_ESTIMATES) \
	$(BASE_PARTIAL)
	Rscript $< \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--output $@

$(4_RESULTS_FIGURES_DIR)/hyperparameter-table.tex: \
	$(4_RESULTS_SRC_DIR)/hyperparameter-table.R \
	$(SAMPLES_LNLGIS) \
	$(BASE_PARTIAL)
	Rscript $< \
		--samples $(SAMPLES_LNLGIS) \
		--output $@

## Intermediates

$(SIX_YEAR_AVERAGE): \
	$(4_RESULTS_SRC_DIR)/six-year-average.R \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGIS) \
	$(UTILS_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--samples $(SAMPLES_LNLGIS) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--output $@

$(TREND_GRID): \
	$(4_RESULTS_SRC_DIR)/trend-grid.R \
	$(CONTROL_EMISSIONS) \
	$(REGION_GRID) \
	$(SAMPLES_LNLGIS) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	$(UTILS_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--region-grid $(REGION_GRID) \
		--samples $(SAMPLES_LNLGIS) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
		--output $@

$(CLIMATOLOGY_BY_REGION): \
	$(4_RESULTS_SRC_DIR)/climatology-by-region.R \
	$(CONTROL_EMISSIONS) \
	$(REGION_GRID) \
	$(SAMPLES_LNLGIS) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	$(LANDSCHUTZER_CLIMATOLOGY_2X25) \
	$(UTILS_PARTIAL)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--region-grid $(REGION_GRID) \
		--samples $(SAMPLES_LNLGIS) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
		--landschutzer-climatology $(LANDSCHUTZER_CLIMATOLOGY_2X25) \
		--output $@

$(PERTURBATIONS_AUGMENTED): \
	$(4_RESULTS_SRC_DIR)/perturbations-augmented.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(REGION_SF): \
	$(4_RESULTS_SRC_DIR)/region-sf.R \
	$(REGION_GRID)
	Rscript $< \
		--region-grid $(REGION_GRID) \
		--output $@

$(REGION_GRID): \
	$(4_RESULTS_SRC_DIR)/region-grid.R \
	$(TRANSCOM_MASK_GEOS_2X25)
	Rscript $< \
		--transcom-grid $(TRANSCOM_MASK_GEOS_2X25) \
		--output $@

$(LANDSCHUTZER_CLIMATOLOGY_2X25): \
	$(GEOS_2X25_GRID) \
	$(LANDSCHUTZER_CLIMATOLOGY)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(LANDSCHUTZER_CLIMATOLOGY) $@
	ncks -A -v variable $(LANDSCHUTZER_CLIMATOLOGY) $@

$(SIB4_CLIMATOLOGY_ASSIM_2X25): \
	$(GEOS_2X25_GRID) \
	$(SIB4_CLIMATOLOGY_ASSIM)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(SIB4_CLIMATOLOGY_ASSIM) $@
	ncks -A -v variable $(SIB4_CLIMATOLOGY_ASSIM) $@

$(SIB4_CLIMATOLOGY_RESP_TOT_2X25): \
	$(GEOS_2X25_GRID) \
	$(SIB4_CLIMATOLOGY_RESP_TOT)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(SIB4_CLIMATOLOGY_RESP_TOT) $@
	ncks -A -v variable $(SIB4_CLIMATOLOGY_RESP_TOT) $@
