$(shell mkdir -p 2_matching/intermediates)
$(shell mkdir -p 2_matching/figures)

2_matching: 2_matching/figures/validate-base.pdf

2_matching/intermediates/runs/RUNS: \
	2_matching/src/prepare-matching.R
	Rscript $< \
		--input 1_transport/intermediates/runs \
		--oco2-soundings $(OCO2_OBSERVATIONS) \
		--tccon-sounding-directory $(TCCON_SOUNDING_DIRECTORY) \
		--obspack-directory $(OBSPACK_DIRECTORY) \
		--lauder $(LAUDER_DATA) \
		--template-directory 2_matching/src/templates \
		--output 2_matching/intermediates/runs

2_matching/intermediates/runs-r10-r15-rNZ/RUNS: \
	2_matching/src/prepare-matching.R
	Rscript $< \
		--input 1_transport/intermediates/runs-r10-r15-rNZ \
		--oco2-soundings $(OCO2_OBSERVATIONS) \
		--tccon-sounding-directory $(TCCON_SOUNDING_DIRECTORY) \
		--obspack-directory $(OBSPACK_DIRECTORY) \
		--lauder $(LAUDER_DATA) \
		--template-directory 2_matching/src/templates \
		--output 2_matching/intermediates/runs-r10-r15-rNZ

2_matching/figures/validate-base.pdf: 2_matching/src/validate-base.Rmd
	Rscript -e "rmarkdown::render('$<', output_file = '$(shell basename $@)', output_dir = '$(shell dirname $@)', knit_root_dir = '../../', params = list(oco2_soundings = 'data/OCO2_b10c_10sec_GOOD_r5.nc4', tccon_sounding_directory = 'data/downloaded_20210708', obspack_directory = 'data/obspack_co2_1_OCO2MIP_v3.2_2021-05-20/data/daily', base_path = '2_matching/intermediates/runs/base'))"
