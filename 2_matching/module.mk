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
	Rscript -e "rmarkdown::render('$<', output_file = '$(shell basename $@)', output_dir = '$(shell dirname $@)', knit_root_dir = '../../', params = list(oco2_soundings = '$(OCO2_OBSERVATIONS)', tccon_sounding_directory = '$(TCCON_SOUNDING_DIRECTORY)', obspack_directory = '$(OBSPACK_DIRECTORY)', base_path = '2_matching/intermediates/runs/base'))"
