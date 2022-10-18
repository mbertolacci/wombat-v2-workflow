$(shell mkdir -p 3_inversion/intermediates)

CONTROL_EMISSIONS = 3_inversion/intermediates/control-emissions.fst
PERTURBATIONS = 3_inversion/intermediates/perturbations.fst
OBSERVATIONS = 3_inversion/intermediates/observations.fst
BASIS_VECTORS = 3_inversion/intermediates/basis-vectors.fst
CONSTRAINTS = 3_inversion/intermediates/constraints.rds
PRIOR = 3_inversion/intermediates/prior.rds
SENSITIVITIES_BASE_PART_1 = 3_inversion/intermediates/sensitivities-base-oco2-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-0-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-2-part-1.fst \
	3_inversion/intermediates/sensitivities-base-tccon-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-base-oco2-daily-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-daily-assim-0-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-daily-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-daily-assim-2-part-1.fst \
	3_inversion/intermediates/sensitivities-base-tccon-daily-part-1.fst
SENSITIVITIES_R10_R15_RNZ_PART_1 = 3_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-0-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-2-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-0-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-2-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-daily-part-1.fst

H_IS = 3_inversion/intermediates/H-IS.mat.lz4
H_LNLG = 3_inversion/intermediates/H-LNLG.mat.lz4

RESIDUAL_1ST_STAGE = 3_inversion/intermediates/residual-1st-stage.fst
HYPERPARAMETER_ESTIMATES = 3_inversion/intermediates/hyperparameter-estimates.fst

SAMPLES_IS = 3_inversion/intermediates/samples-IS.rds
SAMPLES_LNLG = 3_inversion/intermediates/samples-LNLG.rds
SAMPLES_LNLGIS = 3_inversion/intermediates/samples-LNLGIS.rds

$(SAMPLES_IS): \
	3_inversion/src/samples.R \
  	$(OBSERVATIONS) \
  	$(BASIS_VECTORS) \
  	$(HYPERPARAMETER_ESTIMATES) \
  	$(CONSTRAINTS) \
  	$(PRIOR) \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
    $(H_IS)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode IS \
		--control \
		    2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--component-name IS \
		--component-parts IS \
		--component-transport-matrix $(H_IS) \
		--output $@

$(SAMPLES_LNLG): \
	3_inversion/src/samples.R \
  	$(OBSERVATIONS) \
  	$(BASIS_VECTORS) \
  	$(HYPERPARAMETER_ESTIMATES) \
  	$(CONSTRAINTS) \
  	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
    $(H_LNLG)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
		--component-name LNLG \
		--component-parts "LN|LG" \
		--component-transport-matrix $(H_LNLG) \
		--output $@

$(SAMPLES_LNLGIS): \
	3_inversion/src/samples.R \
  	$(OBSERVATIONS) \
  	$(BASIS_VECTORS) \
  	$(HYPERPARAMETER_ESTIMATES) \
  	$(CONSTRAINTS) \
  	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
    $(H_LNLG) \
    $(H_IS)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG IS \
		--control \
		    2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		    2_matching/intermediates/runs/base/oco2-hourly.fst \
		--component-name LNLG IS \
		--component-parts "LN|LG" IS \
		--component-transport-matrix $(H_LNLG) $(H_IS) \
		--output $@

# Hyperparameter estimates
$(HYPERPARAMETER_ESTIMATES): \
	3_inversion/src/hyperparameter-estimates.R \
  	$(OBSERVATIONS) \
  	$(RESIDUAL_1ST_STAGE)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--residuals $(RESIDUAL_1ST_STAGE) \
		--output $@

## Inversions (1st stage) to get a residual for hyperparameter estimates
$(RESIDUAL_1ST_STAGE): \
	3_inversion/src/residual.R \
    2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
    2_matching/intermediates/runs/base/oco2-hourly.fst \
  	$(OBSERVATIONS) \
  	$(BASIS_VECTORS) \
  	$(PRIOR) \
	$(CONSTRAINTS) \
    $(H_IS) \
    $(H_LNLG)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--overall-observation-mode LN LG IS \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--component-name LNLG IS \
		--component-parts "LN|LG" IS \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
		--output $@

## Transport matrix (H)

TRANSPORT_MATRIX_DEPS = 3_inversion/src/transport-matrix.R \
	3_inversion/intermediates/sensitivities-base-oco2-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-base-oco2-daily-part-1.fst \
	3_inversion/intermediates/sensitivities-base-obspack-daily-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-base-oco2-hourly-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily-part-1.fst \
	3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1-part-1.fst \
	$(BASIS_VECTORS) \
	$(OBSERVATIONS)
TRANSPORT_MATRIX_CALL = Rscript 3_inversion/src/transport-matrix.R \
	--basis-vectors $(BASIS_VECTORS) \
	--observations $(OBSERVATIONS) \
	--sensitivities-base \
		3_inversion/intermediates/sensitivities-base \
		3_inversion/intermediates/sensitivities-r10-r15-rNZ

$(H_IS): \
	$(TRANSPORT_MATRIX_DEPS)
	$(TRANSPORT_MATRIX_CALL) \
		--overall-observation-mode IS \
		--control 2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--output $@

$(H_LNLG): \
	$(TRANSPORT_MATRIX_DEPS)
	$(TRANSPORT_MATRIX_CALL) \
		--overall-observation-mode LN LG \
		--control 2_matching/intermediates/runs/base/oco2-hourly.fst \
		--output $@

## Sensitivities

$(SENSITIVITIES_R10_R15_RNZ_PART_1) &: \
	3_inversion/src/sensitivities.R
	Rscript $< \
		--input \
			oco2-hourly.fst \
			obspack-hourly-assim-1.fst \
			obspack-hourly-assim-2.fst \
			tccon-hourly.fst \
			oco2-daily.fst \
			obspack-daily-assim-1.fst \
			obspack-daily-assim-2.fst \
			tccon-daily.fst \
			obspack-hourly-assim-0.fst \
			obspack-daily-assim-0.fst \
		--resolution hourly hourly hourly hourly daily daily daily daily \
			hourly daily \
		--runs 1_transport/intermediates/runs-r10-r15-rNZ \
		--matched-runs 2_matching/intermediates/runs-r10-r15-rNZ \
		--output-base 3_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-hourly \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1 \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-2 \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-hourly \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1 \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-2 \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-daily \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-0 \
			3_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-0

$(SENSITIVITIES_BASE_PART_1) &: \
	3_inversion/src/sensitivities.R
	Rscript $< \
		--input \
			oco2-hourly.fst \
			obspack-hourly-assim-1.fst \
			obspack-hourly-assim-2.fst \
			tccon-hourly.fst \
			oco2-daily.fst \
			obspack-daily-assim-1.fst \
			obspack-daily-assim-2.fst \
			tccon-daily.fst \
			obspack-hourly-assim-0.fst \
			obspack-daily-assim-0.fst \
		--resolution hourly hourly hourly hourly daily daily daily daily \
			hourly daily \
		--runs 1_transport/intermediates/runs \
		--matched-runs 2_matching/intermediates/runs \
		--output-base 3_inversion/intermediates/sensitivities-base-oco2-hourly \
			3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1 \
			3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-2 \
			3_inversion/intermediates/sensitivities-base-tccon-hourly \
			3_inversion/intermediates/sensitivities-base-oco2-daily \
			3_inversion/intermediates/sensitivities-base-obspack-daily-assim-1 \
			3_inversion/intermediates/sensitivities-base-obspack-daily-assim-2 \
			3_inversion/intermediates/sensitivities-base-tccon-daily \
			3_inversion/intermediates/sensitivities-base-obspack-hourly-assim-0 \
			3_inversion/intermediates/sensitivities-base-obspack-daily-assim-0

## Preliminaries

$(CONSTRAINTS): \
	3_inversion/src/constraints.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(PRIOR): \
	3_inversion/src/prior.R \
	$(PERTURBATIONS) \
	$(BASIS_VECTORS)
	# NOTE(mgnb): the GpGp package fails if the number of threads of greater
	# than one
	OMP_NUM_THREADS=1 Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(BASIS_VECTORS): \
	3_inversion/src/basis-vectors.R \
	$(PERTURBATIONS)
	Rscript $< \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(OBSERVATIONS): \
	3_inversion/src/observations.R \
	$(OCO2_OBSERVATIONS)
	Rscript $< \
		--oco2-observations $(OCO2_OBSERVATIONS) \
		--obspack-directory $(OBSPACK_DIRECTORY) \
		--tccon-sounding-directory $(TCCON_SOUNDING_DIRECTORY) \
		--start-date $(INVERSION_START_DATE) \
		--end-date $(INVERSION_END_DATE) \
		--output $@

$(PERTURBATIONS): \
	3_inversion/src/perturbations.R
	Rscript $< \
		--runs 1_transport/intermediates/runs 1_transport/intermediates/runs-r10-r15-rNZ \
		--matched-runs 2_matching/intermediates/runs 2_matching/intermediates/runs-r10-r15-rNZ \
		--output $@

$(CONTROL_EMISSIONS): \
	3_inversion/src/control-emissions.R
	Rscript $< \
		--matched-runs 2_matching/intermediates/runs \
		--output $@
