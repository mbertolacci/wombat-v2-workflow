build_base_run <- function(
  ext_data,
  fossil_path,
  fires_path,
  ocean_climatology_path,
  ocean_residual_path,
  sib4_assim_climatology_path,
  sib4_assim_residual_path,
  sib4_resp_tot_climatology_path,
  sib4_resp_tot_residual_path
) {
  main_settings <- geoschem_main_settings(
    simulation_menu = geoschem_simulation_menu(
      start = ymd('2014-09-01'),
      end = ymd('2021-04-01'),
      root_data_directory = ext_data,
    ),
    timestep_menu = geoschem_timestep_menu(
      transport_convection_timestep = 600,
      chemistry_emissions_timestep = 1200
    ),
    advected_species_menu = geoschem_advected_species_menu(
      type_of_simulation = 12,
      species = 'CO2'
    ),
    meteorology = 'merra2'
  )

  hemco_diagn <- hemco_diagnostics(
    hemco_diagnostic(
      name = 'EmisCO2_Total',
      species = 'CO2',
      extension_number = -1,
      category = -1,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_total_emissions'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_FossilFuel',
      species = 'CO2',
      extension_number = 0,
      category = 1,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_anthropogenic_emissions'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Ocean_Climatology',
      species = 'CO2',
      extension_number = 0,
      category = 2,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_ocean_emissions_climatology'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Ocean_Residual',
      species = 'CO2',
      extension_number = 0,
      category = 3,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_ocean_emissions_residual'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Biofuel',
      species = 'CO2',
      extension_number = 0,
      category = 4,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_biofuel_emissions'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Assim_Climatology',
      species = 'CO2',
      extension_number = 0,
      category = 5,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_assim_climatology'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Assim_Residual',
      species = 'CO2',
      extension_number = 0,
      category = 6,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_assim_residual'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Resp_Tot_Climatology',
      species = 'CO2',
      extension_number = 0,
      category = 7,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_resp_tot_climatology'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Resp_Tot_Residual',
      species = 'CO2',
      extension_number = 0,
      category = 8,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_resp_tot_residual'
    ),
    hemco_diagnostic(
      name = 'EmisCO2_Biomass',
      species = 'CO2',
      extension_number = 0,
      category = 9,
      hierarchy = -1,
      dimensions = 2,
      output_unit = 'kg/m2/s',
      long_name = 'CO2_biomass_burning_emissions'
    )
  )

  history <- geoschem_history(
    file_exp_id = './output/GEOSChem',
    collections = list(
      geoschem_history_collection(
        name = 'Restart',
        frequency = '00000100 000000',
        duration = '00000100 000000',
        mode = 'instantaneous',
        fields = c(
          'SpeciesRst_?ALL?', 'Met_DELPDRY', 'Met_PS1WET', 'Met_PS1DRY',
          'Met_SPHU1', 'Met_TMPU1'
        )
      ),
      geoschem_history_collection(
        name = 'SpeciesConcHourly',
        frequency = '00000000 010000',
        duration = '00000001 000000',
        mode = 'time-averaged',
        fields = 'SpeciesConc_?ADV?',
      ),
      geoschem_history_collection(
        name = 'LevelEdgeDiagsHourly',
        frequency = '00000000 010000',
        duration = '00000001 000000',
        mode = 'time-averaged',
        fields = c('Met_PEDGE', 'Met_PEDGEDRY')
      ),
      geoschem_history_collection(
        name = 'StateMetHourly',
        frequency = '00000000 010000',
        duration = '00000001 000000',
        mode = 'time-averaged',
        fields = c('Met_AVGW', 'Met_BXHEIGHT')
      ),
      geoschem_history_collection(
        name = 'SpeciesConcDaily',
        frequency = '00000001 000000',
        duration = '00000100 000000',
        mode = 'time-averaged',
        fields = 'SpeciesConc_?ADV?',
      ),
      geoschem_history_collection(
        name = 'LevelEdgeDiagsDaily',
        frequency = '00000001 000000',
        duration = '00000100 000000',
        mode = 'time-averaged',
        fields = c('Met_PEDGE', 'Met_PEDGEDRY')
      ),
      geoschem_history_collection(
        name = 'StateMetDaily',
        frequency = '00000001 000000',
        duration = '00000100 000000',
        mode = 'time-averaged',
        fields = c('Met_AVGW', 'Met_BXHEIGHT')
      )
    )
  )

  scale_factors <- hemco_scale_factors(
    # qfed2_diurnal = hemco_scale_factor(
    #   id = 4,
    #   name = 'QFED2_TOD',
    #   source_file = '0.1392/0.1392/0.1368/0.1368/0.1368/0.1368/0.1368/0.1368/0.1368/0.48/0.96/1.68/2.4/3.12/3.84/4.08/2.88/1.68/0.96/0.1368/0.1368/0.1368/0.1368/0.1368',
    #   source_variable = NULL,
    #   source_time = NULL,
    #   cre = NULL,
    #   source_dimension = 'xy',
    #   source_unit = '1',
    #   operator = hemco_scale_factor_operator('multiplication')
    # ),
    co2_negative = hemco_scale_factor(
      id = 5,
      name = 'CO2_NEGATIVE',
      source_file = '-1.0',
      source_variable = NULL,
      source_time = '2000/1/1/0',
      cre = 'C',
      source_dimension = 'xy',
      source_unit = '1',
      operator = hemco_scale_factor_operator('multiplication')
    )
  )

  base_emission_fields <- list(
    hemco_base_emission_field(
      extension_number = 0,
      name = 'FOSSILCO2_MIPV10',
      source_file = file.path(
        fossil_path,
        '$YYYY/$MM/fossil_fuel_1x1_$YYYY$MM$DD.nc'
      ),
      source_variable = 'CO2',
      source_time = '2014-2021/1-12/1-31/0-23',
      cre = 'C',
      source_dimension = 'xy',
      source_unit = 'kg/m2/s',
      species = 'CO2',
      scale_factors = NULL,
      category = 1,
      hierarchy = 2
    ),
    # hemco_base_emission_field(
    #   extension_number = 0,
    #   name = 'QFED',
    #   source_file = '$ROOT/QFED/v2018-07/$YYYY/$MM/qfed2.emis_co2.006.$YYYY$MM$DD.nc4',
    #   source_variable = 'biomass',
    #   source_time = '2000-2017/1-12/1-31/0',
    #   cre = 'C',
    #   source_dimension = 'xy',
    #   source_unit = 'kg/m2/s',
    #   species = 'CO2',
    #   scale_factors = scale_factors$qfed2_diurnal,
    #   category = 9,
    #   hierarchy = 1
    # ),
    hemco_base_emission_field(
      extension_number = 0,
      name = 'GFED4',
      source_file = file.path(
        fires_path,
        'gfed4_1deg_$YYYY$MM.nc'
      ),
      source_variable = 'emis_CO2',
      source_time = '2014-2021/1-12/1-31/0-23/+90minutes',
      cre = 'RF',
      source_dimension = 'xy',
      source_unit = 'kg/m2/s',
      species = 'CO2',
      scale_factors = NULL,
      category = 9,
      hierarchy = 1
    ),
    hemco_base_emission_field(
      extension_number = 0,
      name = 'BIOFUEL_CO2',
      source_file = '$ROOT/CO2/v2015-04/BIOFUEL/biofuel_CO2.geos.1x1-1995.nc',
      source_variable = 'CO2',
      source_time = '1995/1/1/0',
      cre = 'C',
      source_dimension = 'xy',
      source_unit = 'kg/m2/s',
      species = 'CO2',
      scale_factors = NULL,
      category = 4,
      hierarchy = 1
    )
  )

  for (part in c(
    'intercept',
    'trend',
    'cos12_1_intercept',
    'sin12_1_intercept',
    'cos12_2_intercept',
    'sin12_2_intercept',
    'cos12_1_trend',
    'sin12_1_trend',
    'cos12_2_trend',
    'sin12_2_trend',
    'residual'
  )) {
    base_emission_fields <- c(base_emission_fields, list(
      hemco_base_emission_field(
        extension_number = 0,
        name = sprintf('OCEAN_LSCHULZ_%s', toupper(part)),
        source_file = if (part == 'residual') {
          ocean_residual_path
        } else {
          ocean_climatology_path
        },
        source_variable = part,
        source_time = if (part == 'residual') {
          '2014-2019/1-12/1/0'
        } else {
          '2014-2021/1-12/1/0'
        },
        cre = if (part == 'residual') {
          'C'
        } else {
          'RF'
        },
        source_dimension = 'xy',
        source_unit = 'kg/m2/s',
        species = 'CO2',
        scale_factors = NULL,
        category = if (part == 'residual') {
          3
        } else {
          2
        },
        hierarchy = 3
      )
    ))
  }

  for (field in c('assim', 'resp_tot')) {
    for (part in c(
      'intercept',
      'trend',
      'cos12_1_intercept',
      'sin12_1_intercept',
      'cos12_2_intercept',
      'sin12_2_intercept',
      'cos12_3_intercept',
      'sin12_3_intercept',
      'cos12_1_trend',
      'sin12_1_trend',
      'cos12_2_trend',
      'sin12_2_trend',
      'cos12_3_trend',
      'sin12_3_trend',
      'residual'
    )) {
      base_emission_fields <- c(base_emission_fields, list(
        hemco_base_emission_field(
          extension_number = 0,
          name = sprintf('SIB4_%s_%s', toupper(field), toupper(part)),
          source_file = if (field == 'assim') {
            if (part == 'residual') {
              sib4_assim_residual_path
            } else {
              sib4_assim_climatology_path
            }
          } else {
            if (part == 'residual') {
              sib4_resp_tot_residual_path
            } else {
              sib4_resp_tot_climatology_path
            }
          },
          source_variable = part,
          source_time = if (part == 'residual') {
            '2014-2020/1-12/1-31/0-23'
          } else {
            '2014-2021/1-12/1-31/0-23'
          },
          cre = if (part == 'residual') {
            'C'
          } else {
            'RF'
          },
          source_dimension = 'xy',
          source_unit = 'kg/m2/s',
          species = 'CO2',
          scale_factors = if (field == 'assim') {
            scale_factors$co2_negative
          } else NULL,
          category = c(
            'assim' = if (part == 'residual') 6 else 5,
            'resp_tot' = if (part == 'residual') 8 else 7
          )[field],
          hierarchy = 1
        )
      ))
    }
  }

  base_emissions <- do.call(hemco_base_emissions, base_emission_fields)

  hemco_set <- hemco_settings(
    root_directory = file.path(ext_data, 'HEMCO'),
    meteorology_directory = file.path(ext_data, 'GEOS_2x2.5/MERRA2'),
    log_file = 'output/HEMCO.log',
    diagnostics_file = 'HEMCO_Diagn.rc',
    diagnostics_prefix = 'output/HEMCO_diagnostics',
    diagnostics_frequency = '00000000 010000',
    diagnostics_time_stamp = 'Mid',
    wildcard = '*',
    separator = '/',
    unit_tolerance = 1,
    negative_values = 2,
    only_unitless_scale_factors = FALSE,
    verbose = 0,
    warnings = 1
  )

  extension_switches <- hemco_extension_switches(
    hemco_extension_switch(
      0,
      'Base',
      TRUE,
      '*',
      c('GC_RESTART' = TRUE, 'MODIS_XLAI' = TRUE, 'Yuan_XLAI' = FALSE)
    )
  )

  configuration <- hemco_configuration(
    hemco_set,
    extension_switches,
    base_emissions = base_emissions,
    scale_factors = scale_factors
  )

  geoschem_run_configuration(
    main_settings,
    history,
    configuration,
    hemco_diagn,
    extra_directories = 'output'
  )
}
