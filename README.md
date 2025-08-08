This readme file was generated on 2025-02-13 by Tara Jacobsen

General Information
------------------------------------------------------------------------------

Title: NASA's Confirmed Exoplanets:

Author: JamesKYChoi on Kaggle https://www.kaggle.com/datasets/jameskychoi/confirmed-exoplanet-latest-update-dataset 

Additional cleaning by: Tara Jacobsen

Raw Data Source: NASA Exoplanet Archive https://exoplanetarchive.ipac.caltech.edu/index.html

Data of access: 2025-01-30

Data Overview
------------------------------------------------------------------------------

As of December 10, 2024, when this dataset was last updated, a total of 5,788 exoplanets have been confirmed. This dataset contains detailed information about these exoplanets organized into over 200 columns. Data is regularly updated, but is not a live dataset. 

Data is originally sourced from the NASA Exoplanet Archive. "The NASA Exoplanet Archive collects and maintains many sets of parameters for planets, planetary orbits, and host systems as they are published in the refereed literature." (https://exoplanetarchive.ipac.caltech.edu/docs/pscp_about.html)


Methodological Information
------------------------------------------------------------------------------

Data was originally collected by NASA and stored by the NASA Exoplanet Archive. It was extracted from the NASA CSV file into a SQL Database where a view was created to extract over 200 columns.

The data was later cleaned again with Python. All categorical data was transformed into lowercase characters. All columns were ensured to be lowercase with '_' for spaces. Data was checked for duplicated. Few choice columns with large numbers of NA values were removed, including occultation depth-related columns. Discovery publication date was changed to a datetime. All strings within the dataset were made lowercase. 

"Data validation (DV) summaries are provided by the Kepler Mission that show the results of the data validation tests conducted for the object." (https://exoplanetarchive.ipac.caltech.edu/docs/ICEexohelp.html)

Data is public domain. 

Units are included in the column description, where applicable. 

Column Descriptions
------------------------------------------------------------------------------

| Column | The column indicates... |
|----------|----------|
| planet_name | Planet name most commonly used in the literature | 
| host_star_name  | Stellar name most commonly used in the literature | 
| planet_letter | Letter assigned to the planetary component of a planetary system | 
| hd_id | Name of the star as given by the Henry Draper Catalog |
| hip_id | Name of the star as given by the Hipparcos Catalog |
| tic_id | Name of the star as given by the TESS Input Catalog|
| gaia_id | Name of the star as given by the Gaia Catalog|
|num_stars_in_system | Number of stars in the planetary system|
| num_planets_in_system|Number of confirmed planets in the planetary system |
| num_moons_in_system|Number of moons in the planetary system |
| is_circumbinary| Flag indicating whether the planet orbits a binary system (1=yes, 0=no)|
| discovery_method| Method by which the planet was first identified|
| discovery_year| Year the planet was discovered|
|discovery_publication_date  |Publication Date of the planet discovery referee publication	(datetime) |
|discovery_locale |Location of observation of planet discovery (ground or space) |
|discovery_facility |Name of facility of planet discovery observations |
| discovery_telescope |Name of telescope of planet discovery observations |
|discovery_instrument | Name of instrument of planet discovery observations|
| detected_by_radial_velocity|Flag indicating if the planet host star exhibits radial velocity variations due to the planet (1=yes, 0=no)	 |
|detected_by_pulsar_timing |Boolean flag indicating if the planet host star exhibits pulsar timing variations due to the planet (1=yes, 0=no)	 |
| detected_by_pulsation_timing_variations|Boolean flag indicating if the planet host star exhibits pulsation timing variations due to the planet (1=yes, 0=no)|
| detected_by_transit| Flag indicating if the planet transits its host star (1=yes, 0=no)	|
| detected_by_astrometry| Flag indicating if the planet host star exhibits astrometrical variations due to the planet (1=yes, 0=no)	|
| detected_by_orbital_brightness_modulation|Flag indicating whether the planet exhibits orbital modulations on the phase curve (1=yes, 0=no)	 |
|detected_by_microlensing | Boolean flag indicating if the planetary system acted as a lens during an observed microlensing event (1=yes, 0=no)	|
| detected_by_eclipse_timing_variations| Flag indicating whether a circumbinary planet that orbits an eclipsing binary induces eclipse timing variations (ETVs) in the binary pair (1=yes, 0=no)	|
| detected_by_imaging| Flag indicating if the planet has been observed via imaging techniques (1=yes, 0=no)	|
|detected_by_disk_kinematics | Boolean flag indicating if the presence of the planet was inferred due to its kinematic influence on the protoplanetary disk of its host star (1=yes, 0=no)	|
| is_controversial|Boolean flag indicating if the exoplanet discovery or characteristics are controversial (1=yes, 0=no) |
|orbital_period_days |Time the planet takes to make a complete orbit around the host star or system (days)	 |
| orbital_period_upper_unc_days| Upper uncertainty in the time (in days) the planet takes to complete an orbit around its host star or system |
| orbital_period_lower_unc_days| Lower uncertainty in the time (in days) the planet takes to complete an orbit around its host star or system |
|orbital_period_limit_flag | Flag indicating whether the orbital period is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| orbit_semi_major_axis_au| The longest radius of an elliptic orbit, or, for exoplanets detected via gravitational microlensing or direct imaging, the projected separation in the plane of the sky	|
|orbit_semi_major_axis_upper_unc_au | Upper uncertainty in the orbit semi-major axis (au) |
| orbit_semi_major_axis_lower_unc_au| Lower uncertainty in the orbit semi-major axis (au)|
| orbit_semi_major_axis_limit_flag| Flag indicating whether the orbit semi-major axis is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| angular_separation_mas| Angular separation on the sky between the star and planet in milliarcseconds (mas). |
| angular_separation_limit_flag|Flag indicating whether the angular seperation is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| planet_radius_earth_radius|Length of a line segment from the center of the planet to its surface, measured in units of radius of the Earth	 |
|planet_radius_upper_unc_earth_radius | Upper uncertainty in the planet radius, measured in units of radius of Earth |
|planet_radius_lower_unc_earth_radius | Lower uncertainty in the planet radius, measured in units of radius of Earth|
| planet_radius_earth_limit_flag|Flag indicating whether the planet radius is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|planet_radius_jupiter_radius |Length of a line segment from the center of the planet to its surface, measured in units of radius of Jupiter |
|planet_radius_upper_unc_jupiter_radius | Upper uncertainty in the planet radius, measured in units of radius of Jupiter|
|planet_radius_lower_unc_jupiter_radius |Lower uncertainty in the planet radius, measured in units of radius of Jupiter |
| planet_radius_jupiter_limit_flag|Flag indicating whether the planet radius is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|planet_mass_earth_mass | Amount of matter contained in the planet, measured in units of masses of the Earth|
| planet_mass_upper_unc_earth_mass|Upper uncertainty in the planet mass, measured in units of masses of the Earth |
| planet_mass_lower_unc_earth_mass|Lower uncertainty in the planet mass, measured in units of masses of the Earth |
|planet_mass_earth_limit_flag |Flag indicating whether the planet mass is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| planet_mass_jupiter_mass| Amount of matter contained in the planet, measured in units of masses of Jupiter	|
|planet_mass_upper_unc_jupiter_mass |Upper uncertainty in the planet mass, measured in units of masses of the Jupiter |
|planet_mass_lower_unc_jupiter_mass |Lower uncertainty in the planet mass, measured in units of masses of the Jupiter |
| planet_mass_jupiter_limit_flag|Flag indicating whether the planet mass is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| planet_mass_provenance|Provenance of the measurement of the best mass. Options are: Mass, M*sin(i)/sin(i), and M*sini	 |
| planet_density_gcm3| Amount of mass per unit of volume of the planet (g/cm**2)|
|planet_density_upper_unc_gcm3 |Upper uncertainty in the planet density (g/cm**2)  |
| planet_density_lower_unc_gcm3|Lower uncertainty in the planet density (g/cm**2)  |
|planet_density_limit_flag |Flag indicating whether the planet density is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|eccentricity | Amount by which the orbit of the planet deviates from a perfect circle	|
|eccentricity_upper_unc | Upper uncertainty in the eccentricity|
|eccentricity_lower_unc |Lower uncertainty in the eccentricity |
| eccentricity_limit_flag| Flag indicating whether the eccentricity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| insolation_flux_earth_flux| Insolation flux is another way to give the equilibrium temperature. It's given in units relative to those measured for the Earth from the Sun.	|
| insolation_flux_upper_unc_earth_flux| Upper uncertainty in the insolition flux |
|insolation_flux_lower_unc_earth_flux |Lower uncertainty in the insolition flux  |
|insolation_flux_limit_flag |Flag indicating whether the insolation flux is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|equilibrium_temperature_k |The equilibrium temperature of the planet as modeled by a black body heated only by its host star, or for directly imaged planets, the effective temperature of the planet required to match the measured luminosity if the planet were a black body (K) |
|equilibrium_temperature_upper_unc_k | Upper uncertainty in the equilibrium temperature|
| equilibrium_temperature_lower_unc_k| Lower uncertainty in the equilibrium temperature|
| equilibrium_temperature_limit_flag|Flag indicating whether the equilibrium temperature flux is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|inclination_deg | Angle of the plane of the orbit relative to the plane perpendicular to the line-of-sight from Earth to the object (degrees)|
|inclination_upper_unc_deg | Upper uncertainty in the inclination degree|
|inclination_lower_unc_deg | Upper uncertainty in the inclination degree|
| inclination_limit_flag|Flag indicating whether the inclination degree is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|transit_midpoint_days | The time given by the average of the time the planet begins to cross the stellar limb and the time the planet finishes crossing the stellar limb.	|
|transit_midpoint_upper_unc_days |Upper uncertainty in the transit midpoint |
|transit_midpoint_lower_unc_days | Lower uncertainty in the transit midpoint |
|transit_midpoint_limit_flag |Flag indicating whether the transit midpoint is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| transit_midpoint_time_system| Time system used to record the transit midpoint|
| data_show_transit_timing_variations| Flag indicating if the planet orbit exhibits transit timing variations from another planet in the system (1=yes, 0=no).|
| impact_parameter| The sky-projected distance between the center of the stellar disc and the center of the planet disc at conjunction, normalized by the stellar radius|
|impact_parameter_upper_unc |Upper uncertainty in the impact parameter |
| impact_parameter_lower_unc|Lower uncertainty in the impact parameter |
|impact_parameter_limit_flag |Flag indicating whether the impact parameter is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|transit_depth_percent |The size of the relative flux decrement caused by the orbiting body transiting in front of the star (%) |
|transit_depth_upper_unc_percent | Upper uncertainty in the transit depth %|
|transit_depth_lower_unc_percent |Lower uncertainty in the transit depth % |
|transit_depth_limit_flag |Flag indicating whether the transit depth is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| transit_duration_hours|The length of time from the moment the planet begins to cross the stellar limb to the moment the planet finishes crossing the stellar limb (hours) |
|transit_duration_upper_unc_hours | Upper uncertainty in the transit duration|
|transit_duration_lower_unc_hours | Lower uncertainty in the transit duration|
|transit_duration_limit_flag | Flag indicating whether the transit duration is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|ratio_semi_major_axis_to_stellar_radius|The distance between the planet and the star at mid-transit divided by the stellar radius. For the case of zero orbital eccentricity, the distance at mid-transit is the semi-major axis of the planetary orbit. |
|ratio_semi_major_axis_to_stellar_radius_upper_unc |Upper uncertainty in the semi major axis ratio |
|ratio_semi_major_axis_to_stellar_radius_lower_unc |Lower uncertainty in the semi major axis ratio  |
|ratio_semi_major_axis_to_stellar_radius_limit_flag |Flag indicating whether the semi major axis ratio is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|ratio_planet_to_stellar_radius| The planet radius divided by the stellar radius|
|ratio_planet_to_stellar_radius_upper_unc | Upper uncertainty in the ratio of the planet to stellar radius |
|ratio_planet_to_stellar_radius_lower_unc | Lower uncertainty in the ratio of the planet to stellar radius|
| ratio_planet_to_stellar_radius_limit_flag| Flag indicating whether the ratio of the planet to stellar radius is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|epoch_of_periastron_days | The time of the planet's periastron passage	(days)|
|epoch_of_periastron_upper_unc_days | Upper uncertainty in the time of the planet's periastron passage (days)|
| epoch_of_periastron_lower_unc_days|Lower uncertainty in the time of the planet's periastron passage (days) |
|epoch_of_periastron_limit_flag | Flag indicating whether the time of the planet's periastron passage is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|epoch_of_periastron_time_system | |
| argument_of_periastron_deg| The angular separation between the orbit's ascending node and periastron. Note: there are a varying conventions in the exoplanet literature regarding argument of periastron (or periapsis). For example, some publications refer the planet's orbit, others to the host star's reflex orbit, which differs by 180 deg. The values in the Exoplanet Archive are not corrected to a standardized system, but are as-reported for each publication.	|
|argument_of_periastron_upper_unc_deg | Upper uncertainty in the angular seperation between the orbit's ascending node and periastron|
|argument_of_periastron_lower_unc_deg |Lower uncertainty in the angular seperation between the orbit's ascending node and periastron |
|argument_of_periastron_limit_flag |Flag indicating whether the angular seperation between the orbit's ascending node and periastron is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|radial_velocity_amplitude_ms |Half the peak-to-peak amplitude of variability in the stellar radial velocity (m/s) |
| radial_velocity_amplitude_upper_unc_ms| Upper uncertainty in the radial velocity amplitude |
|radial_velocity_amplitude_lower_unc_ms |Lower uncertainty in the radial velocity amplitude |
|radial_velocity_amplitude_limit_flag |Flag indicating whether the radial velocity amplitude is based on a limit (-1=yes, 0=no) |
|projected_obliquity_deg |The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet, projected into the plane of the sky. Depending on the choice of coordinate system, projected obliquity is represented in the literature as either lambda (λ) or beta (β), where λ is defined as the negative of β (i.e., λ = -β). Since λ is reported more often than β, all values of β have been converted to λ.|
|projected_obliquity_upper_unc_deg | Upper uncertainty in the projected obliquity degree |
| projected_obliquity_lower_unc_deg|Lower uncertainty in the projected obliquity degree |
| projected_obliquity_limit_flag|Flag indicating whether the projected obliquity degree is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|true_obliquity_deg | The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet (degrees)|
|true_obliquity_upper_unc_deg |Upper uncertainty in the true obliquity degree |
|true_obliquity_lower_unc_deg |Lower uncertainty in the true obliquity degree |
|true_obliquity_limit_flag |Flag indicating whether the true obliquity degree is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|spectral_type | Classification of the star based on their spectral characteristics following the Morgan-Keenan system	|
|stellar_effective_temp_k | Temperature of the star as modeled by a black body emitting the same total amount of electromagnetic radiation	 (K)|
| stellar_effective_temp_upper_unc_k| Upper uncertainty in the effective stellar temperature (K)|
|stellar_effective_temp_lower_unc_k | Lower uncertainty in the effective stellar temperature (K)|
| stellar_effective_temp_limit_flag|Flag indicating whether the effective stellar temperature is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|stellar_radius_solar_radius | Length of a line segment from the center of the star to its surface, measured in units of radius of the Sun	|
|stellar_radius_upper_unc_solar_radius | Upper uncertainty in the stellar radius, measured in units of radius of the Sun |
|stellar_radius_lower_unc_solar_radius | Lower uncertainty in the stellar radius, measured in units of radius of the Sun|
|stellar_radius_limit_flag | Flag indicating whether the stellar radius is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|stellar_mass_solar_mass |Amount of matter contained in the star, measured in units of masses of the Sun |
| stellar_mass_upper_unc_solar_mass| Upper uncertainty in the stellar mass|
|stellar_mass_lower_unc_solar_mass | Lower uncertainty in the stellar mass |
|stellar_mass_limit_flag |Flag indicating whether the stellar mass is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
|stellar_metallicity_dex |Measurement of the metal content of the photosphere of the star as compared to the hydrogen content (dex) |
| stellar_metallicity_upper_unc_dex| Upper uncertainty in the stellar metallicity|
| stellar_metallicity_lower_unc_dex|Lower uncertainty in the stellar metallicity |
|stellar_metallicity_limit_flag | Flag indicating whether the stellar metallicity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|stellar_metallicity_ratio |Ratio for the Metallicity Value ([Fe/H] denotes iron abundance, [M/H] refers to a general metal content)	 |
| stellar_luminosity_log_solar|Amount of energy emitted by a star per unit time, measured in units of solar luminosities	 |
|stellar_luminosity_upper_unc_log_solar |Upper uncertainty in the stellar luminosity |
|stellar_luminosity_lower_unc_log_solar | Lower uncertainty in the stellar luminosity|
| stellar_luminosity_limit_flag| Flag indicating whether the stellar luminosity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| stellar_surface_gravity_log10_cms2|Gravitational acceleration experienced at the stellar surface (log10(cm/s**2) |
|stellar_surface_gravity_upper_unc_log10_cms2 |Upper uncertainty in the stellar surface gravity |
|stellar_surface_gravity_lower_unc_log10_cms2 |Lower uncertainty in the stellar surface gravity |
|stellar_surface_gravity_limit_flag | Flag indicating whether the stellar surface gravity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
|stellar_age_gyr | The age of the host star	(gigayear)|
|stellar_age_upper_unc_gyr| Upper uncertainty in the stellar age|
|stellar_age_lower_unc_gyr |Lower uncertainty in the stellar age |
|stellar_age_limit_flag | Flag indicating whether the stellar age is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| stellar_density_gcm3 |Amount of mass per unit of volume of the star (g/cm**3) |
| stellar_density_upper_unc_gcm3 |Upper uncertainty in the stellar density |
| stellar_density_lower_unc_gcm3 |Lower uncertainty in the stellar density |
| stellar_density_limit_flag | Flag indicating whether the stellar density is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)|
| stellar_rotational_velocity_kms |Rotational velocity at the equator of the star multiplied by the sine of the inclination (km/s) |
| stellar_rotational_velocity_upper_unc_kms |Upper uncertainty in the stellar rotational velocity |
| stellar_rotational_velocity_lower_unc_kms |Lower uncertainty in the stellar rotational velocity  |
| stellar_rotational_velocity_limit_flag |Flag indicating whether the stellar rotational velocity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| stellar_rotational_period_days |The time required for the planet host star to complete one rotation, assuming it is a solid body	(days) |
| stellar_rotational_period_upper_unc_days | Upper uncertainty in the stellar rotational period |
| stellar_rotational_period_lower_unc_days |Lower uncertainty in the stellar rotational period |
| stellar_rotational_period_limit_flag |Flag indicating whether the stellar rotational period is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit)  |
| systemic_radial_velocity_kms | Velocity of the star in the direction of the line of sight	(km/s)|
| systemic_radial_velocity_upper_unc_kms | Upper uncertainty in the systemic radial velocity|
| systemic_radial_velocity_lower_unc_kms | Lower uncertainty in the systemic radial velocity|
| systemic_radial_velocity_limit_flag |Flag indicating whether the systemic radial velocity is based on a limit (-1 = below detection limit, 0 = no limit/constraint, 1 = above detection limit) |
| ra_str |Right Ascension of the planetary system in sexagesimal format	 |
| ra_deg |Right Ascension of the planetary system in decimal degrees	 |
| dec_str |Declination of the planetary system in sexagesimal notation	 |
| dec_deg |Declination of the planetary system in decimal degrees	 |
| galactic_latitude_deg | Galactic latitude of the planetary system in units of decimal degrees	|
| galactic_longitude_deg |Galactic longitude of the planetary system in units of decimal degrees	 |
| ecliptic_latitude_deg |Ecliptic latitude of the planetary system in units of decimal degrees	 |
| ecliptic_longitude_deg |Ecliptic longitude of the planetary system in units of decimal degrees	 |
| total_proper_motion_mas_yr |Angular change in position over time as seen from the center of mass of the Solar System (mas/yr)|
| total_proper_motion_upper_unc_mas_yr |Upper uncertainty in the total proper motion |
| total_proper_motion_lower_unc_mas_yr |Lower uncertainty in the total proper motion |
| proper_motion_ra_mas_yr | Angular change in right ascension over time as seen from the center of mass of the Solar System	(mas/yr)|
| proper_motion_ra_upper_unc_mas_yr | Upper uncertainty in the proper motion (RA)|
| proper_motion_ra_lower_unc_mas_yr | Lower uncertainty in the proper motion (RA)|
| proper_motion_dec_mas_yr |Angular change in declination over time as seen from the center of mass of the Solar System (mas/yr) |
| proper_motion_dec_upper_unc_mas_yr |Upper uncertainty in the proper motion (Dec) |
| proper_motion_dec_lower_unc_mas_yr |Lower uncertainty in the proper motion (Dec) |
| distance_pc | Distance to the planetary system in units of parsecs (pc)|
| distance_upper_unc_pc | Upper uncertainty in the distance|
| distance_lower_unc_pc |Lower uncertainty in the distance |
| parallax_mas |Difference in the angular position of a star as measured at two opposite positions within the Earth's orbit (mas) |
| parallax_upper_unc_mas | Upper uncertainty in the parallax|
| parallax_lower_unc_mas |Lower uncertainty in the parallax |
| b_magnitude_johnson |Brightness of the host star as measured using the B (Johnson) band in units of magnitudes |
| b_magnitude_johnson_upper_unc | Upper uncertainty in the B (Johnson) Magnitude|
| b_magnitude_johnson_lower_unc | Lower uncertainty in the B (Johnson) Magnitude|
| v_magnitude_johnson |Brightness of the host star as measured using the V (Johnson) band in units of magnitudes	 |
| v_magnitude_johnson_upper_unc |Upper uncertainty in the V (Johnson) Magnitude |
| v_magnitude_johnson_lower_unc | Lower uncertainty in the V (Johnson) Magnitude|
| j_magnitude_2mass |Brightness of the host star as measured using the J (2MASS) band in units of magnitudes |
| j_magnitude_2mass_upper_unc |Upper uncertainty in the J (2MASS) band magnitude |
| j_magnitude_2mass_lower_unc |Lower uncertainty in the J (2MASS) band magnitude |
| h_magnitude_2mass |Brightness of the host star as measured using the H (2MASS) band in units of magnitudes	 |
| h_magnitude_2mass_upper_unc |Upper uncertainty in the K (2MASS) band magnitude |
| h_magnitude_2mass_lower_unc |Lower uncertainty in the K (2MASS) band magnitude |

Column descriptions courtesy of NASA: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html 
