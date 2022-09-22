# CVI - Climate Vulnerability Index calculations, analyses, and visualizations

Data and analyses code for

Lewis PGT, Chiu WA, Nasser E, Proville J, Barone A, Danforth C, Kim B, Prozzi J, Craft E. 2022.  Characterizing Vulnerabilities to Climate Change Across the United States (in review)

Data can be visualized on an ArcGIS dashboard  available at https://arcg.is/1bvKjK0.

To recreate all analyses and figures, run the scripts in the following order. NOTE: *Script (1) cannot be run without access to the Dropbox that has all the original data files.* However, all other scripts can be run.
 
(1) Source get_datasets.R to load all data from individual files into one master file. Creates the following files:
- CVI_master.csv - all details as to each indicator, including source file information and notes
- Checkoutput.txt - output capture during data loading, used for diagnostics
- CVI_indicators_current.csv - most important meta-data for each indicator
- CVI_data_current.csv - indicator values (columns) for each census tract (rows)

(2) Source get_cvi_toxpi.R to convert to ToxPi indices.
This script conducts data post-processing (e.g., imputation using median, applying transformations for adverse direction, calculating percentiles), and then calculates ToxPI index values. Creates the following files:
- Diagnostics/CheckDist.pdf - shows distribution of each indicator
- Diagnostics/CVI-corr.pdf - rank correlations across indicators
- CVI-pct/CVI_data_pct.csv - indicator data transformed to percentiles
- CVI-pct/CVI-pct.pdf - Boxplots showing distribution of all indicators
- CVI-pct/CVI-pct-allinone.csv - ToxPI treating all indicators as equally weighted (not used)
- CVI-pct/ToxPi-pct-allinone.pdf - Example ToxPI profiles (top 10 and bottom 10), overall distribution of scores, and correlations across category domains
- CVI-pct/CVI-pct-cat-XXXX.csv - Category domain-specific ToxPI files for category XXXX
- CVI-pct/ToxPi-pct-subcat.pdf - Example category-specific ToxPI profiles (top 10 and bottom 10), overall distribution of scores, and correlations across subcategory domains
- CVI-pct/CVI-pct-comb.csv - overall CVI ToxPI
- CVI-pct/ToxPi-pct-subcat-comb.pdf - Example category-specific ToxPI profiles (top 10 and bottom 10), overall distribution of scores, and correlations across category domains
- CVI-pct/CVI-pct-comb-baseline.csv - ToxPI for only Baseline vulnerability categories
- CVI-pct/ToxPi-pct-subcat-comb-baseline.pdf - Example category-specific ToxPI profiles for Baseline vulnerability (top 10 and bottom 10), overall distribution of scores, and correlations across category domains
- CVI-pct/CVI-pct-comb-climate.csv - ToxPI for only Climate change risk categories
- CVI-pct/ToxPi-pct-subcat-comb-climate.pdf - Example category-specific ToxPI profiles for Climate Change risk (top 10 and bottom 10), overall distribution of scores, and correlations across category domains
Note - Files with ".gis.csv" are reformatted versions for use by ToxPI-gis python scripts to create ESRI GIS maps. Files with ".Rdata" are Rdata files with the ToxPI slices, model, and results stored as R objects.

(3) Source get_cvi_county_toxpi.R to make county median ToxPI indices. Creates the following files:
- CVI-county-pct/CVI-county_data_pct.csv - Median indicator percentiles aggregated by county
- CVI-county-pct/CVI-county-pct-ZZZZ.csv - same files as (2) but aggregated by county median
- CVI-county-pct/ToxPi-county-pct-ZZZZ.pdf - same files as (2) but aggregated by county median

(4) Source plot_cvi_county_pct.R to do county maps for QA. Creates the following files:
- CVI-county-pct/CVI-county_data_pct.pdf - Maps of county medians for each indicator
- CVI-county-pct/CVI-county_data_maxpct.csv - Max indicator percentiles aggregated by county
- CVI-county-pct/CVI-county_data_fracNA.csv - Fraction NA of county for each indicator
- CVI-county-pct/CVI-county_data_maxpct.pdf - Maps of county max for each indicator
- CVI-county-pct/CVI-county_data_fracNA.pdf - Maps of county fraction NA for each indicator

(5) Source State_county_summaries.R to make county maps for paper. Creates the following files:

- Scores_summary.csv - summary of scores for each census tract, including top category, used to make figure
- Figures/CVI_maps_8.pdf - county level median maps with overall CVI and individual categories
- Figures/CVI_maps_3.pdf - county level median maps with overall CVI, baseline and climate themes
- Figures/CVI_maps_10.pdf - counmty level median maps with overall CVI, baseline and climate themes, and individual categories
- SuppFigures/CVI_maps_max.pdf - county level max maps with overall CVI and individual categories
- SuppFigures/CVI_maps_relvar.pdf - county level relative variance (county variance/overall variance) maps with overall CVI and individual categories
- Figures/StateSummaryBoxplots.pdf - Summary boxplots by state for overall CVI and individual categories
- Figures/GeoScale.pdf - Bar chart of geographic scale of indicators by category
- Figures/R2.pdf - Bar chart of fraction of variance due to state or county by category
- Figures/Top.Categories.pdf - Bar chart of frequency of top category by category
- Figures/State_County_summary_fig.pdf - Combination of StateSummaryBoxplots, GeoScale, R2, Top.Categories into single figure

(6) Source get_cvi_state_toxpi.R to make state median ToxPi indices

- CVI-state-pct/CVI-state_data_pct.csv - Median indicator percentiles aggregated by state
- CVI-state-pct/CVI-state-pct-ZZZZ.csv - same files as (2) but aggregated by state median
- CVI-state-pct/ToxPi-state-pct-ZZZZ.pdf - same files as (2) but aggregated by state median

(7) Source Merge_shapefile_CVI.R to make tract-level shapefiles for choropleths. Creates shapfiles in the following folders:
- Shapefiles/CVI Data Tracts - shapefile with census tract-level indicators in attribute table
- Shapefiles/CVI Pct Tracts - shapefile with census tract-level indicator percentiles in attribute table
- Shapefiles/CVI ToxPi Tracts - shapefile with census tract-level ToxPi CVI scores in attribute table
- Shapefiles/CVI Baseline ToxPi Tracts - shapefile with census tract-level ToxPi CVI Baseline scores in attribute table
- Shapefiles/CVI Climate ToxPi Tracts - shapefile with census tract-level ToxPi CVI Climate scores in attribute table
- Shapefiles/CVI ToxPi.XXXX - shapefile with census tract-level category XXXX and subcategory scores in attribute table

Uses shapefiles in Data/2010 Tracts for 2010 census tracts

(8) Source State_county_subcat_summaries.R to make category maps. Creates following PDF files:
- SuppFigures/CVI_category-XXXX-map.pdf - county level choropleth maps by category XXXX
- SuppFigures/CVI_category-XXXX-summary.pdf - county level summary boxplots, geographic scale, fraction of variance due to state or county, and correlation among subcategories by category XXXX

(9) Source State_county_parameter_summaries.R to make subcategory maps. Creates following PDF files:
- SuppFigures/CVI_subcategorymap-XXXX.pdf - county level choropleth maps by subcategory XXXX
- SuppFigures/CVI_subcategorymap-XXXX-summary.pdf - county level summary boxplots, geographic scale, fraction of variance due to state or county, and correlation among indicators by subcategory XXXX

(10) Source CVI_analysis.R to make k-means clustering figures. Creates the following files:
- CVI-pct/CVI-pct-comb-clusters.csv - overall CVI ToxPI results with cluster assignments
- Figures/k-means clustering.pdf - Figure showing distribution of CVI scores by cluster, ToxPI propfile for mean of each cluster, contributions to variance by state or county for each cluster, and maps of individual clusters
- SuppFigures/PCA.pdf - PCA plots for first three components (PC1 vs PC2 and PC1 vs PC3) colored by cluster

(11) Source CVI_county_zoom.R to make Harris county zoom in figures. Creates the following files:
- Figures/Harris County Zoom.pdf - figure showing CVI score choropleth map for Harris county, callouts with ToxPI profiles for top three census tract, and category ToxPI profiles for top scoring census tract
- SuppFigures/Harris County TopTract scores.pdf - figure showing bar charts of individual indicators for top scoring census tract in Harria county

(12) Source CVI_CEJST_comparsion.R to do comparison between CVI and CEJST disadvantaged communities. Creates the following files:
- SuppFigures/CVI_CEJST_boxplot.pdf - boxplot comparing CVI scores for CEJST disadvantaged vs. not disadvantaged communities

(13) Source plot_explanatory_vars.R to run correlation analyses between CVI overall categories (total score, baseline and climate change), and a series of explanatory variables (land use, per-capita income, population density). In addition, plot results summarized by two datasets of interest: HOLC grade, and EPA region. 

Uses CEJST.csv file which contains the census tracts with whether it is CEJST disadvantaged 

Note: Other files in Data were used to create "internal points" file for census tracts (since Gaz files sometimes had internal points on or outside the boundary of the census tract)

