# Invasive-Grasses-supporting-documents
Short description and code for invasive grasses document

Create seasonal NDVI: Javascript in Google Earth Engine (GEE) was used to create seasonal min/max Normalized Difference Vegetation Index (NDVI) for southern California. 
Script: Generate SoCal NDVI

Generating Plot Locations: R was used to create random plot locations within selected NLCD classes. Stratification based on NDVI seasonal difference and NLCD land cover ensured sufficient inclusion of less dominant land cover types. Polygons were created for each plot that were coincident with the raster cells, in order that future vegetation indices & visual plot summaries encompass identical coordinates. 
Script: ShrubGrass_RandomPoints.R

Imagery Evaluation: For each point, percent grasses was measured within a 30m Landsat cell surrounding each point. Percent grasses was recorded for training years 2009, 2011, 2013, 2015, 2017, & 2018

Statistics of Seasonal Summaries of Vegetation Indices: Plot point locations were uploaded into GEE, and Javascript was used to generate vegetation indices (VI) from Landsat images for 3 seasons of 3 month periods each for each training year at each plot location. Indices were averaged for each season (spring, summer, fall) and statistics (mean, min, max, median, & stdev) were generated for each seasonal measurement for each index. Vegetation indices collected included NDII, NDVI,  ARG (Arctangent of Red Green), NBR (Normalized Burn Ratio), and Band 3.  Spatial records were converted to tabular format before exporting. 
Script: ExtractVIWPts_moreIndices

Develop Model:  In R, a table was assembled containing training point ID, and for each training year; percent grasses, time since most recent fire, most recent fire year (from fire history), and the seasonal VI statistical values.
For each training year 3 statistical methods were applied (linear model, betareg, and random forest) and models were evaluated based on how well they predicted percent herbaceous cover, in order to determine the best predictive model with the least error that could be applied to all years. Covariates used were fire data and the independent variables developed in Earth Engine.  
Script: shrubgrass_findmodel_20200623.R

Create spatial independent variables: Javascript in GEE was used to generate rasters representing the independent variables for the 2009 model, for each training year. 
Script: ExtractVIWPts_moreIndices 

Predict Herbaceous Cover: R was used to apply the best model (in terms of least error and temporal transferability), which was from 2009 training year, to each year's corresponding rasters and predict a herbaceous cover raster for each training year. 
Script: shrubgrass_findmodel_20200623.R 


