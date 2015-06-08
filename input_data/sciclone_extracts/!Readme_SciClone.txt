The files in this directory are the extracted datasets produced
using the SciClone computing cluster at William and Mary.

Because the underlying datasets to produce this data are terrabytes in
scale, they are not included in this distribution.

The process used to produce the extracts followed a three-step process:

(A) Datasets were downloaded and converted to a GeoTiff format, keeping
native projections and raster resolutions.

(B) For each community boundary, the average pixel value was calculated
for a given raster.  This followed the "small" option in the raster
extraction toolkit within R.

(C) These values were saved to a CSV for joining into the master database,
done via the script "KFW_dataMerge.r"
