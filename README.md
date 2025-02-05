The GeoStats Python package - `geostats`
===========

<a href="https://pypi.python.org/pypi/hmi/">![PyPiVersion](https://img.shields.io/pypi/v/hmi.svg)</a> <a href="">![Pyversions](https://img.shields.io/pypi/pyversions/hmi.svg)</a> <a href="https://hmi.readthedocs.io/en/latest/">![ReadTheDocs](https://readthedocs.org/projects/hmi/badge/?version=latest&style=plastic)</a> <a>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14291903.svg)](https://doi.org/10.5281/zenodo.14291903)</a>

The geostats package is a python module that provides an interface to compute spatial statistics based on a shapefile for various datasets. The package uses a user-provided shapefile or `GeoDataFrame` to compute statistics for each polygon. 

Statistics
----

Default statistics computed are `'min', 'max', 'median', 'mean', 'majority', 'sum','std','count'`. Additional statistics can be computed by defining a function to compute it and passing the function to the package using the `addstats` keyword. E.g., the package includes the additional statistic `gini` defined by

```python
def gini(self):
    '''
    gini(geo)

    Return computed Gini coefficient.
    '''
    if self.count()>1:
        xsort = sorted(self.data[self.mask==False].flatten()) # increasing order
        y = np.cumsum(xsort)
        B = sum(y) / (y[-1] * len(xsort))
        return 1 +  1./len(xsort) - 2*B
    else:
        return 1
```
which is then passed using `addstats={'gini':gini}` keyword.

Usage example
------

```python
# Import packages
import geopandas
import geostats
import requests
import io

# Download countries shapefile from NaturalEarth
headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36', 'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'}

url = 'https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_countries.zip'
r = requests.get(url, headers=headers)
countries = gp.read_file(io.BytesIO(r.content))

# Compute statistics for all measures
Stats = geostats.geostats(countries, stats=['mean','max','min','std', 'sum'])
Stats.geostats()
print(Stats.df)

# Compute statistics for only one or a few measures 
mymeasures = ['CSI', 'Sea100']
Stats = geostats.geostats(countries, stats=['mean','max','min','std', 'sum'], measures=mymeasures)
Stats.geostats()
print(Stats.df)

```

`Stats.df` is a `GeoDataFrame` that contains the  information from the original `GeoDataFrame` or `shapefile` and all the requested statistics for each measure. Naming conventions are such that statistic `s` for measure `m` appears as column `ms`. Below are the names for each measure in the database.

Datasets
------
The package computes spatial statistics using the following datasets:


| Variable Type | Source/Citation	| License	| geostats name	|	Location	|
|--------|--------|-------|--------------|---------|
| <tr><td colspan="5" align="center"><strong>Geography, Climate, and Agricultural Data</strong></td></tr>  | | | | |
| Mobility| Özak, Ömer. “[Distance to the pre-industrial technological frontier and economic development.](http://rdcu.be/I4YI)” Journal of Economic Growth 23.2 (2018): 175-221. Özak, Ömer. “[The voyage of homo-economicus: Some economic measures of distance.](http://omerozak.com/pdf/Ozak_voyage.pdf)” Department of Economics, Brown University (2010). [(Data)](https://doi.org/10.5281/zenodo.14285746) | This data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License.](https://creativecommons.org/licenses/by-sa/4.0/)| `HMI`, `HMI10`, `HMISea`, `HMISea10`| `/geostats-data/HMI/tifs/`| 
|GAEZ Crop Yield| IIASA/FAO, 2012. Global Agro‐ecological Zones (GAEZ v3.0). IIASA, Laxenburg, Austria and FAO, Rome, Italy. Fischer, G., Nachtergaele, F.O., Prieler, S. , Teixeira, E. , Toth, G., van Velthuizen, H., Verelst, L., & Wiberg, D. (2012). [Global Agro-ecological Zones (GAEZ v3.0)- Model Documentation.](https://pure.iiasa.ac.at/id/eprint/13290/1/GAEZ_Model_Documentation.pdf) IIASA, Laxenburg, Austria and FAO, Rome, Italy. , Laxenburg, Austria; Rome, Italy. | [License](https://www.fao.org/fileadmin/user_upload/gaez/docs/User_Agreement_and_Disclaimer_EN.pdf) states "Reproduction and dissemination of material contained in GAEZ v3.0 or educational, research, personal or other noncommercial purposes are authorized without any prior written permission from the copyright holders, provided FAO and IIASA are fully acknowledged." | Crop name followed by input level given by `hi`, `med` or `lo`| `/geostats-data/CropsYield/tifs/`| 
|Crop Caloric Suitability| Oded Galor and Ömer Özak, 2016. “[The Agricultural Origins of Time Preference](http://dx.doi.org/10.1257/aer.20150020),” American Economic Review, 2016, 106(10): 3064–3103. Oded Galor and Ömer Özak, 2015. “[Land Productivity and Economic Development: Caloric Suitability vs. Agricultural Suitability](http://papers.ssrn.com/abstract=2625180),” Brown University Working Paper. Özak, Ö. (2015). [Caloric Suitability Index - Data (v1.0) [Data set]](https://doi.org/10.5281/zenodo.14714917). Zenodo. |This data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License.](https://creativecommons.org/licenses/by-sa/4.0/)| `CSI` followed by crop name and input level given by `hi`, `med` or `lo`| `/geostats-data/CSICrops/tifs/`| 
|Caloric Suitability Index| Oded Galor and Ömer Özak, 2016. “[The Agricultural Origins of Time Preference](http://dx.doi.org/10.1257/aer.20150020),” American Economic Review, 2016, 106(10): 3064–3103. Oded Galor and Ömer Özak, 2015. “[Land Productivity and Economic Development: Caloric Suitability vs. Agricultural Suitability](http://papers.ssrn.com/abstract=2625180),” Brown University Working Paper. Özak, Ö. (2015). [Caloric Suitability Index - Data (v1.0) [Data set]](https://doi.org/10.5281/zenodo.14714917). Zenodo. |This data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License.](https://creativecommons.org/licenses/by-sa/4.0/)| `PeriodMeasure` where `Period` is `pre1500`: Pre-1500 crop-based data; `pre15002`: Pre-1500 crop-based data excluding Asian crops in Africa; `post1500`: Post-1500 crop-based data; `dif`: Difference between post-1500 and pre-1500; `dif2`: Difference between post-1500 and pre-15002; `Measure` is `Average`, `Maximum`, `Total`, or other statistic of caloric suitability across crops | `/geostats-data/CSI/tifs/`| 
|Caloric Plow Suitability| Galor, Oded, Ömer Özak and Assaf Sarid, “[Geographical Origins and Economic Consequences of Language Structures](http://ssrn.com/abstract=2820889)” Brown University Working Paper, 2016. Oded Galor and Ömer Özak, 2016. “[The Agricultural Origins of Time Preference](http://dx.doi.org/10.1257/aer.20150020),” American Economic Review, 2016, 106(10): 3064–3103. Oded Galor and Ömer Özak, 2015. “[Land Productivity and Economic Development: Caloric Suitability vs. Agricultural Suitability](http://papers.ssrn.com/abstract=2625180),” Brown University Working Paper. Özak, Ö. (2015). [Caloric Suitability Index - Data (v1.0) [Data set]](https://doi.org/10.5281/zenodo.14714917). Zenodo. |This data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License.](https://creativecommons.org/licenses/by-sa/4.0/)| `CropTypePeriodMeasure` where `CropType` is `plowpos` (plow positive), `plowneg` (plow negative), `plowpot` (plow potential: difference between positive and negative); `Period` is `pre1500`: Pre-1500 crop-based data; `pre15002`: Pre-1500 crop-based data excluding Asian crops in Africa; `post1500`: Post-1500 crop-based data; `dif`: Difference between post-1500 and pre-1500; `dif2`: Difference between post-1500 and pre-15002;  `Measure` is `en`, `enmean`, `entot`, or other statistic of caloric suitability across crops | `/geostats-data/CSIPlow/tifs/`| 
|Caloric Suitability and Growth Cycle of best crop| Oded Galor and Ömer Özak, 2016. “[The Agricultural Origins of Time Preference](http://dx.doi.org/10.1257/aer.20150020),” American Economic Review, 2016, 106(10): 3064–3103. Oded Galor and Ömer Özak, 2015. “[Land Productivity and Economic Development: Caloric Suitability vs. Agricultural Suitability](http://papers.ssrn.com/abstract=2625180),” Brown University Working Paper. Özak, Ö. (2015). [Caloric Suitability Index - Data (v1.0) [Data set]](https://doi.org/10.5281/zenodo.14714917). Zenodo. |This data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License.](https://creativecommons.org/licenses/by-sa/4.0/)| `PeriodMethodMeasure` where `Period` is `pre1500`: Pre-1500 crop-based data; `pre15002`: Pre-1500 crop-based data excluding Asian crops in Africa; `post1500`: Post-1500 crop-based data; `dif`: Difference between post-1500 and pre-1500; `dif2`: Difference between post-1500 and pre-15002; `Methods` is how the best crop was chosen: `en` (calories), `inv` (return), `cyc` (cycle), `hvst` (harvests), `winter` (unproductive period), etc. followed by `statistic`.| `/geostats-data/CSICycle/tifs/`| 
|Elevation|NOAA National Centers for Environmental Information. 2022: [ETOPO 2022 15 Arc-Second Global Relief Model.](https://www.ncei.noaa.gov/products/etopo-global-relief-model) NOAA National Centers for Environmental Information. DOI: 10.25921/fd45-gt74. Accessed January 16, 2025. [Documentation](https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/docs/1.2%20ETOPO%202022%20User%20Guide.pdf)| [Creative Commons Zero 1.0 Universal Public Domain Dedication (CC0-1.0).](https://creativecommons.org/publicdomain/zero/1.0/deed.en) according to [their website](https://data.noaa.gov/metaview/page?xml=NOAA/NESDIS/NGDC/MGG/DEM//iso/xml/etopo_2022.xml&view=getDataView&header=none)| `etopo`, `etopocyl`| `/geostats-data/ETOPO/`|
|Elevation|GLOBE Task Team and others (Hastings, David A., Paula K. Dunbar, Gerald M. Elphingstone, Mark Bootz, Hiroshi Murakami, Hiroshi Maruyama, Hiroshi Masaharu, Peter Holland, John Payne, Nevin A. Bryant, Thomas L. Logan, J.-P. Muller, Gunter Schreier, and John S. MacDonald), eds., 1999. [Global land one-kilometer base elevation (GLOBE)](https://www.ngdc.noaa.gov/mgg/topo/report/globedocumentationmanual.pdf) Digital Elevation Model, Version 1.0. National Oceanic and Atmospheric Administration, National Geophysical Data Center, 325 Broadway, Boulder, Colorado 80303, U.S.A. Digital data base on the World Wide Web (URL: http://www.ngdc.noaa.gov/mgg/topo/globe.html) and CDROMs. | Can be freely downloaded from the [GLOBE website](https://www.ngdc.noaa.gov/mgg/topo/gltiles.html) and redistributed according to [documentation](https://www.ngdc.noaa.gov/mgg/topo/report/globedocumentationmanual.pdf)|`globe`, `globecyl`|`/geostats-data/GLOBE/`|
|Ruggedness| Based on GLOBE data following methodology in Riley, S. J., DeGloria, S. D., & Elliot, R. (1999). [Index that quantifies topographic heterogeneity.](http://download.osgeo.org/qgis/doc/reference-docs/Terrain_Ruggedness_Index.pdf) intermountain Journal of sciences, 5(1-4), 23-27.| [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en) | `rix`, `rixcyl`| `/geostats-data/RIX/`|
|Uban area, crop area, gras area| Klein Goldewijk, K., Beusen A., van Drecht, G., de Vos, M., 2011. [The HYDE 3.1 spatially explicit database of human induced land use change over the past 12,000 years](https://doi.org/10.1111/j.1466-8238.2010.00587.x), Global Ecology and Biogeography 20(1): 73-86. Klein Goldewijk, K., Beusen, A., Janssen, P., 2010. [Long term dynamic modeling of global population and built-up area in a spatially explicit way, HYDE 3 .1.](https://doi.org/10.1177/0959683609356587) The Holocene, 20(4), 565-573. | Freely available on the web [here](https://public.yoda.uu.nl/geo/UU01/8K9D7F.html). License is [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)|`uopp_XXXXX`, `grasXXXX`, `cropXXXX`, where `XXXX=10000BC, 9000BC, ..., 1000AD, ..., 2005AD`|`/geostats-data/HYDE/tifs/`|
| Climate (1901-2012) | Harris, I. P. D. J., Jones, P., Osborn, T., & Lister, D. (2014). [Updated high-resolution grids of monthly climatic observations-the CRU TS3. 10 Dataset.](https://doi.org/10.1002/joc.3711) International journal of climatology, 34, 623-642. | [Contains public sector information licensed under the Open Government Licence v3.0.](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) [(License)](https://dap.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_3.21/00README_catalogue_and_licence.txt?download=1)| Rasters have been aggregated to show the mean, volatility, and spatial correlation (mean and min) for each cell across all years. The names of the different climatic variables are `pre` - Precipitation (mm/month), `tmp` - Mean Temperature (°C), `tmn` - Minimum Temperature (°C), `tmx` - Maximum Temperature (°C), `dtr` - Diurnal Temperature Range (°C), `reh` - Relative Humidity (%), `pet` - Potential Evapotranspiration (mm/month), `cld` - Cloud Cover (%), `frs` - Frequency of Frost Days (days/month), `vap` - Vapor Pressure (hPa), `wet` - Wet Days (days/month). So, e.g., for `tmp` the measure names are `tmpmean`, `tmpspatcorr`, `tmpvolatility` followed by the statistic name and `cyl` if it is uses the cylindrival equal area projection. |`/geostats-data/CRU/tif/`|
|Malaria Ecology|Kiszewski, A., Mellinger, A., Spielman, A., Malaney, P., Sachs, S. E., & Sachs, J. (2004). [A global index representing the stability of malaria transmission.](https://www.academia.edu/download/94497521/486.pdf) American Journal of tropical medicine and hygiene, 70(5), 486-498. [https://doi.org/10.4269/ajtmh.2004.70.486](https://doi.org/10.4269/ajtmh.2004.70.486)| Freely available on the web [here](https://www.dropbox.com/scl/fi/t4je3iqo795gdu18v4swr/ME_raster.zip?rlkey=921dnfkahuxh4g4jy0uh8mmdl&e=2&dl=0). | `stxv5`, `stxv5cyl`|`/geostats-data/Malaria/malaria_ecology/`|
| Tse-Tse (Africa, presence probability)| Wint, W. & Rogers, D. (2000) [Predicted Distributions of Tsetse in Africa. Consultancy Report for the Animal Health Service of the Animal Production and Health Division of the Food and Agriculture Organization of the United Nations. FAO, Rome.](https://openknowledge.fao.org/items/956f7aad-64e2-4bff-af3b-623b2215587c) [Accessed on Jan. 14, 2015] [(Data)](https://data.apps.fao.org/map/catalog/srv/eng/catalog.search#/metadata/f8a4e330-88fd-11da-a88f-000d939bc5d8)| FAO shares data under various versions of the Creative Commons Attribution-NonCommercial-ShareAlike. See, e.g., [OPEN DATA LICENSING FOR STATISTICAL DATABASES](https://openknowledge.fao.org/server/api/core/bitstreams/6fe33207-e53d-4060-8d81-396b323789d5/content) or this [related data](https://data.apps.fao.org/map/catalog/srv/eng/catalog.search#/metadata/377805c0-1e4d-11dc-abdf-000d939bc5d8). | `Fuscacyl`, `Morsitanscyl`, `Palpaliscyl`, `tsetsecyl` | `/geostats-data/tsetse/` |
|Ecological Diversity (based on 827 ecological zones around year 2000)| Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V., Underwood, E. C., D’amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C. et al. (2001). [Terrestrial ecoregions of the world: A new map of life on earth a new global map of terrestrial ecoregions provides an innovative tool for conserving biodiversity](https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2), BioScience 51(11): 933–938. [(Old website)](https://www.worldwildlife.org/biomes) [(Data)](https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world)| PENDING CLARIFICATION: If it cannot be shared, then either user will have to download [this file](https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip) and put in right place, or the script will have to download and unzip [this file](https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip)| `ecodiversitywwf`, `ecopolarizationwwf`| `/geostats-data/Ecological_Zones/`|
|Ecological Diversity (based on 16 biomes around year 2000)| Aggregation of ecological zones in Olson et al. (2001) | PENDING CLARIFICATION: If it cannot be shared, then either user will have to download [this file](https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip) and put in right place, or the script will have to download and unzip [this file](https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip) and then aggregate biomes. | `ecodiversity`, `ecopolarization`| `/geostats-data/Ecodiversity/`|
|Ecological Diversity (during Last Glacial Maximum LGM)|Ray, Nicolas, Adams, Jonathan. [A GIS-based vegetation map of the world at the Last Glacial Maximum (25,000-15,000 BP)](https://www.ncei.noaa.gov/access/paleo-search/study/6220). In: *Internet archaeology*, 2001, vol. 11. (Accessed on Jun. 5, 2017) [(Paper)](https://intarch.ac.uk/journal/issue11/rayadams_toc.html) [(Data)](https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/ray2001/)| NOOA Website says "Use Constraints: Please cite original publication, online resource, dataset and publication DOIs (where available), and date accessed when using downloaded data. If there is no publication information, please cite investigator, title, online resource, and date accessed." *Internet Archaeology* is Open-Access and [it states](https://intarch.ac.uk/about/index.html) that "Except where otherwise noted, content from this work may be used under the terms of the [Creative Commons Attribution 3.0 (CC BY) Unported licence](http://creativecommons.org/licenses/by/3.0/), which permits unrestricted use, distribution, and reproduction in any medium, provided that attribution to the author(s), the title of the work, the Internet Archaeology journal and the relevant URL/DOI are given."| `ecodiversityLGM`, `ecopolarizationLGM`| `/geostats-data/EcoDiversityLGM/`
|Coast Length|Made with [Natural Earth](https://www.naturalearthdata.com). Free vector and raster map data @ naturalearthdata.com.|[Creative Commons Public Domain License](https://www.naturalearthdata.com/about/terms-of-use/).| `coastlen` |`/geostats-data/Waters-Coasts/`|
|Land within 100kms of Sea|Made with [Natural Earth](https://www.naturalearthdata.com). Free vector and raster map data @ naturalearthdata.com.|[Creative Commons Public Domain License](https://www.naturalearthdata.com/about/terms-of-use/).| `sea100`, `sea100pct` |`/geostats-data/Waters-Coasts/`|
| Agricultural Suitability| Ramankutty, N., J.A. Foley , J. Norman, and K. McSweeney. [The global distribution of cultivable lands: current patterns and sensitivity to possible climate change.](https://doi.org/10.1046/j.1466-822x.2002.00294.x) Submitted to Global Ecology and Biogeography, March 2001 | Originally data freely available at [http://atlas.sage.wisc.edu/](http://atlas.sage.wisc.edu/). Updated website (with limited data) [https://sage.nelson.wisc.edu/data-and-models/atlas-of-the-biosphere/mapping-the-biosphere/land-use/suitability-for-agriculture/](https://sage.nelson.wisc.edu/data-and-models/atlas-of-the-biosphere/mapping-the-biosphere/land-use/suitability-for-agriculture/)|`climfac`, `climsuit`, `soilfac`, `suit`, `climfaccyl`, `climsuitcyl`, `soilfaccyl`, `suitcyl` | `/geostats-data/Ramankutty/` |
| <tr><td colspan="5" align="center"><strong>Lights, Population, and Other</strong></td></tr>  | | | | |
|Lights| F. C. Hsu, K. Baugh, T. Ghosh, M. Zhizhin, and C. Elvidge, “[DMSP-OLS Radiance Calibrated Nighttime Lights Time Series with Intercalibration](https://eogdata.mines.edu/products/dmsp/),” Remote Sensing, vol. 7, pp. 1855–1876, 2015. Image and Data processing by [NOAA's National Geophysical Data Center](https://www.ngdc.noaa.gov/eog/index.html). DMSP data collected by the US Air Force Weather Agency. | Original data can be downloaded from [this website](https://www.ngdc.noaa.gov/eog/dmsp/downloadV4composites.html). Data now housed at the [Earth Observation Group](https://payneinstitute.mines.edu/eog/). License is [creative commons 4](https://eogdata.mines.edu/files/EOG_products_CC_License.pdf), and can be redistributed. | `FXXYEAR`, `FXXYEARcyl` where `XX=10,...,18` and `YEAR=1992,...,2012` | `/geostats-data/Lights/tifs/`|
|Lights| Li, X., Zhou, Y., Zhao, M., & Zhao, X. (2020). [A harmonized global nighttime light dataset 1992–2018](https://doi.org/10.1038/s41597-020-0510-y). Scientific data, 7(1), 168. | Original data can be downloaded from [this website](https://doi.org/10.6084/m9.figshare.9828827.v2). License is [creative commons 4](https://figshare.com/articles/dataset/Harmonization_of_DMSP_and_VIIRS_nighttime_light_data_from_1992-2018_at_the_global_scale/9828827/2), and can be redistributed. | `HLDYEAR` where `YEAR=1992,...,2020` | `/geostats-data/HLD/tifs/`|
|Population|Center For International Earth Science Information Network-CIESIN-Columbia University, United Nations Food And Agriculture Programme-FAO, & Centro Internacional De Agricultura Tropical-CIAT. (2005). [Gridded Population of the World, Version 3 (GPWv3): Population Count Grid (Version 3.00)](https://doi.org/10.7927/H4639MPP) [Data set]. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H4639MPP](https://doi.org/10.7927/H4639MPP) (Accessed Mar. 25, 2015)| Users are prohibited from any commercial, non-free resale, or redistribution. | `glupXXag`, `glupXXg`, `glpXXag`, `glpXXg`, where `XX=90, 95, 00` is the year | `/geostats-data/GPW/v3/population/`|
|Population|Center For International Earth Science Information Network-CIESIN-Columbia University. (2018). [Gridded Population of the World, Version 4 (GPWv4): Population Count, Revision 11 (Version 4.11)](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-popcount-r11-4.11) [Data set]. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H4JW8BX5](https://doi.org/10.7927/H4JW8BX5) (Accessed Jan. 17, 2025)| This work is licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0). Users are free to use, copy, distribute, transmit, and adapt the work for commercial and non-commercial purposes, without restriction, as long as clear attribution of the source is provided. | `gpwv4popcXXXX` where `XXXX=2000, 2005, 2010, 2015, 2020` is the year | `/geostats-data/GPW/v4/population/`|
|Population|Center For International Earth Science Information Network-CIESIN-Columbia University. (2018). [Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of UN WPP Country Totals, Revision 11 (Version 4.11)](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-apct-wpp-2015-r11-4.11) [Data set]. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H4PN93PB](https://doi.org/10.7927/H4PN93PB) (Accessed Jan. 17, 2025)|This work is licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0). Users are free to use, copy, distribute, transmit, and adapt the work for commercial and non-commercial purposes, without restriction, as long as clear attribution of the source is provided. | `gpwv4popcadjXXXX` where `XXXX=2000, 2005, 2010, 2015, 2020` is the year | `/geostats-data/GPW/v4/population/`|
|Population Density|Center For International Earth Science Information Network-CIESIN-Columbia University, & Centro Internacional De Agricultura Tropical-CIAT. (2005). [Gridded Population of the World, Version 3 (GPWv3): Population Density Grid (Version 3.00)](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv3-popdens-3.00) [Data set]. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H4XK8CG2](https://doi.org/10.7927/H4XK8CG2) (Accessed Mar. 25, 2015)| Users are prohibited from any commercial, non-free resale, or redistribution. | `gludsXXag`, `gludsXXg`, `gldsXXag`, `gldsXXg`, where `XX=90, 95, 00` is the year | `/geostats-data/GPW/v3/popdens/`|
|Population Density|Center For International Earth Science Information Network-CIESIN-Columbia University. (2017). [Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11 (Version 4.11)](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-popdens-r11-4.11) [Data set]. Palisades, NY: Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H49C6VHW](https://doi.org/10.7927/H49C6VHW) (Accessed Jan. 17, 2025)|This work is licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0). Users are free to use, copy, distribute, transmit, and adapt the work for commercial and non-commercial purposes, without restriction, as long as clear attribution of the source is provided. | `gpwv4popdXXXX` where `XXXX=2000, 2005, 2010, 2015, 2020` is the year | `/geostats-data/GPW/v4/popdens/`|
|Population Density|Center For International Earth Science Information Network-CIESIN-Columbia University. (2018). [Gridded Population of the World, Version 4 (GPWv4): Population Density Adjusted to Match 2015 Revision UN WPP Country Totals, Revision 11 (Version 4.11)](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-apdens-wpp-2015-r11-4.11) [Data set]. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). [https://doi.org/10.7927/H4F47M65](https://doi.org/10.7927/H4F47M65) (Accessed Jan. 17, 2025)|This work is licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0). Users are free to use, copy, distribute, transmit, and adapt the work for commercial and non-commercial purposes, without restriction, as long as clear attribution of the source is provided. | `gpwv4popdadjXXXX` where `XXXX=2000, 2005, 2010, 2015, 2020` is the year | `/geostats-data/GPW/v4/popdens/`|
|Population (counts, rural/ruban, density)| Klein Goldewijk, K., Beusen A., van Drecht, G., de Vos, M., 2011. [The HYDE 3.1 spatially explicit database of human induced land use change over the past 12,000 years](https://doi.org/10.1111/j.1466-8238.2010.00587.x), Global Ecology and Biogeography 20(1): 73-86. Klein Goldewijk, K., Beusen, A., Janssen, P., 2010. [Long term dynamic modeling of global population and built-up area in a spatially explicit way, HYDE 3 .1.](https://doi.org/10.1177/0959683609356587) The Holocene, 20(4), 565-573. | Freely available on the web [here](https://public.yoda.uu.nl/geo/UU01/8K9D7F.html). License is [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)|`popc_XXXXX`, `urbc_XXXX`, `rurc_XXXX`, `popd_XXXX`, where `XXXX=10000BC, 9000BC, ..., 1000AD, ..., 2005AD`|`/geostats-data/HYDE/tifs/`|
|Population Density (Africa) | [UNESCO (1987) through UNEP/GRID-Sioux Falls](https://na.unep.net/siouxfalls/datasets/datalist.php) (Accessed Jan. 22, 2015)| Most of the UNEP/GRID data sets are freely available for download by users via the Internet. UNEP/GRID does not place any restrictions on the use of this data, but does request that users cite UNEP/GRID. See [license](https://na.unep.net/siouxfalls/datasets/datapolicy.php). | `afpopdXX`, where `XX=60, 70, 80, 90, 00` is the year. | `/geostats-data/UNEP/`|
|Population (counts, 2000-2023)| Bright, E. and P. R. Coleman, 2001. [Landscan 2000-2023](https://landscan.ornl.gov/about), Oak Ridge, TN: Oak Ridge National Laboratory. (Accessed Jan. 17, 2025) [(Data)](https://landscan.ornl.gov/)| [Creative Commons Attribution 4.0 International License](https://landscan.ornl.gov/licensing).| `lsXXXX` where `XXXX=2000,...,20023` |`/geostats-data/Landscan/`|


Install
-------

Given the requirements of the package, the best way to install is to first create a `mamba/conda`  environment as follows (this creates an `python-3.11` environment and adds some basic packages that are needed).

```bash
mamba create --name GeoStats --override-channels -c conda-forge python=3.11 pip geopandas georasters jupyterlab jupyter seaborn geoplot pysal p7zip
```

activate the environment

```bash
mamba activate GeoStats
```

then 

```bash
 pip install geostats
```
 or

```bash
 pip install git+git://github.com/ozak/geostats.git
```

Example Usage: GeoStats
-------------------------

``` python
 import geostats 
``` 

## Citation

If you use the package please cite:

```
Özak, Ömer. 2014. "The GeoStats Python package - `geostats`".
```

Also, make sure to cite all 


 
Issues
------

Find a bug? Report it via github issues by providing

- a link to download the smallest possible raster and vector dataset necessary to reproduce the error
- python code or command to reproduce the error
- information on your environment: versions of python, gdal and numpy and system memory


# Copyright 

&copy; Ömer Özak (2014)

This code and data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License](https://creativecommons.org/licenses/by-sa/4.0/) and [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).
![](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-sa.svg) ![](https://www.gnu.org/graphics/gplv3-127x51.png)
