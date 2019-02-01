# hydronetworks

**hydronetworks** is a network-aware hydrological modelling tool using GeoPandas and rasterio, aimed at country-scale scriptable hydrological modelling.

The Jupyter notebook (`example.ipynb`) is split into several sections:

1. Create network arcs and nodes based on input GIS river network
2. Apply Shreve stream order to the network
3. Use [rasterio](https://github.com/mapbox/rasterio) to extract attributes from several raster layers
4. Calculate rainfall runoff from precipitatin
5. Use the network to calculate discharge contributions downstream
6. Estimate mini-hydropower potential in the network based on some parameters
7. Save the results as several geopackages

There's some more detail and pretty pictures in my blog post: [Modelling hydrological networks at massive scale](https://rdrn.me/modelling-hydrological-networks/)


## Installation

**Requirements**

hydronetworks is confirmed to work with Python >= 3.6 with the following packages installed:
 - `numpy` >= 1.15.4
 - `pandas` >= 0.23.4
 - `shapely` >= 1.6.4
 - `rasterio` >= 1.0.13
 - `geopandas` >= 0.4.0
 - `pyproj` >= 1.9.5.1
 - `matplotlib` >= 3.0.2

**Install from GitHub**

Download or clone the repository and install the required packages (preferably in a virtual environment):

```
git clone https://github.com/carderne/hydronetworks.git
cd hydronetworks
pip install -r requirements.txt
```
You can run ```make test``` in the directory, which will do an entire run through using the test data and confirm whether everything is set up properly.