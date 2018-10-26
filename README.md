# hydro-networks

*Note: hydromodeller.py is deprecated. Use [hydro-networks.ipynb](https://github.com/carderne/hydro-networks/blob/master/hydro-networks.ipynb) instead.*

**hydro-networks** is a network-aware hydrological modelling tool using GeoPandas and rasterio, aimed at country-scale scriptable hydrological modelling.

The Jupyter notebook is split into several sections:

1. Create network arcs and nodes based on input GIS river network
2. Apply Shreve stream order to the network
3. Use [rasterio](https://github.com/mapbox/rasterio) to extract attributes from several raster layers
4. Calculate rainfall runoff from precipitatin
5. Use the network to calculate discharge contributions downstream
6. Estimate mini-hydropower potential in the network based on some parameters
7. Save the results as several geopackages

There's some more detail and pretty pictures in my blog post: [Modelling hydrological networks at massive scale](https://rdrn.me/modelling-hydrological-networks/) 