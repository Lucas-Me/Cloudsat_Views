# imports
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader # Import shapefiles

# Recorte para a area de estudo
lat_min = -60.00
lat_max = 18.00
lon_min = -93.0
lon_max = -25.00

# Convertendo intervalo de longitude e recorta para a area de estudo.
multi = xr.open_dataset('./Dados/fdn_multi.nc')
multi = multi.assign_coords(dict(longitude = (((multi.longitude.values + 180) % 360) - 180)))
multi = multi.sel(dict(longitude = slice(lon_min, lon_max), latitude = slice(lat_max, lat_min)))

# filtrando a variavel
umidade = multi['q'].sel(dict(
    time = "2019-02-24T18:00:00",
    level = 850
))*1000

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [lon_min, lon_max, lat_min, lat_max] # South America
       
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(15,15))

# Use the Cilindrical Equidistant projection in cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Add a shapefile
#https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
shapefile = list(shpreader.Reader('./Dados/BR_UF_2019.shp').geometries())
# or
#https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip
#shapefile = list(shpreader.Reader('ne_10m_admin_1_states_provinces.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',facecolor='none', linewidth=0.5)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Add a land mask
ax.add_feature(cfeature.LAND)

# Plot the image
img = ax.imshow(umidade, extent=extent, cmap='RdYlBu')

# Add a colorbar
plt.colorbar(img, label='Umidade específica', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

# Add a title
plt.title(f'Umidade Específica 850 hPa', fontweight='bold', fontsize=13, loc='center')

# Save the image
plt.savefig('./Figuras/Teste.jpg')
plt.close()