# imports
import cartopy.io.shapereader as shpreader # Import shapefiles
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import xarray as xr
import numpy as np

# abrindo netcdf
multi = xr.open_dataset('./Dados/fdn_multi.nc')

# Constantes
lat_min = -60.00
lat_max = 18.00
lon_min = -93.0
lon_max = -25.00
level = 850
time = '2019-02-24T18:00:00'
var_ = 'q' # umidade especifica no netcdf

# Convertendo intervalo de longitude e recorta para a area de estudo.
multi = multi.assign_coords(dict(longitude = (((multi.longitude.values + 180) % 360) - 180)))
multi = multi.sel(dict(
  longitude = slice(lon_min, lon_max),
  latitude = slice(lat_max, lat_min) # array latitude esta em ordem decrescente
)) 

# filtrando a variavel e selecioando o instante de tempo e nivel vertical.
umidade = multi[var_].sel(dict(
    time = time,
    level = level
))*1000 # conversao de kg/kg para g/kg

# Selecionar a extensao do plot [min. lon, max. lon, min. lat, max. lat]
extent = [lon_min, lon_max, lat_min, lat_max] # South America
       
# Tamanho do plot (largura x altura)
# por padrao o dpi = 100, logo a resolucao da figura sera (largura * dpi, altura * dpi)
plt.figure(figsize=(15,15))

# Selecionar a projecao do Cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Adicionar um shapefile
'''
Estados Brasil -> link:
  https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2019/Brasil/BR/br_unidades_da_federacao.zip
'''
shapefile = list(shpreader.Reader('./Dados/BR_UF_2019.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',facecolor='none', linewidth=0.5)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Add a land mask
ax.add_feature(cfeature.LAND)

# Plota a imagem
img = ax.imshow(umidade, extent=extent, cmap='RdYlBu')

# Adiciona uma colorbar
plt.colorbar(img, label='Umidade específica', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

# Adiciona um titulo
plt.title(f'Umidade Específica 850 hPa', fontweight='bold', fontsize=13, loc='center')

# Salva a imagem
plt.savefig('./Figuras/Teste.jpg')
