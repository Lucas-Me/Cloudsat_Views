import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr


#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
input = './Dados'
output = './Figuras'

# propriedades dos arquivos
time_ = "2019-02-24T18:00:00"
var = 'z' # variavel a ser usada = geopotencial
g = 9.8 # aceleracao da gravidade

# nome dos arquivos
pnmm_name = 'fdn_pnmm.nc'
multi_name = 'fdn_multi.nc'

# recorte da area de estudo
#lat_min = -70.00
#lat_max = 20.00
#lon_min = -100.00
#lon_max = -10.00

# recorte zoom da area de estudo
lat_min = -60.00
lat_max = -20.00
lon_min = -70.00
lon_max = -30.00

extent = [lon_min, lon_max, lat_min, lat_max] # South America

#---Manipulação dos dados-------------------------------------------------------------------------
pnmm = xr.open_dataset(os.path.join(input, pnmm_name))
pnmm = pnmm.assign_coords(dict(longitude = (((pnmm.longitude.values + 180) % 360) - 180)))
pnmm = pnmm.sortby(pnmm.longitude)

multi = xr.open_dataset(os.path.join(input, multi_name))[var]
multi = multi.assign_coords(dict(longitude = (((multi.longitude.values + 180) % 360) - 180)))
multi = multi.sortby(multi.longitude)

# latitudes e longitudes
lats = multi.latitude.values
lons = multi.longitude.values
xx, yy = np.meshgrid(lons, lats)

# convertendo as variaveis para unidades de hPa
mslp = pnmm['msl'].sel(dict(time = time_)) / 100

# dados nos niveis de pressao
hgt_1000 = multi.sel(dict(level = 1000, time = time_))/g # converte para altura geopotencial
hgt_500 = multi.sel(dict(level = 500, time = time_))/g # converte para altura geopotencial

# calcula e suaviza a espessura entre 1000 e 500
espessura_1000_500 = (hgt_500 - hgt_1000)

#---Projecao do mapa e dos dados-------------------------------------------------------------------------
mapproj = ccrs.PlateCarree()

dataproj = ccrs.PlateCarree()

#---Criar figura e plotar os dados-------------------------------------------------------------------------
fig = plt.figure(1, figsize = (11, 11), dpi = 200)
ax = plt.subplot(111, projection=mapproj)

# Adiciona o shapefile do brasil
shapefile = list(shpreader.Reader('./Dados/BR_UF_2019.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black',facecolor='none', linewidth=0.5)

# Adiciona coastlines, bordas e linhas de grade
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 10), ylocs=np.arange(-90, 90, 10), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Adiciona uma mascara de continente
ax.add_feature(cfeature.LAND)
ax.set_extent(extent) 

# define os intervalos da legenda
step = 60
clevs = np.arange(4680, 6120 + step, step)

# lista de cores, em ordem crescete. RGBA
colors = np.array([
    [67,99,155,255],
    [102,164,239,255],
    [167,230,254,255],
    [255,255,255,255], 
    [249,219,132,255],
    [252,137,98,255],
    [165,63,66,255]
])

# cria um novo cmap a partir do pre-existente
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors/255, clevs.shape[0] - 1)

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(clevs, cmap.N)
kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
              'rightside_up': True, 'use_clabeltext': True}

# plota o mapa de espessura
cs = ax.contourf(lons, lats, espessura_1000_500, levels=clevs, cmap = cmap, 
                norm = norm, transform=dataproj)

# adiciona uma colorbar
plt.colorbar(cs, label='Espessura 1000-500 hPa (m)', extend='both', orientation='vertical',
            shrink = 0.95, pad=0.05, fraction=0.05)

# plota o PNMM
clevmslp = np.arange(800., 1120., 2)
cs2 = ax.contour(lons, lats, mslp, clevmslp, colors='k', linewidths=1.25,
                 linestyles='solid', transform=dataproj)
plt.clabel(cs2, **kw_clabels)

# Put on some titles
plt.title('PNMM (hPa) & Espessura 1000-500 hPa (m)', loc='left')
plt.title('{}'.format(time_), loc='right')

plt.savefig(os.path.join(output, f'pnmm_espessura_{time_}_zoom.jpeg'))
plt.close()