#---IMPORTS------------------------------------------------------------------------------------------------

from osgeo import gdal
from netCDF4 import Dataset
from matplotlib import cm
from datetime import datetime
import cartopy
import cartopy.crs as ccrs  # for plotting in cartographic projection
import cartopy.io.shapereader as shpreader # Import shapefiles
import matplotlib.pyplot as plt  # the plotting interface
import numpy as np  # for numerical array computations
import os

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

from cloudsat_read import get_geodata, read_data
from INPE_utilities import download_CMI, reproject, loadCPT

#---VARIAVEIS-&-PREPARATIVOS----------------------------------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
# input_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Dados'
# output_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Figuras'
input_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
output_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'

# arquivos a serem abertos
cloudsat_fname = os.path.join(input_, 'frente_oceanica_2B-GEOPROF.h5')
shapefile = list(shpreader.Reader(os.path.join(input_, 'BR_UF_2019.shp')).geometries())

# recorte da area de estudo (caso continental)
# lat_min = -40
# lat_max = -20
# lon_min = -65
# lon_max = -40

# recorte da area de estudo (caso frente oceanica)
lat_min = -55
lat_max = -15
lon_min = -60
lon_max = -15

extent = [lon_min, lon_max, lat_min, lat_max] # South America

# data da imagem de satelite
# continental
# yyyymmdd = "20190224" #
# hhmn = "1745"

# oceanica
yyyymmdd = "20190701" #
hhmn = "1600"
_datetime = yyyymmdd + hhmn

#---PREPARATIVOS DA IMAGEM DE SATELITE-----------------------------------------------------------------------------------------------------
# file_ir = download_CMI(_datetime, 13, input_) # baixa a imagem, se nao existir
file_ir =  'OR_ABI-L2-CMIPF-M6C13_G16_s20191821600278_e20191821609597_c20191821610071' # pra frente oceanica
var = 'CMI'

# abre o arquivo
img = gdal.Open(f'NETCDF:{input_}/{file_ir}.nc:' + var)

# Read the header metadata
metadata = img.GetMetadata()
scale = float(metadata.get(var + '#scale_factor'))
offset = float(metadata.get(var + '#add_offset'))
undef = float(metadata.get(var + '#_FillValue'))
dtime = metadata.get('NC_GLOBAL#time_coverage_start')

# Load the data
ds_cmi = img.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize).astype(float)

# Apply the scale, offset and convert to celsius
ds_cmi = (ds_cmi * scale + offset) - 273.15

# Reproject the file
filename_ret = f'{input_}/IR_{_datetime}.nc'
reproject(filename_ret, img, ds_cmi, extent, undef)

# Open the reprojected GOES-R image
file = Dataset(filename_ret)

# Get the pixel values
data = file.variables['Band1'][:]

#---PREPARATIVOS DO CLOUDSAT-----------------------------------------------------------------------------------------------------

# Le a passagem do cloudsat e a elevacao ao longo do caminho.
cloudsat_lons, cloudsat_lats, cloudsat_height, cloudsat_time, elev = get_geodata(
     cloudsat_fname, return_list = True
)   

#---PLOTAGEM-----------------------------------------------------------------------------------------------------

fig = plt.figure(figsize=(15,15))
ax = plt.axes(projection=ccrs.PlateCarree())

# Converte um arquivo CPT para ser usado no python
cpt = loadCPT(os.path.join(input_, 'IR4AVHRR6.cpt'))
cmap = cm.colors.LinearSegmentedColormap('cpt', cpt) 
    
# Plot the image
img1 = ax.imshow(data, origin='upper', vmin=-103.0, vmax=84, extent=extent, cmap=cmap, alpha=1.0)

# plotas as coastlines, fronteiras e gridline.
ax.coastlines(resolution='50m', color='white', linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.5)
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='white',facecolor='none', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Plot da trajetoria
ax.plot(cloudsat_lons, cloudsat_lats, color = "tab:red", linewidth=4, transform=ccrs.PlateCarree());
ax.set_extent(extents = extent, crs = ccrs.PlateCarree())

# Plota a colorbar
plt.colorbar(img1, label='Brightness Temperatures (°C)', extend='both', orientation='horizontal', pad=0.03, fraction=0.05)

# Extrai a data
date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))

# Adiciona um titulo
plt.title('GOES-16 Band 13 ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Reg.: ' + str(extent) , fontsize=10, loc='right')

# Salva a imagem
plt.savefig(os.path.join(output_, 'Satelite_oceanica_cloudsat_trajetoria.jpg'), bbox_inches='tight', dpi=300)
plt.close()
