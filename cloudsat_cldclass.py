#---IMPORTS----------------------------------------------------------------------------------------
from configparser import Interpolation
import os
from tracemalloc import stop

import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import matplotlib.colors as mcolors
import numpy as np  # for numerical array computations
from scipy.interpolate import NearestNDInterpolator

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

from cloudsat_read import get_geodata, read_data, read_cldclass
import cloudsat_utils

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
# input_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Dados'
# output_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Figuras'
input_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
output_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'

# nome do arquivo cldclass
cldclass_fname = 'frente_continental_2B-CLDCLASS_P1_R05.h5'

# recorte da area de estudo
lat_min = -35
lat_max = -27.5
lon_min = -65
lon_max = -40
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#---Cloudsat Classificacao de Nuvens--------------------------------------------------------------------------

# CldClass do cloudsat
cldclass, cloud_codes = read_cldclass(os.path.join(input_, cldclass_fname))

# dimensoes do dado
cldst_lons, cldst_lats, cldst_height, cldst_time, cldst_elev = get_geodata(
    os.path.join(input_, cldclass_fname), return_list= True
)

# array unidimensional neste produto
cldst_elev = cldst_elev * 1e-3 # convertendo elevacao para km.

# Encontrar indices do array onde o recorte da area esta localizado
ii = np.where(
    (cldst_lons > lon_min)
    & (cldst_lons < lon_max)
    & (cldst_lats > lat_min)
    & (cldst_lats < lat_max)
)[0]
i1, i2 = ii[0], ii[-1]

# indexando a classificacao das nuvens
cldclass = cldclass[i1:i2, :]

# grade na horizontal
cloudsat_x = np.arange(i1, i2, dtype = np.float32)
XX = np.tile(cloudsat_x, (cldclass.shape[1], 1)).T

# grade na vertical
cloudsat_z0 = 0  # km
cloudsat_z1 = 20  # km
cloudsat_nz = 100 # Number of pixels (levels) in the vertical.
cloudsat_z = (cldst_height[i1:i2, :] * 1e-3).astype(np.float32)
new_z = np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz) # novo Z

# preaprando para interpolar -> metodo nearest
x_coords = np.ravel(XX)
z_coords = np.ravel(cloudsat_z)
values = np.ravel(cldclass)

# interpolador
interp = NearestNDInterpolator(list(zip(x_coords, z_coords)), values)
XX, ZZ = np.meshgrid(cloudsat_x, new_z)
cld_class_interpolado = interp(XX, ZZ)

#---Plotando os dados--------------------------------------------------------------------------
print(cloud_codes)
# configurando um colormap, e normalizando o mesmo
fig, ax = plt.subplots(figsize=(22, 6), dpi = 200)

# plot da classifciacao, colorbar ->
clevs = np.linspace(0, 8, 9).astype(np.int32)

# lista de cores, em ordem crescete. [R, G, B, A]
colors = np.array([
    [1, 1, 1, 255], # valor 0 -> clear
    [0, 57, 255, 255], # valor 1 -> Ci
    [0, 162, 255,255], # valor 2 -> As
    [0, 239, 231,255], # Valor 3 -> Ac
    [69, 228, 109, 255], # Valor 4 -> St
    [253, 248, 0, 255], # Valor 5 -> Sc
    [233, 99, 0, 255], # valor 6 -> Cu
    [251, 2, 0, 255], # valor 7 -> Ns
    [177, 0, 5, 255] # valor 8 -> Dc 
])

# cria um novo cmap, para valores inteiros
cld_cmap = mcolors.ListedColormap(colors/255)

# plota os pixeis
LATS = np.take(cldst_lats, XX.astype('int64')) # extrai os valores de latitude respectivo a cada indice em XX
p = ax.pcolormesh(
    LATS, ZZ, cld_class_interpolado,
    shading = 'auto',
    cmap = cld_cmap
)

# colorbar
cbar = fig.colorbar(p, location = 'right', pad = 0.01)
cbar.ax.set_yticks(clevs)
cbar.ax.set_yticklabels(list(cloud_codes.values()))

# titulos e eixos
ax.set_xlabel("Latitude (°)")
ax.set_ylabel("Altitude (Km)")
plt.title('Classificacao de nuvens', loc='left')
ax.grid(True, color = "gray")

# salvando a figura
plt.savefig(os.path.join(output_, 'cloudsat_classificacao_continental.png'), bbox_inches='tight')
plt.close()
