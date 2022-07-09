#---IMPORTS----------------------------------------------------------------------------------------

import copy
import os
from pathlib import Path  # object-oriented file paths

import cartopy.crs as ccrs  # for plotting in cartographic projection
import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import numpy as np  # for numerical array computations

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

import cloudsat_utils
from cloudsat_read import get_geodata, read_data

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
input = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Trabalho/Dados'
output = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Trabalho/Figuras'

# nome do arquivo geoprof
geoprof_fname = '2019055170406_68325_CS_2B-GEOPROF_GRANULE_P_R04_E08_F01.h5'

# recorte da area de estudo
lat_min = -35
lat_max = -27
lon_min = -65
lon_max = -40
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#---Cloudsat refletividade--------------------------------------------------------------------------

cloudsat_lons, cloudsat_lats, cloudsat_height, cloudsat_time, elev = get_geodata(
    os.path.join(input, geoprof_fname), return_list=True
)
elev = elev * 1e-3 # convertendo elevacao para km.

cldst_radar = read_data(os.path.join(input, geoprof_fname))

# Encontrar indices do array onde o recorte da area esta localizado
ii = np.where(
    (cloudsat_lons > lon_min)
    & (cloudsat_lons < lon_max)
    & (cloudsat_lats > lat_min)
    & (cloudsat_lats < lat_max)
)[0]
i1, i2 = ii[0], ii[-1]

# grade na horizontal
cloudsat_x = np.arange(i1, i2, dtype = np.float32)

# grade na vertical
cloudsat_z0 = 0  # km
cloudsat_z1 = 18  # km
cloudsat_nz = 1000  # Number of pixels (levels) in the vertical.
#cloudsat_z = np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz)
cloudsat_z = (cloudsat_height * 0.001).astype(np.float32)

# indexando a refletividade do radar
cldst_radar = cldst_radar[i1:i2, :]

# interpolar a refletividade do radar para niveis verticais regulares
cldst_radar = cloudsat_utils.cc_interp2d(
    cldst_radar.filled(np.nan),
    cloudsat_x,
    cloudsat_z,
    i1,
    i2,
    i2 - i1,
    cloudsat_z1,
    cloudsat_z0,
    cloudsat_nz,
).T[::-1, :]

#---Plotando os dados--------------------------------------------------------------------------

# configurando um colormap, e normalizando o mesmo
radr_cmap = copy.copy(mpl.cm.jet)
radr_cmap.set_bad("w")
radr_cmap.set_under("w")
radr_norm = mpl.colors.BoundaryNorm(np.linspace(-20, 30, 11), radr_cmap.N)
radr_kw = dict(cmap=radr_cmap, norm=radr_norm, rasterized=True, shading="auto")


# Plot da refletividade
fig, ax = plt.subplots(figsize=(20, 6))
p = ax.pcolormesh(
    cloudsat_lats[i1:i2],
    np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz),
    cldst_radar,
    **radr_kw
)
fig.colorbar(p)

# plot da elevacao
ax.fill_between(
    cloudsat_lats[i1:i2],
    elev[i1:i2],
    color = "black"
)

# titulos e eixos
ax.set_xlabel("Latitude (°)")
ax.set_ylabel("Altitude (Km)")
ax.set_title("Relfetividade do radar")

# salvando a figura
plt.savefig(os.path.join(output, 'cloudsat_refletividade.png'), bbox_inches='tight', pad_inches=0, dpi=100)
plt.close()