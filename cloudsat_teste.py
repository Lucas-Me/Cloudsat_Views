# Script responsavel pelos plots de concentracao, raio efetivo e conteudo de agua na nuvem
# Manipulacao de dados do CLOUDSAT

#---IMPORTS----------------------------------------------------------------------------------------
import os

import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import matplotlib.colors as mcolors
import numpy as np  # for numerical array computations
import matplotlib.tri as tri


#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

from cloudsat_functions import get_hdf_geodata, get_hdf_data
import cloudsat_utils

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
input_ = '/mnt/f/lucas/conteudo/fisica das nuvens e precipitacao/Dados'
output = '/mnt/f/lucas/conteudo/fisica das nuvens e precipitacao/Figuras'

# nome do arquivo CWC-RO e do AUX-ECMWF
cwc_fname = '2007124170226_05408_CS_2B-CWC-RO_GRANULE_P1_R05_E02_F00.hdf'

# recorte da area de estudo
lat_min = -40.6
lat_max = -29.1
lon_min = -55.3
lon_max = -52
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#--------------------------------------------Cloudsat CWC-RO --------------------------------------------------------------------------

# Temperatura retirado do ECMWF e interpolado na trajetoria do cloudsat
radius = get_hdf_data(os.path.join(input_, cwc_fname), 'RO_liq_effective_radius') # micrometros
concentracao = get_hdf_data(os.path.join(input_, cwc_fname), 'RO_liq_number_conc') # cm^{-3}
conteudo = get_hdf_data(os.path.join(input_, cwc_fname), 'RO_liq_water_content') # mg/m^3

print(np.count_nonzero(~np.isnan(radius)))
print(np.count_nonzero(~np.isnan(concentracao)))
print(np.count_nonzero(~np.isnan(conteudo)))

# demais dimensoes do dado do CWC-RO, height esta em metros
cwc_lons, cwc_lats, cwc_height, cwc_time, cwc_elev = get_hdf_geodata(
        os.path.join(input_, cwc_fname), 
        ["Longitude", "Latitude", "Height", "Profile_time", "DEM_elevation"]
        )

# resgatando indices da area de estudo
jj = np.where(
    (cwc_lons > lon_min)
    & (cwc_lons < lon_max)
    & (cwc_lats > lat_min)
    & (cwc_lats < lat_max)
)[0]
j1, j2 = jj[0], jj[-1]

# grade na horizontal
cwc_x = np.arange(j1, j2, dtype = np.int64)
nx = j2 - j1

# grade na vertical
cwc_z0 = 0  # altura inicial [km]
cwc_z1 = 16  # altura final [km]
cwc_nz = 125  # Numero de pixels na vertical
cwc_z = (cwc_height * 1e-3)[j1:j2, :] # converte de m para km
zj = np.linspace(cwc_z0, cwc_z1, cwc_nz)

# coordenadas e valores
XX = np.tile(cwc_x, (125, 1)).T
x_coords = np.ravel(XX)
y_coords = np.ravel(cwc_z)
coords = np.column_stack((x_coords, y_coords))

# interpolar as variaveis para o nivel de referencia
triang = tri.Triangulation(coords[:, 0], coords[:, 1])
Xj, Yj = np.meshgrid(cwc_x, zj)

# indexando as variaveis
radius = radius[j1:j2, :]
concentracao = concentracao[j1:j2, :]
conteudo = conteudo[j1:j2, :]*1e-3 # conversao para g*m^{-3}

vars_ = [radius, concentracao, conteudo]
for i in range(len(vars_)):
    # eliminando valores negativos
    print(np.count_nonzero(~np.isnan(vars_[i])))
    vars_[i][vars_[i] < 0] = np.nan

    interpolator = tri.LinearTriInterpolator(
        triang,
        np.ravel(vars_[i])
    )
    vars_[i] = interpolator(Xj, Yj)

#---Plotando os dados--------------------------------------------------------------------------

# Criando a figura e algumas constantes
fig, ax = plt.subplots(figsize=(22, 22), dpi = 200, nrows = 3, ncols = 1)
labels = [
    'Raio efetivo da água líquida (um) e Temperatura (°C)',
    'Concentração de partículas de água líquida (1/cm³) e Temperatura (°C)',
    'Conteúdo de água líquida (g/m³) e Temperatura (°C)'
]

# intervalos de contorno das variaveis
clevs = [
    np.linspace(0, 180, 10), # radius, micrometros
    np.linspace(0, 400, 10), # concentracao, 1/cm³
    np.linspace(0, 2, 10) # conteudo, g/m^3
]

# lista de cores, em ordem crescete.
colors = np.array([
    'blue',
    'dodgerblue',
    'cyan',
    'palegreen',
    'gold', 
    'darkorange',
    'red'
])

for row in range(3):

    # cria um novo cmap a partir do pre-existente
    var_cmap = mcolors.LinearSegmentedColormap.from_list(
        'Custom cmap', colors, clevs[row].shape[0] - 1)
    var_cmap.set_bad('w')
    var_cmap.set_under('darkblue')
    var_cmap.set_over('darkred')

    # nromaliza com base nos intervalos
    var_norm = mpl.colors.BoundaryNorm(clevs[row], var_cmap.N, clip = False)
    var_kw = dict(cmap=var_cmap, norm=var_norm)

    # plot contourf
    p = ax[row].contourf(
        np.take(cwc_lats, Xj),
        Yj,
        vars_[row],
        levels = clevs[row],
        extend = 'max',
        **var_kw
    )
    cbar = fig.colorbar(p, location = 'right', pad = 0.02, ax = ax[row], extend = 'max')

    # titulos e eixos
    ax[row].set_xlabel("Latitude (°)")
    ax[row].set_ylabel("Altitude (Km)")
    ax[row].set_title(labels[row], loc='left')

    # plot da elevacao
    ax[row].fill_between(
        cwc_lats[j1:j2],
        cwc_elev[j1:j2]*1e-3,
        color = "black"
    )

    # limites do plot
    ax[row].set_ylim(bottom = 0)

# salvando a figura
plt.savefig(os.path.join(output, 'propriedades_agua_liquida.png'), bbox_inches='tight')
plt.close()
