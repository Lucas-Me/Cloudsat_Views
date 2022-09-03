# Script responsavel pelos plots de concentracao, raio efetivo e conteudo de agua na nuvem

#---IMPORTS----------------------------------------------------------------------------------------
import os

import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import matplotlib.colors as mcolors
import numpy as np  # for numerical array computations
import matplotlib.tri as tri

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

from cloudsat_read import get_geodata, read_data
import cloudsat_utils

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
# input_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Dados'
# output_ = '/mnt/f/Lucas/Conteudo/Fisica das nuvens e precipitacao/Figuras'
input_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
output_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'

'''Caso Frente Fria Continental'''
# # nome dos arquivos CWC-RO e do AUX-ECMWF (caso continental)
# cwc_fname = 'frente_continental_2B-CWC-RO_P1_R05.h5'
# ecmwf_fname = 'frente_continental_ECMWF-AUX_P_R05.h5'

# # recorte da area de estudo (caso continental)
# lat_min = -35
# lat_max = -27.5
# lon_min = -65
# lon_max = -40
# extent = [lon_min, lon_max, lat_min, lat_max] # South America

'''Caso Frente Fria Oceanica'''
# nome do arquivo geoprof 
cwc_fname = 'frente_oceanica_2B-CWC-RO.h5'
ecmwf_fname = 'frente_oceanica_ECMWF-AUX.h5'

# recorte da area de estudo
lat_min = -40
lat_max = -29
lon_min = -30
lon_max = -25
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#--------------------------------------------Cloudsat Temperatura --------------------------------------------------------------------------

# Temperatura retirado do ECMWF e interpolado na trajetoria do cloudsat
ecmwf_temp = read_data(os.path.join(input_, ecmwf_fname),
                        data_field = 'Temperature', fillmask= True)

# demais dimensoes do dado do ecmwf interpolado no cloudsat
ecmwf_lons, ecmwf_lats, ecmwf_height, ecmwf_time, ecmwf_elev = get_geodata(
    os.path.join(input_, ecmwf_fname), return_list = True,
    varnames = ['Longitude', 'Latitude', "EC_height", "Profile_time", "DEM_elevation"]
)

# Encontrar indices do array onde o recorte da area esta localizado
# diferente pro dado interpolado e pro extraido diretamente do cloudat (secao anterior)
ii = np.where(
    (ecmwf_lons > lon_min)
    & (ecmwf_lons < lon_max)
    & (ecmwf_lats > lat_min)
    & (ecmwf_lats < lat_max)
)[0]
i1, i2 = ii[0], ii[-1]

# grade na horizontal
ecmwf_x = np.arange(i1, i2, dtype = np.float32)

# grade na vertical
ecmwf_z0 = 0  # km
ecmwf_z1 = 16  # km
ecmwf_nz = 1000  # Number of pixels (levels) in the vertical.
ecmwf_z = (ecmwf_height * 1e-3).astype(np.float32)

# indexando a variavel temperatura potencial
temperature = ecmwf_temp [i1:i2, :]

# interpolar a temperatura potencial para os niveis de referencia
temperature = cloudsat_utils._interp2d_ecmwf(
    temperature,
    ecmwf_x,
    ecmwf_z,
    i1,
    i2,
    i2 - i1,
    ecmwf_z1,
    ecmwf_z0,
    ecmwf_nz,
).T[::-1, :]

#--------------------------------------------Cloudsat Conteudo de gelo----------

# variaveis da goticula, retirada do cloudsat
radius = read_data(os.path.join(input_, cwc_fname),'RO_ice_effective_radius',  fillmask= True) # micrometro
concentracao = read_data(os.path.join(input_, cwc_fname), 'RO_ice_number_conc',  fillmask= True) # cm^{-3}
conteudo = read_data(os.path.join(input_, cwc_fname), 'RO_ice_water_content',  fillmask= True) # mg/m^3

# demais dimensoes do dado do ecmwf interpolado no cloudsat
cwc_lons, cwc_lats, cwc_height, cwc_time, cwc_elev = get_geodata(
    os.path.join(input_, cwc_fname), return_list = True
)

# Encontrar indices do array onde o recorte da area esta localizado
# diferente pro dado interpolado e pro extraido diretamente do cloudat (secao anterior)
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
cwc_z0 = 0  # km
cwc_z1 = 16  # km
cwc_nz = 125  # Number of pixels (levels) in the vertical.
cwc_z = (cwc_height * 1e-3).astype(np.float32)[j1:j2, :]
zi = np.linspace(cwc_z0, cwc_z1, cwc_nz)

# coordenadas e valores
XX = np.tile(cwc_x, (125, 1)).T
x_coords = np.ravel(XX)
y_coords = np.ravel(cwc_z)
coords = np.column_stack((x_coords, y_coords))

# interpolar a temperatura potencial para os niveis de referencia
triang = tri.Triangulation(coords[:, 0], coords[:, 1])
Xi, Yi = np.meshgrid(cwc_x, zi)

# indexando as variaveis
radius = radius[j1:j2, :]
concentracao = concentracao[j1:j2, :]
conteudo = conteudo[j1:j2, :]*1e-3

vars_ = [radius, concentracao, conteudo]
for i in range(len(vars_)):
    interpolator = tri.LinearTriInterpolator(
        triang,
        np.ravel(vars_[i])
    )
    vars_[i] = interpolator(Xi, Yi)

#---Plotando os dados--------------------------------------------------------------------------

# Criando a figura e algumas constantes
fig, ax = plt.subplots(figsize=(22, 22), dpi = 200, nrows = 3, ncols = 1)
labels = [
    'Raio efetivo do gelo (um) e Temperatura (°C)',
    'Concentração de partículas de gelo (1/L) e Temperatura (°C)',
    'Conteúdo de gelo (g/m³) e Temperatura (°C)'
]

# intervalos de contorno das variaveis
clevs = [
    np.linspace(0, 180, 10), # radius, micrometros
    np.linspace(0, 400, 10), # concentracao, 1/L
    np.linspace(0, 2, 10) # conteudo, g/m^3
]

# lista de cores, em ordem crescete. RGBA
colors = np.array([
    'blue',
    'dodgerblue',
    'cyan',
    'palegreen',
    'gold', 
    'darkorange',
    'red'
])

# eixo horizontal
xticks, xticklabels = cloudsat_utils.xticks_alongtrack(
    x_indices = cwc_x,
    x_lons = cwc_lons[j1:j2],
    x_lats = cwc_lats[j1:j2],
    bins = 6
)

for row in range(3):
    # plot das isolinhas de temperatura
    kw_clabels = {'fontsize': 12, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                'rightside_up': True, 'use_clabeltext': True}
    temp_contour = ax[row].contour(
        ecmwf_lats[i1:i2],
        np.linspace(ecmwf_z0, ecmwf_z1, ecmwf_nz),
        temperature - 273.15,
        np.linspace(-60, 20, 5),
        colors='green',
        linewidths=1.25,
        linestyles='dashed'
    )
    ax[row].clabel(temp_contour, **kw_clabels) 

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
        np.take(cwc_lats, Xi.astype('int64')),
        Yi,
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
    #ax[row].set_xticks(xticks)
    #ax[row].set_xticklabels(xticklabels)

    # plot da elevacao
    ax[row].fill_between(
        cwc_lats[j1:j2],
        cwc_elev[j1:j2]*1e-3,
        color = "black"
    )

    # limite vertical do plot
    ax[row].set_ylim(bottom = 0)

# salvando a figura
plt.savefig(os.path.join(output_, 'oceanico_propriedades_gelo.png'), bbox_inches='tight')
plt.close()
