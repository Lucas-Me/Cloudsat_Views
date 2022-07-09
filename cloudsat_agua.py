# Script responsavel pelos plots de concentracao, raio efetivo e conteudo de agua na nuvem

#---IMPORTS----------------------------------------------------------------------------------------
import os

import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import matplotlib.colors as mcolors
import numpy as np  # for numerical array computations

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

import cloudsat_utils
from cloudsat_read import get_geodata, read_data

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
input = './Dados'
output = './Figuras'

# nome do arquivo geoprof
cwc_fname = '2019055170406_68325_CS_2B-CWC-RO_GRANULE_P1_R05_E08_F03.h5'
ecmwf_fname = '2019055170406_68325_CS_ECMWF-AUX_GRANULE_P_R05_E08_F03.h5'

# recorte da area de estudo
lat_min = -35
lat_max = -27
lon_min = -65
lon_max = -40
extent = [lon_min, lon_max, lat_min, lat_max] # South America


#--------------------------------------------Cloudsat Temperatura --------------------------------------------------------------------------

# Temperatura retirado do ECMWF e interpolado na trajetoria do cloudsat
ecmwf_temp = read_data(os.path.join(input, ecmwf_fname), data_field = 'Temperature')

# demais dimensoes do dado do ecmwf interpolado no cloudsat
ecmwf_lons, ecmwf_lats, ecmwf_height, ecmwf_time, ecmwf_elev = get_geodata(
    os.path.join(input, ecmwf_fname), return_list=True,
    varnames = ["Longitude", "Latitude", "EC_height", "Profile_time", "DEM_elevation"]
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
    temperature.filled(np.nan),
    ecmwf_x,
    ecmwf_z,
    i1,
    i2,
    i2 - i1,
    ecmwf_z1,
    ecmwf_z0,
    ecmwf_nz,
).T[::-1, :]

#--------------------------------------------Cloudsat Conteudo de agua liquda--------------------------------------------------------------------------

# variaveis da goticula, retirada do cloudsat
radius = read_data(os.path.join(input, cwc_fname), data_field = 'RO_liq_effective_radius') # micrometro
concentracao = read_data(os.path.join(input, cwc_fname), data_field = 'RO_liq_number_conc') # cm^{-3}
conteudo = read_data(os.path.join(input, cwc_fname), data_field = 'RO_liq_water_content') # mg/m^3

# demais dimensoes do dado do ecmwf interpolado no cloudsat
cwc_lons, cwc_lats, cwc_height, cwc_time, cwc_elev = get_geodata(
    os.path.join(input, cwc_fname), return_list=True,
    varnames = ["Longitude", "Latitude", "Height", "Profile_time", "DEM_elevation"]
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
cwc_x = np.arange(j1, j2, dtype = np.float32)

# grade na vertical
cwc_z0 = 0  # km
cwc_z1 = 16  # km
cwc_nz = 1000  # Number of pixels (levels) in the vertical.
cwc_z = (cwc_height * 1e-3).astype(np.float32)

# indexando as variaveis
radius = radius[j1:j2, :]
concentracao = concentracao[j1:j2, :]
conteudo = conteudo[j1:j2, :]*1e-3


# interpolar a temperatura potencial para os niveis de referencia
vars_ = [radius, concentracao, conteudo]
for i in range(len(vars_)):
    vars_[i] = cloudsat_utils._interp2d(
        vars_[i].filled(np.nan),
        cwc_x,
        cwc_z,
        j1,
        j2,
        j2 - j1,
        cwc_z1,
        cwc_z0,
        cwc_nz,
    ).T[::-1, :]
    
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
    np.linspace(0, 40, 10), # radius, micrometros
    np.linspace(0, 100, 10), # concentracao, cm^{-3}
    np.linspace(0, 4, 10) # conteudo, g/m^3
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

    # nromaliza com base nos intervalos
    var_norm = mpl.colors.BoundaryNorm(clevs[row], var_cmap.N)
    var_kw = dict(cmap=var_cmap, norm=var_norm)

    # plot contourf
    p = ax[row].contourf(
        cwc_lats[j1:j2],
        np.linspace(cwc_z0, cwc_z1, cwc_nz),
        vars_[row],
        levels = clevs[row],
        **var_kw
    )
    cbar = fig.colorbar(p, location = 'right', pad = 0.02, ax = ax[row])

    # titulos e eixos
    ax[row].set_xlabel("Latitude (°)")
    ax[row].set_ylabel("Altitude (Km)")
    ax[row].set_title(labels[row], loc='left')

# salvando a figura
plt.savefig(os.path.join(output, 'propriedades_agua.png'), bbox_inches='tight')
plt.close()
