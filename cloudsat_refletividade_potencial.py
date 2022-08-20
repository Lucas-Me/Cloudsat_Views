#---IMPORTS----------------------------------------------------------------------------------------
import os

import matplotlib as mpl
import matplotlib.pyplot as plt  # the plotting interface
import matplotlib.colors as mcolors
import numpy as np  # for numerical array computations
import matplotlib.tri as tri
import metpy.calc.thermo as thermo # calculo da propriedades termodinamicas
from metpy.units import units

#---IMPORTS LOCAIS----------------------------------------------------------------------------------------

from cloudsat_read import get_geodata, read_data
import cloudsat_utils

#---VARIAVEIS E PREPARATIVOS--------------------------------------------------------------------------
# Diretorios de entrada e saida de arquivos
# input_ = '/mnt/f/lucas/conteudo/fisica das nuvens e precipitacao/Dados'
# output = '/mnt/f/lucas/conteudo/fisica das nuvens e precipitacao/Figuras'
input_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
output_ = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'

# nome do arquivo geoprof
geoprof_fname = 'frente_continental_2B-GEOPROF_P1_R05.h5'
ecmwf_fname = 'frente_continental_ECMWF-AUX_P_R05.h5'

# recorte da area de estudo
lat_min = -35
lat_max = -27.5
lon_min = -65
lon_max = -40
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#---Cloudsat Refletividade--------------------------------------------------------------------------

# Refletividade do cloudsat
cldst_radar = read_data(os.path.join(input_, geoprof_fname),
                        "Radar_Reflectivity", fillmask = True)

# dimensoes do dado
cldst_lons, cldst_lats, cldst_height, cldst_time, cldst_elev = get_geodata(
    os.path.join(input_, geoprof_fname), return_list=True,
    varnames = ["Longitude", "Latitude", "Height", "Profile_time", "DEM_elevation"]
)
cldst_elev = cldst_elev * 1e-3 # convertendo elevacao para km.

# Encontrar indices do array onde o recorte da area esta localizado
ii = np.where(
    (cldst_lons > lon_min)
    & (cldst_lons < lon_max)
    & (cldst_lats > lat_min)
    & (cldst_lats < lat_max)
)[0]
i1, i2 = ii[0], ii[-1]

# grade na horizontal
cloudsat_x = np.arange(i1, i2, dtype = np.float32)

# grade na vertical
cloudsat_z0 = 1  # km
cloudsat_z1 = 16  # km
cloudsat_nz = 1000  # Number of pixels (levels) in the vertical.
cloudsat_z = (cldst_height * 1e-3).astype(np.float32)
zi = np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz)

# coordenadas e valores
XX = np.tile(cloudsat_x, (cldst_radar.shape[1], 1)).T
x_coords = np.ravel(XX)
y_coords = np.ravel(cloudsat_z[i1:i2, :])
coords = np.column_stack((x_coords, y_coords))

# interpolar a temperatura potencial para os niveis de referencia
triang = tri.Triangulation(coords[:, 0], coords[:, 1])
Xi, Yi = np.meshgrid(cloudsat_x, zi)

# indexando a refletividade do radar
cldst_radar = cldst_radar[i1:i2, :]

# interpolando 
interpolator = tri.LinearTriInterpolator(
    triang,
    np.ravel(cldst_radar)
)
cldst_radar = interpolator(Xi, Yi)


#---Cloudsat Temperatura Potencial--------------------------------------------------------------------------

# Temperatura e Pressao retirado do ECMWF e interpolado na trajetoria do cloudsat
ecmwf_temp = read_data(os.path.join(input_, ecmwf_fname), 'Temperature', fillmask = True) # em kelvin
ecmwf_press = read_data(os.path.join(input_, ecmwf_fname), 'Pressure', fillmask = True)*1e-2 # converte de Pa para hPa
ecmwf_specific_umidity =  read_data(os.path.join(input_, ecmwf_fname), 'Specific_humidity', fillmask = True) # em kg/kg

ecmwf_dewpoint = thermo.dewpoint_from_specific_humidity(
    pressure = ecmwf_press * units.mbar,
    temperature = ecmwf_temp * units.kelvin,
    specific_humidity = ecmwf_specific_umidity * units("dimensionless")
).to(units.kelvin)

theta_e = thermo.equivalent_potential_temperature(
    pressure = ecmwf_press * units.mbar,
    temperature = ecmwf_temp * units.kelvin,
    dewpoint = ecmwf_dewpoint
)


# demais dimensoes do dado do ecmwf interpolado no cloudsat
ecmwf_lons, ecmwf_lats, ecmwf_height, ecmwf_time, ecmwf_elev = get_geodata(
    os.path.join(input_, ecmwf_fname), return_list= True,
    varnames = ["Longitude", "Latitude", "EC_height", "Profile_time", "DEM_elevation"]
)

# Encontrar indices do array onde o recorte da area esta localizado
# diferente pro dado interpolado e pro extraido diretamente do cloudat (secao anterior)
jj = np.where(
    (ecmwf_lons > lon_min)
    & (ecmwf_lons < lon_max)
    & (ecmwf_lats > lat_min)
    & (ecmwf_lats < lat_max)
)[0]
j1, j2 = jj[0], jj[-1]

# grade na horizontal
ecmwf_x = np.arange(j1, j2, dtype = np.float32)

# grade na vertical
ecmwf_z0 = 0  # km
ecmwf_z1 = 16  # km
ecmwf_nz = 1000  # Number of pixels (levels) in the vertical.
ecmwf_z = (ecmwf_height * 1e-3).astype(np.float32)

# indexando a variavel temperatura potencial
theta_e = np.array(theta_e[j1:j2, :])

# interpolar a temperatura potencial para os niveis de referencia
theta_e = cloudsat_utils._interp2d_ecmwf(
    theta_e,
    ecmwf_x,
    ecmwf_z,
    j1,
    j2,
    j2 - j1,
    ecmwf_z1,
    ecmwf_z0,
    ecmwf_nz,
).T[::-1, :]

# suavizando o resultando com o metodo de interpolacao Spline, fator zoom = 3

#---Plotando os dados--------------------------------------------------------------------------

# configurando um colormap, e normalizando o mesmo
fig, ax = plt.subplots(figsize=(22, 6), dpi = 200)

# plot das isolinhas de temperatura potencial
kw_clabels = {'fontsize': 12, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
              'rightside_up': True, 'use_clabeltext': True, 'colors' : "darkgreen"}
clevtheta = np.arange(250, 400, 5)
theta_contour = ax.contour(
    ecmwf_lats[j1:j2],
    np.linspace(ecmwf_z0, ecmwf_z1, ecmwf_nz),
    theta_e,
    clevtheta,
    colors='green',
    linewidths=1.25,
    linestyles='solid'
)
plt.clabel(theta_contour, **kw_clabels) 

# plot da refletividade, colorbar ->
clevs = np.linspace(-20, 30, 11)

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
radr_cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors/255, clevs.shape[0] - 1)
radr_cmap.set_bad("w")
radr_cmap.set_under("w")
radr_cmap.set_over("darkred")

# nromaliza com base nos intervalos
radr_norm = mpl.colors.BoundaryNorm(clevs, radr_cmap.N)
radr_kw = dict(cmap=radr_cmap, norm=radr_norm)

# plot contourf
p = ax.contourf(
    np.take(cldst_lats, Xi.astype(np.int64)),
    Yi,
    cldst_radar,
    levels = clevs,
    extend = "max",
    **radr_kw
)
cbar = fig.colorbar(p, location = 'right', pad = 0.03, extend = "max")

# plot da elevacao
ax.fill_between(
    cldst_lats[i1:i2],
    cldst_elev[i1:i2],
    color = "black"
)

# titulos e eixos
ax.set_xlabel("Latitude (°)")
ax.set_ylabel("Altitude (Km)")
ax.set_ylim(bottom = 0)
plt.title('Refletividade do Radar (dBZe) e Temperatura Potencial Equivalente (K)', loc='left')

# salvando a figura
plt.savefig(os.path.join(output_, 'cloudsat_refletividade_continental.png'), bbox_inches='tight')
plt.close()
