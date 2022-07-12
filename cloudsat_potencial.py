#---IMPORTS----------------------------------------------------------------------------------------
import os
import scipy.ndimage

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
# geoprof_fname = '2019055170406_68325_CS_2B-GEOPROF_GRANULE_P_R04_E08_F01.h5'
geoprof_fname = '2010225163215_22839_CS_2B-GEOPROF_GRANULE_P1_R05_E03_F00.h5'

# ecmwf_fname = '2019055170406_68325_CS_ECMWF-AUX_GRANULE_P_R05_E08_F03.h5'
ecmwf_fname = '2010225163215_22839_CS_ECMWF-AUX_GRANULE_P_R05_E03_F00.h5'

# recorte da area de estudo
# lat_min = -35
# lat_max = -27
# lon_min = -65
# lon_max = -40
lat_min = -40.6
lat_max = -29.1
lon_min = -47.6
lon_max = -44.4
extent = [lon_min, lon_max, lat_min, lat_max] # South America

#---Cloudsat Refletividade--------------------------------------------------------------------------

# Refletividade do cloudsat
cldst_radar = read_data(os.path.join(input, geoprof_fname))

# dimensoes do dado
cldst_lons, cldst_lats, cldst_height, cldst_time, cldst_elev = get_geodata(
    os.path.join(input, geoprof_fname), return_list=True
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

#---Cloudsat Temperatura Potencial--------------------------------------------------------------------------

# Temperatura e Pressao retirado do ECMWF e interpolado na trajetoria do cloudsat
ecmwf_temp = read_data(os.path.join(input, ecmwf_fname), data_field = 'Temperature') # em kelvin
ecmwf_press = read_data(os.path.join(input, ecmwf_fname), data_field = 'Pressure')*1e-2 # converte de Pa para hPa
ecmwf_specific_umidity =  read_data(os.path.join(input, ecmwf_fname), data_field = 'Specific_humidity') # em kg/kg

# calculo de theta
R = 0.286
theta = ecmwf_temp*(1000/ecmwf_press)**(R)

# calculo da pressao de vapor de saturacao
A = 2.53*(10**8)*10 # kPam converte pra hPa
B = 5.42*(10**3) # K

# calculo da theta_e (temperatura potencial equivalente)
Td = B/np.log(A*0.622/(ecmwf_press*ecmwf_specific_umidity)) # em kelvin
Tncl = 1/(1/(Td-56)+np.log(ecmwf_temp/Td)/800)+56
L = 2.5e6
Cp = 1005
theta_e = theta*np.exp(ecmwf_specific_umidity*L/(Cp*Tncl))

# demais dimensoes do dado do ecmwf interpolado no cloudsat
ecmwf_lons, ecmwf_lats, ecmwf_height, ecmwf_time, ecmwf_elev = get_geodata(
    os.path.join(input, ecmwf_fname), return_list=True,
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
theta_e = theta_e[j1:j2, :]

# interpolar a temperatura potencial para os niveis de referencia
theta_e = cloudsat_utils._interp2d_ecmwf(
    theta_e.filled(np.nan),
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
              'rightside_up': True, 'use_clabeltext': True}
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

# nromaliza com base nos intervalos
radr_norm = mpl.colors.BoundaryNorm(clevs, radr_cmap.N)
radr_kw = dict(cmap=radr_cmap, norm=radr_norm)

# plot contourf
p = ax.contourf(
    cldst_lats[i1:i2],
    np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz),
    cldst_radar,
    levels = clevs,
    **radr_kw
)
cbar = fig.colorbar(p, location = 'right', pad = 0.03)

# plot da elevacao
ax.fill_between(
    cldst_lats[i1:i2],
    cldst_elev[i1:i2],
    color = "black"
)

# titulos e eixos
ax.set_xlabel("Latitude (°)")
ax.set_ylabel("Altitude (Km)")
plt.title('Refletividade do Radar (DbZ) e Temperatura Potencial Equivalente (K)', loc='left')

# salvando a figura
plt.savefig(os.path.join(output, 'cloudsat_refletividade.png'), bbox_inches='tight')
plt.close()
