#---IMPORTS----------------------------------------------------------------------------------------
import os

import matplotlib.pyplot as plt  # the plotting interface
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
cldclass_fname_1 = 'frente_continental_2B-CLDCLASS_P1_R05.h5'
cwc_fname_1 = 'frente_continental_2B-CWC-RO_P1_R05.h5'

# recorte da area de estudo
lat_min_1 = -35
lat_max_1 = -27.5
lon_min_1 = -65
lon_max_1 = -40

# nome dos arquivos CWC-RO e do AUX-ECMWF (caso oceanico)
cwc_fname_2 = 'frente_oceanica_2B-CWC-RO.h5'
cldclass_fname_2 = 'frente_oceanica_2B-CLDCLASS.h5'

# recorte da area de estudo
lat_min_2 = -40
lat_max_2 = -29
lon_min_2 = -30
lon_max_2 = -25

coords = [
    [lon_min_1, lon_max_1, lat_min_1, lat_max_1], # continental
    [lon_min_2, lon_max_2, lat_min_2, lat_max_2] # oceanica
]

#---Cloudsat Classificacao de Nuvens--------------------------------------------------------------------------
cldclass_casos = []
casos = [os.path.join(input_, cldclass_fname_1), os.path.join(input_, cldclass_fname_2)]

for i in range(2): # 2 casos
    cldclass, cloud_codes = read_cldclass(casos[i])

    # dimensoes do dado
    cldst_lons, cldst_lats, cldst_height, cldst_time, cldst_elev = get_geodata(
        os.path.join(input_, casos[i]), return_list= True
    )

    # Encontrar indices do array onde o recorte da area esta localizado
    ii = np.where(
    (cldst_lons > coords[i][0])
    & (cldst_lons < coords[i][1])
    & (cldst_lats > coords[i][2])
    & (cldst_lats < coords[i][3])
    )[0]
    i1, i2 = ii[0], ii[-1]

    # indexando a classificacao das nuvens
    cldclass = cldclass[i1:i2, :]
    cldclass_casos.append(cldclass)


    # # grade na vertical
    # cloudsat_z0 = 0  # km
    # cloudsat_z1 = 20  # km
    # cloudsat_nz = 100 # Number of pixels (levels) in the vertical.
    # cloudsat_z = (cldst_height[i1:i2, :] * 1e-3).astype(np.float32)
    # new_z = np.linspace(cloudsat_z0, cloudsat_z1, cloudsat_nz) # novo Z


    # # interpolador
    # interp = NearestNDInterpolator(list(zip(x_coords, z_coords)), values)
    # XX, ZZ = np.meshgrid(cloudsat_x, new_z)
    # cld_class_interpolado = interp(XX, ZZ)

#--------------------------------------------Cloudsat Conteudo de gelo----------

casos = [os.path.join(input_, cwc_fname_1), os.path.join(input_, cwc_fname_2)]
cwc_radius = []
cwc_concentracao = []
cwc_conteudo = []

for i in range(2):
    # variaveis da goticula, retirada do cloudsat
    radius = read_data(casos[i],'RO_ice_effective_radius',  fillmask= True) # micrometro
    concentracao = read_data(casos[i], 'RO_ice_number_conc',  fillmask= True) # L^{-1}
    conteudo = read_data(casos[i], 'RO_ice_water_content',  fillmask= True) # mg/m^3

    # demais dimensoes do dado do ecmwf interpolado no cloudsat
    cwc_lons, cwc_lats, cwc_height, cwc_time, cwc_elev = get_geodata(
        casos[i], return_list = True
    )

    # Encontrar indices do array onde o recorte da area esta localizado
    jj = np.where(
        (cwc_lons > coords[i][0])
        & (cwc_lons < coords[i][1])
        & (cwc_lats > coords[i][2])
        & (cwc_lats < coords[i][3])
    )[0]
    j1, j2 = jj[0], jj[-1]

    # # grade na horizontal
    # cwc_x = np.arange(j1, j2, dtype = np.int64)
    # nx = j2 - j1

    # # grade na vertical
    # cwc_z0 = 0  # km
    # cwc_z1 = 16  # km
    # cwc_nz = 125  # Number of pixels (levels) in the vertical.
    # cwc_z = (cwc_height * 1e-3).astype(np.float32)[j1:j2, :]
    # zi = np.linspace(cwc_z0, cwc_z1, cwc_nz)

    # # coordenadas e valores
    # XX = np.tile(cwc_x, (125, 1)).T
    # x_coords = np.ravel(XX)
    # y_coords = np.ravel(cwc_z)
    # coords = np.column_stack((x_coords, y_coords))

    # indexando as variaveis
    cwc_radius.append(radius[j1:j2, :])
    cwc_concentracao.append(concentracao[j1:j2, :])
    cwc_conteudo.append(conteudo[j1:j2, :]*1e-3)

#---Processamento  e Plot ---------------------------------------------------- #
nuvens = {'1': 'High clouds', '2': 'Altostratus', '3': 'Altocumulus', '4':'Stratus',
'5': 'Stratocumulus', '6': 'Cumulus', '7': 'Nimbostratus', '8': 'Deep Convective'}

for i in range(1, 9):
    nuvem = nuvens[str(i)]
    mask_cloud_cont = cldclass_casos[0].flatten() == i
    mask_cloud_oce = cldclass_casos[1].flatten() == i
    if np.any(mask_cloud_cont) and np.any(mask_cloud_oce):
        print(nuvem, "existe nos dois casos")
    else:
        continue
        
    # PLota os dados
    fig, ax = plt.subplots(ncols = 3, figsize = (15, 5), dpi = 200)
    masks = [mask_cloud_cont, mask_cloud_oce]
    colors = ['sienna', 'teal']
    origem = ["continental", "oceanica"]
    for j in range(2):
        radius_values = cwc_radius[j].flatten()[masks[j].flatten()]
        conc_values = cwc_concentracao[j].flatten()[masks[j].flatten()]
        conteudo_values = cwc_conteudo[j].flatten()[masks[j].flatten()]
        p = ax[0].scatter(radius_values, conc_values, marker = "o", edgecolor = 'none', facecolors = colors[j], alpha = 0.3, label = origem[j])
        ax[1].scatter(radius_values, conteudo_values, marker = "o", edgecolor = 'none', facecolors = colors[j], alpha = 0.3, label = origem[j])
        ax[2].scatter(conc_values, conteudo_values, marker = "o", edgecolor = 'none', facecolors = colors[j], alpha = 0.3, label = origem[j])
        #
        ax[0].set_xlabel("Raio efetivo do gelo (μm)")
        ax[0].set_ylabel('Concentração de partículas de gelo (1/L)')
        ax[1].set_xlabel("Raio efetivo do gelo (μm)")
        ax[1].set_ylabel('Conteúdo de gelo (g/m³)')
        ax[2].set_xlabel('Concentração de partículas de gelo (1/L)')
        ax[2].set_ylabel('Conteúdo de gelo (g/m³)')

    letra = ['a', 'b', 'c']
    for n in range(3):
        bottom, top = ax[n].get_ylim()
        left, right = ax[n].get_xlim()
        ax[n].set_xlim(0, right)
        ax[n].set_ylim(0, top)
        ax[n].text(right*0.05, top*0.95, f'({letra[n]})')
        ax[n].grid(True, axis = 'both')

    fig.text(0.45, 0.9, nuvem, fontsize = "x-large")
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol = 2, loc = (0.75, 0.91))
    fig.savefig(os.path.join(output_, f"Microfisica gelo {nuvem}.png"), bbox_inches='tight')
    plt.close()

