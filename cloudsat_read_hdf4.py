
# Script responsavel por ler o hdf do cloudsat, e extrair informacoes e dados

#---IMPORTS----------------------------------------------------------------------------------------
from pyhdf.SD import SD, SDC 
from pyhdf.HDF import *
from pyhdf.VS import *
import numpy as np

#---FUNCOES----------------------------------------------------------------------------------------

def get_hdf_geodata(hdfname : str, varnames : list):
    file = HDF(hdfname, SDC.READ) # leo hdf
    vs = file.vstart() # inicia a interface vdata

    data_info_list = vs.vdatainfo()
    values_list = []
    for var in varnames:
        if var == "Height":
            array = get_hdf_data(hdfname, var)
        else:
            vdata = vs.attach(var) # extrai o vdata da variavel
            array = vdata[:] # extrai os valores do array
            array = np.array(array).reshape((1, len(array)))[0] # padronizar o array
            
            vdata.detach() # "fecha" o vdata
    
        if var == "Profile_time":
            pass

        values_list.append(array)

    # termina a interface vdata e fecha o arquivo aberto
    vs.end()
    file.close()

    # retorna os arrays
    if len(values_list) == 1:
        values_list = values_list[0]

    return values_list


def get_hdf_data(hdfname : str, varname : str):
    file = SD(hdfname, SDC.READ) # le o arquivo

    sds_obj = file.select(varname) # seleciona o sds
    data = sds_obj.get().astype('float32') # extrai o array

    attrs = sds_obj.attributes() # atributos do objeto
    
    # checa se existe um atributo fill value
    if '_FillValue' in list(attrs.keys()):
        data[data == attrs['_FillValue']] = np.nan

    # eliminando valores negativos
    return data


    