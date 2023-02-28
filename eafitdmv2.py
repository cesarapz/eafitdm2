import os, sys
import shutil
from eafitdm_library import *
import numpy as np
import math

def eafitdm(inputFile,outdir,moleculeDB):
    #eafitdm version python es la evolucion de eafitdm, creado incialmente  en korn shell.
    #eafitdm es un script que permite la creacion del conjunto de archivos inciales
    #para ejecutar una simulacion por dinamica molecular.
    #
    #Ideado por Cesar Augusto Perez  y Jorge David Caro.
    #programador: Cesar Augusto Perez
    #Version python iniciada en marzo de 2022.
    #######################################################################################
    #LECTURA

    #MOLEDIR = '/home/cesar/bin/MoleculeDataBase/'
    MOLEDIR = moleculeDB
    MOLEDIR = 'MoleculeDataBase/'
    list_files = []
    

    print(f"esta es la base de datos {moleculeDB}")
    
    
    
    #Lectura del input
    with open(inputFile, 'r') as f:
        lines = f.readlines()
        lines = [line.replace("\n","").replace("\t"," ") for line in lines]

    #leer tipo de sistema:
    systemType = lines[1].split(" ")[0]

    #leer numero de componentes
    nComponents = int(lines[2].split(" ")[0])

    #lectura de dimensiones
    dimensions = list(map(int,lines[3].split(" ")[0].split(",")))
    print ("Leyendo entrada:")
    print (f"Tipo de Sistema:{systemType}")
    print (f"#componentes:   {nComponents}")
    print (f"Dimensiones:    {dimensions}")

    #lectura de parametros de componentes:
    #Verificacion de si un componente no esta en la biblioteca, se cancela la ejecución.
    component_props = {}
    i = 0
    inputInter = []

    lista_componentes = []
    for icomp in range(4,4+nComponents):
        line = lines[icomp].split(",")
        name = line[1]
        print(f"buscando componente  {name} en la MoleculeDB")

        component_props[name] = line
        line.insert(0,i)

        if os.path.isfile(f"{MOLEDIR}{line[2]}.lmp"):

            print(f"La molecula {name} existe en la biblioteca de eafitdm")
            shutil.copyfile(f"{MOLEDIR}{line[2]}.lmp",f"{line[2]}.lmp")
            size = file_lines_number(f"{MOLEDIR}{line[2]}.lmp")
            line.insert(1,size)
            inputInter.append(f"{line[2]} {name}")
            lista_componentes.append(name)
        else:
            sys.exit(f"""No se encuentra referencia de La molécula <{name}> en la biblioteca de eafitdm
                se cancela la ejecucion, revise por favor el manual para alimentar la biblioteca""")
        i += 1

    #verificacion de existencia de archivo
    if os.path.exists("system.lt"):
        os.remove("system.lt")

    ###################################################################################
    #Creaccion de las coordenadas del sistema con packmol.

    list_epsilon = []    #lista de los valores de epsilon de cada atomo en el orden de aparicion.
    list_sigma = []      #lista de los valores de sigma de cada atomo en el orden de aparicion.
    list_pairs = []      #lista de las etiquetas para formar los pares en el system.lt

    for molecula in component_props.keys():

        #lee el archivo lammps del componente, descargado desde ligpargen
        mole = read_lmp(f"{molecula}.lmp")
        N_Atoms = len(mole['atomic'])
        i=1
        for iatom in range(N_Atoms):
            list_epsilon.append(mole['atomic'][iatom][2])
            list_sigma.append(mole['atomic'][iatom][3])
            Z = mole['atomic'][iatom][0]
            atomSymb = SymZ[Z]
            list_pairs.append(f"@atom:{molecula}/{atomSymb}{i}")
            i += 1

        #Crear archivos lt de las moleculas de interes
        create_lt(molecula,mole)
        #Crea archivos xyz de las moleculas de interes
        create_xyz(molecula,mole)

        list_files.append(f"{molecula}.lt")
        list_files.append(f"{molecula}.xyz")
        list_files.append(f"{molecula}_.xyz")
    #Llamado a creacion de archivo de entrada a packmol
    systemSize = create_packmol_input(systemType,component_props,dimensions)

    #Llamado a creacion de archivo system.lt (la primera parte)
    create_moltemplate_input(component_props,systemSize)

    #Complemento de archivo system.lt: mezcla de epsilon y sigma para los atomos
    N = len(list_epsilon)
    with open('system.lt',"a") as syslt:
        k = 0
        for i in range(N):
            text1 = list_pairs[i]
            for j in range(N):
                if j >= i:
                    text2 = list_pairs[j]
                    epsilon = math.sqrt(list_epsilon[i]*list_epsilon[j])
                    sigma   = math.sqrt(list_sigma[i]*list_sigma[j])
                    k += 1
                    syslt.write(f"\tpair_coeff\t{text1}\t{text2}\t\t{epsilon:.5f}\t{sigma:.5f}\n")
        syslt.write("}\n")
        list_files.append('system.lt')
        list_files.append('inputPackmol.inp')


    #crear directorio de ejecucion y copiar los archivos.

    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    else:
        os.mkdir(outdir)

    CWD = os.getcwd()
    for file in list_files:
        print(f'Moviendo el archivo {file}')
        shutil.move(f"{CWD}/{file}",f"{CWD}/{outdir}/{file}")


    shutil.copy(f"{MOLEDIR}/run.in.min",f"{CWD}/{outdir}/run.in.min")
    shutil.copy(f"{MOLEDIR}/run.in.npt",f"{CWD}/{outdir}/run.in.npt")



if __name__ == '__main__':

    ##lectura del input desde la linea de comandos:
    with open(sys.argv[1], 'r') as f:

        inputFile = sys.argv[1]

        #directorio donde se crearán los outputs
        Nreps = int(sys.argv[2])

        #directorio de la base de datos de moleculas
        moleculeDB = sys.argv[3]


    for i in range(Nreps):
        outdir = f'run0{i}'
        print(outdir,moleculeDB)
        eafitdm(inputFile,outdir,moleculeDB)

    os.system("rm *.lmp")











