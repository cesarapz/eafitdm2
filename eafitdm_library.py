texto01="""write_once("In Init") {
#Sientase libre de modificar esto a su gusto!
units             real
atom_style        full
kspace_style 	    pppm 1.0e-4
pair_style        lj/cut/coul/long  10.0	10.0
bond_style        harmonic
angle_style       harmonic
dihedral_style    opls
improper_style    cvff
pair_modify       mix arithmetic
boundary          p p p
}

"""
##########################################################################################
#Diccionario para convertir el numero atomico en el simbolo atómico
SymZ = {1:'H',2:'He',3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',10:'Ne',11:'Na',
        12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',20:'Ca',35:'Br',53:'I'}
##########################################################################################
#Diccionario para convertir la masa atomica en el simbolo atómico o en el numero atomico
mass2id = {'1.008':['H' , 1 ],'4.003' :['He', 2 ],'6.941'  :['Li', 3 ],'9.012'  :['Be', 4 ],
          '10.810':['B' , 5 ],'12.011':['C' , 6 ],'14.007' :['N' , 7 ],'15.999' :['O' , 8 ],
          '18.998':['F' , 9 ],'20.180':['Ne',10 ],'22.990' :['Na',11 ],'24.305' :['Mg',12 ],
          '26.982':['Al',13 ],'28.086':['Si',14 ],'30.974' :['P' ,15 ],'32.066' :['S' ,16 ],
          '35.450':['Cl',17 ],'39.948':['Ar',18 ],'40.008' :['Ca',20 ],'79.904' :['Br',35 ],
          '126.904':['I' ,53 ]}
##########################################################################################
def file_lines_number(filename):
    """
    Funcion que determina el numero de lines de un archivo de texto
    entrada: filename
    salida:  numero de lineas
    Pendiente: programar excepcion.
    """
    with open(filename,"r") as file:
        return len(file.readlines())

def concatenate_files(file1,file2):
    """
    This function append the content of file1 in file2.
    """
    f1=open(file1,'r')
    f2=open(file2,'a')
    f1Content = f1.read().splitlines()

    for line in f1Content:
        f2.write(f"{line}\n")
    f1.close()
    f2.close()
    print ("Concatenacion realizada con exito!")
    return True

def create_packmol_input(systemType,compos,dimensions):
    """
    Funcion que crea el input de packmol, de acuerdo con la opcion de nombre de sistema se
    activa una plantilla.
    """
    pmFile = open("inputPackmol.inp","w")
    pmFile.write("#Archivo de entrada de packmol escrito durante la ejecucion de eafitdm\n")
    pmFile.write("#Codigo diseñado por Jorge David Caro y Cesar Augusto Perez\n")
    pmFile.write("#Codigo escrito por Cesar Augusto Perez, Julio de 2022\n")
    pmFile.write("tolerance	1.0\n")
    pmFile.write("seed -1\n")
    pmFile.write("nloop0 1000\n")
    pmFile.write("filetype xyz\n")
    pmFile.write("output system.xyz\n")
    pmFile.write(f"\n#Sistema tipo {systemType}\n")

    #segun el tipo de sistema se ejecuta una seccion de codigo con el template característico de dicho sistema
    if systemType == 'mip':
        components = list(compos.keys())
        L = dimensions[0]
        lmin = L/2 - L/50
        lmax = L/2 + L/50
        Central_Atom = dimensions[1]
        i = 0
        for icompo in components:

            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {compos[icompo][2]}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{L:.2f}\n")
            if i == 0:      #la sintaxis del archivo de entrada en el modo mip acompaña la dimension del sisteam con la etiqueta del atomo cetral de la platilla
                pmFile.write(f"    atoms {Central_Atom}\n")
                pmFile.write(f"         inside box {lmin:.2f}\t{lmin:.2f}\t{lmin:.2f}\t{lmax:.2f}\t{lmax:.2f}\t{lmax:.2f}\n")
                pmFile.write("     end atoms\n")
            pmFile.write("end structure\n\n")
            i += 1
        inputInter = [L,L,L]
        return inputInter

    elif  systemType == 'solucion':
        components = list(compos.keys())
        for icompo in components:
            L = dimensions[0]
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {compos[icompo][2]}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{L:.2f}\n")
            pmFile.write("end structure\n\n")
        inputInter = [L,L,L]
        return inputInter

    elif  systemType == 'solucion2':
        components = list(compos.keys())
        for icompo in components:
            L = dimensions[0]
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {compos[icompo][2]}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{L:.2f}\n")
            pmFile.write("end structure\n\n")
        inputInter = [L,L,L]
        return inputInter



    elif  systemType == 'micela':
        components = list(compos.keys())
        L = dimensions[0]
        R = dimensions[1]
        for icompo in components[:-1]:
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {compos[icompo][2]}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{L:.2f}\n")
            pmFile.write(f"    atoms  {compos[icompo][4]}\n")
            pmFile.write(f"       inside sphere\t{L/2:.2f}\t{L/2:.2f}\t{L/2:.2f}\t{R:.2f}\n ")
            pmFile.write("    end atoms\n")
            pmFile.write(f"    atoms  {compos[icompo][5]}\n")
            pmFile.write(f"       outside sphere\t{L/2:.2f}\t{L/2:.2f}\t{L/2:.2f}\t{R:.2f}\n ")
            pmFile.write("    end atoms\n")
            pmFile.write("end structure\n")

        solvent = components[-1]
        pmFile.write(f"\nstructure {compos[solvent][3]}.xyz\n")
        pmFile.write(f"    number {compos[solvent][2]}\n")
        pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{L:.2f}\n")
        pmFile.write("end structure\n")
        inputInter = [L,L,L]
        return inputInter

    #En proceso de rellenar estas plantillas.
    elif  systemType == 'bicapaM':   #bicapa fundida (melt)
        components = list(compos.keys())
        H = dimensions[0]
        L = dimensions[1]
        for icompo in components[:-1]:
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {compos[icompo][2]}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t{H/4:.2f}\t{L:.2f}\t{L:.2f}\t{3*H/4:.2f}\n")
            pmFile.write("end structure\n")
        
        
        solvent = components[-1]
        pmFile.write(f"\nstructure {compos[solvent][3]}.xyz\n")
        pmFile.write(f"    number {compos[solvent][2]}\n")
        pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{H:.2f}\n")
        pmFile.write("end structure\n")
        inputInter = [L,L,H]
        pmFile.write(f"structure {compos[icompo][3]}\n")
        return inputInter
    
    elif  systemType == 'bicapa':   #bicapa fundida (melt)
        components = list(compos.keys())
        H  = dimensions[0]
        L  = dimensions[1]
        
                
        for icompo in components[:-1]:
            
            numberMolecules = int(compos[icompo][2])
            Ap = compos[icompo][4]                #Atomo limite entre la cabeza-cola
            lm = float(compos[icompo][5])          #longitud molecular
            #UP
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {int(numberMolecules/2)}\n")
            
            pmFile.write(f"    inside box 0.00\t0.00\t{H/2:.2f}\t{L:.2f}\t{L:.2f}\t{H/2+lm:.2f}\n")
            
            pmFile.write(f"    atoms {Ap}\n")
            pmFile.write(f"       over plane 0.00\t 0.00\t 1.00\t {H/2 + 0.95*lm:.2f}\n")
            pmFile.write(f"    end atoms\n")
            pmFile.write(f"end structure\n")
            #DOWN
            
            pmFile.write(f"structure {compos[icompo][3]}.xyz\n")
            pmFile.write(f"    number {int(numberMolecules/2)}\n")
            pmFile.write(f"    inside box 0.00\t0.00\t{H/2-lm:.2f}\t{L:.2f}\t{L:.2f}\t{H/2:.2f}\n")
            pmFile.write(f"    atoms {Ap}\n")
            pmFile.write(f"       below plane 0.00\t 0.00\t 1.00\t {H/2-0.95*lm:.2f}\n")
            pmFile.write(f"    end atoms\n")
            pmFile.write(f"end structure\n")
        
        solvent = components[-1]
        
        pmFile.write(f"\nstructure {compos[solvent][3]}.xyz\n")
        pmFile.write(f"    number {int(int(compos[solvent][2])/2)}\n")
        pmFile.write(f"    inside box 0.00\t0.00\t0.00\t{L:.2f}\t{L:.2f}\t{H/4:.2f}\n")
        pmFile.write("end structure\n")
        
        pmFile.write(f"\nstructure {compos[solvent][3]}.xyz\n")
        pmFile.write(f"    number {int(int(compos[solvent][2])/2)}\n")
        pmFile.write(f"    inside box 0.00\t0.00\t{3*H/4:.2f} \t{L:.2f}\t{L:.2f}\t{H:.2f}\n")
        pmFile.write("end structure\n")
        
        inputInter = [L,L,H]
        #pmFile.write(f"structure {compos[icompo][3]}\n")
        return inputInter
    
    
    
    elif systemType == 'surfactante':
        pmFile.write(f"structure {compos[icompo][3]}\n")
        return True
    elif systemType == 'interface':
        pmFile.write(f"structure {compos[icompo][3]}\n")
        return True
    elif systemType == 'interfaceM':
        pmFile.write(f"structure {compos[icompo][3]}\n")
        return True
    elif systemType == 'surfactante':
        pmFile.write(f"structure {compos[icompo][3]}\n")
        return True
    else:
        return False
    pmFile.close()

def create_moltemplate_input(component_props,systemSize):
    """
    Codigo python que crea el system.lt con una composicion específica. La primera parte.
    """
    with open('system.lt',"w") as syslt:
        syslt.write(f"{60*'#'}\n")
        syslt.write("#Archivo de entrada para moltemplate creado con eafitdm\n")
        for component in component_props.keys():
            #print(component)
            syslt.write(f"import\t'{component}.lt'\n")
        syslt.write(f"{texto01}\n")
        syslt.write(f"write_once('"'Data Boundary'"') {\n")
        syslt.write(f"\t0.0\t{systemSize[0]:.2f}\txlo\t xhi\n")
        syslt.write(f"\t0.0\t{systemSize[1]:.2f}\tylo\t yhi\n")
        syslt.write(f"\t0.0\t{systemSize[2]:.2f}\tzlo\t zhi\n")
        syslt.write("}\n\n")

        i = 1
        for component in component_props.keys():
            syslt.write(f"molecs{i} = new {component}[{component_props[component][2]}]\n")
            i += 1
        syslt.write(f"write_once('"'In Settings'"') {\n")

##############################################################################
#el siguiente script lee un archivo lmp generado en LIGPARGEN y carga la información de los componentes
#del sistema.
def read_lmp (file):
    """
    Esta función lee el *.lmp de una molecula de trabajo y crea un diccionario con esta informacion:
        claves: atomic, bond, angle, dihedral, improper
        valores: topologia y parametros.
    """
    #Cargar el contenido de lmp como un lista, linea por elemento.
    with open (file,"r") as lmpfile:
        content = lmpfile.readlines()

    #Se inicia el diccionario donde se cargarán los datos lmp de la molecula
    #esta estructura se retorna.
    molecule_parameters = {}

    #lectura de parametros de la molecula
    N_atoms =     int(content[2].split()[0])
    N_bonds =     int(content[3].split()[0])
    N_angles =    int(content[4].split()[0])
    N_dihedrals = int(content[5].split()[0])
    N_impropers = int(content[6].split()[0])
    #captura de marcadores de bloques de datos
    net_content = [line.replace('\n','') for line in content if line != '\n']
    content = net_content
    net_content = []
    for line in content:
        text = line.replace(' ','')
        net_content.append(text)

    #marcadores de parametros: Numeros de linea donde comieza cada bloque de informacion
    Mark_Atoms = net_content.index('Atoms') + 1
    Mark_Mass  = net_content.index('Masses') + 1
    Mark_PairC = net_content.index('PairCoeffs') + 1
    Mark_BondC = net_content.index('BondCoeffs') + 1
    Mark_Bonds = net_content.index('Bonds') + 1

    #Se lee información de masa, epsilon, sigma, carga y coordenadas
    atomic_parameters  = []
    j = Mark_PairC
    k = Mark_Atoms
    ii = 0
    jj = 1
    listSym = []

    for i in range(Mark_Mass,Mark_Mass+N_atoms):
        element = mass2id[content[i].split()[1]]
        temp0 = [element[1]]
        temp1 = [float(content[i].split()[1])]
        temp2 = [float(x) for x in content[j].split()[1:]]
        temp3 = [float(x) for x in content[k].split()[1:]]
        temp4 = [f"{element[0]}{jj}"]
        listSym.append(f"{element[0]}{jj}")
        atomic_parameters.append(temp0+temp1+temp2+temp3+temp4)
        j  += 1
        k  += 1
        ii += 1
        jj += 1
    molecule_parameters['atomic']=atomic_parameters

    #cargar la informacion de enlaces
    bond_parameters = []
    k = Mark_Bonds
    for i in range(Mark_BondC,Mark_BondC+N_bonds):
        temp1 = [float(x) for x in content[i].split()[1:]]
        temp2 = [int(x) for x in content[k].split()]
        Sym1 = listSym[temp2[2]-1]
        Sym2 = listSym[temp2[3]-1]
        temp2[2],temp2[3] = Sym1, Sym2

        bond_parameters.append(temp1+temp2)
        k += 1
    molecule_parameters['bond']=bond_parameters

    #cargar la informacion de enlaces angulos de enlace, sii hay
    if 'AngleCoeffs' in net_content:
        Mark_AnglC = net_content.index('AngleCoeffs') + 1
        Mark_Angle = net_content.index('Angles') + 1
        angle_parameters = []
        k = Mark_Angle
        for i in range(Mark_AnglC,Mark_AnglC+N_angles):
            temp1 = [float(x) for x in content[i].split()[1:]]
            temp2 = [int(x) for x in content[k].split()]
            #reemplazar numero de atomos por el simbolo#
            Sym1 = listSym[temp2[2]-1]
            Sym2 = listSym[temp2[3]-1]
            Sym3 = listSym[temp2[4]-1]
            temp2[2],temp2[3],temp2[4] = Sym1, Sym2, Sym3

            angle_parameters.append(temp1+temp2)
            k += 1
        molecule_parameters['angle']=angle_parameters

    #cargar la informacion de enlaces angulos dihedros, sii hay.

    if 'DihedralCoeffs' in net_content:
        Mark_DiheC = net_content.index('DihedralCoeffs') + 1
        Mark_Dihed = net_content.index('Dihedrals') + 1
        dihedral_parameters = []
        k = Mark_Dihed
        for i in range(Mark_DiheC,Mark_DiheC+N_dihedrals):
            temp1 = [float(x) for x in content[i].split()[1:]]
            temp2 = [int(x) for x in content[k].split()]
            #reemplazar numero de atomos por el simbolo#
            Sym1 = listSym[temp2[2]-1]
            Sym2 = listSym[temp2[3]-1]
            Sym3 = listSym[temp2[4]-1]
            Sym4 = listSym[temp2[5]-1]
            temp2[2],temp2[3],temp2[4],temp2[5] = Sym1, Sym2, Sym3, Sym4

            dihedral_parameters.append(temp1+temp2)
            k += 1
        molecule_parameters['dihedral']=dihedral_parameters

    #cargar la informacion de enlaces angulos dihedros, sii hay.
    if 'ImproperCoeffs' in net_content:
        Mark_ImprC = net_content.index('ImproperCoeffs') + 1
        Mark_Impro = net_content.index('Impropers') + 1
        improper_parameters = []
        k = Mark_Impro
        for i in range(Mark_ImprC,Mark_ImprC+N_impropers):
            temp1 = [float(x) for x in content[i].split()[1:]]
            temp2 = [int(x) for x in content[k].split()]
            #reemplazar numero de atomos por el simbolo#
            Sym1 = listSym[temp2[2]-1]
            Sym2 = listSym[temp2[3]-1]
            Sym3 = listSym[temp2[4]-1]
            Sym4 = listSym[temp2[5]-1]
            temp2[2],temp2[3],temp2[4],temp2[5] = Sym1, Sym2, Sym3, Sym4
            improper_parameters.append(temp1+temp2)
            k += 1
        molecule_parameters['improper']=improper_parameters
    return molecule_parameters


def create_lt(molecule, compDic):
    """
    Esta funcion carga el diccionario de datos lmp de la molecula y crea el *.lt  y el *.xyz
    de la molecula
    """
    tt = '{' #Buscar otra forma, esto es muy grosero.
    with open(f"{molecule}.lt","w") as ltfile:
        ltfile.write("#################################################################################\n")
        ltfile.write("#Archivo en formato moltemplate, traducido a partir de la salida tipo lammps\n")
        ltfile.write("#producida por la aplicacion en linea LiParGen\n")
        ltfile.write("#################################################################################\n")
        ltfile.write(f"{molecule} {tt} \n")

        #Escritura de atomos
        ltfile.write(f"     write('"'Data Atoms'"') { \n")
        i = 1
        for atom in compDic['atomic']:
            nAt = atom[0]
            Q   = atom[6]
            X   = atom[7]
            Y   = atom[8]
            Z   = atom[9]
            Sy  = atom[10]
            i  += 1
            ltfile.write(f"\t\t$atom:{Sy}\t$mol:.\t@atom:{Sy}\t{Q:.5f}\t{X:.5f}\t{Y:.5f}\t{Z:.5f}\t\n")
        ltfile.write("     }\n\n")

        #Escritura de masas
        ltfile.write(f"     write_once('"'Data Masses'"') { \n")
        i = 1
        for atom in compDic['atomic']:
            nAt = atom[0]
            M   = atom[1]
            Sy  = atom[10]
            ltfile.write(f"        @atom:{Sy}    {M:.5f}\n")
            i  += 1
        ltfile.write("     }\n")

        #Escritura de DESCRIPTORES DE ENLACES
        ltfile.write(f"     write('"'Data Bonds'"') { \n")
        for bond in compDic['bond']:
            Sy1 =bond[4]
            Sy2 =bond[5]
            ltfile.write(f"\t\t$bond:{Sy1}{Sy2}\t@bond:{Sy1}{Sy2}\t$atom:{Sy1}\t$atom:{Sy2}\n")
        ltfile.write("     }\n")

        #Escritura de PARAMETROS DE ENLACES
        ltfile.write(f"     write_once('"'In Settings'"') { \n")
        for bond in compDic['bond']:
            Sy1 =bond[4]
            Sy2 =bond[5]
            ltfile.write(f"\t\tbond_coeff\t@bond:{Sy1}{Sy2}\t\t{bond[0]:.5f}\t{bond[1]:.5f} \n")
        ltfile.write("     }\n")

        if 'angle' in compDic.keys():
        #Escritura de DESCRIPTORES DE ANGULO
            ltfile.write(f"     write('"'Data Angles'"') { \n")

            for angle in compDic['angle']:
                Sy1 =angle[4]
                Sy2 =angle[5]
                Sy3 =angle[6]
                ltfile.write(f"\t\t$angle:{Sy1}{Sy2}{Sy3}\t\t@angle:{Sy1}{Sy2}{Sy3}  \t$atom:{Sy1}\t$atom:{Sy2}\t$atom:{Sy3}\n")
            ltfile.write("     }\n")
            #Escritura de PARAMETROS DE ANGULO
            ltfile.write(f"     write_once('"'In Settings'"') { \n")

            for angle in compDic['angle']:
                Sy1 =angle[4]
                Sy2 =angle[5]
                Sy3 =angle[6]
                ltfile.write(f"\t\tangle_coeff\t@angle:{Sy1}{Sy2}{Sy3}\t\t{angle[0]:.5f}\t{angle[1]:.5f}\n")
            ltfile.write("     }\n")

        if 'dihedral' in compDic.keys():
        #Escritura de DESCRIPTORES DE ANGULO DIHEDRO
            ltfile.write(f"     write('"'Data Dihedrals'"') { \n")

            for dihedral in compDic['dihedral']:
                Sy1 =dihedral[6]
                Sy2 =dihedral[7]
                Sy3 =dihedral[8]
                Sy4 =dihedral[9]
                ltfile.write(f"\t\t$dihedral:{Sy1}{Sy2}{Sy3}{Sy4}\t@dihedral:{Sy1}{Sy2}{Sy3}{Sy4}\t$atom:{Sy1}\t$atom:{Sy2}\t$atom:{Sy3}\t$atom:{Sy4}\n")
            ltfile.write("     }\n")

            #Escritura de PARAMETROS DE ANGULO DIHEDRO
            ltfile.write(f"     write_once('"'In Settings'"') { \n")

            for dihedral in compDic['dihedral']:
                Sy1 =dihedral[6]
                Sy2 =dihedral[7]
                Sy3 =dihedral[8]
                Sy4 =dihedral[9]

                ltfile.write(f"\t\tdihedral_coeff\t@dihedral:{Sy1}{Sy2}{Sy3}{Sy4}\t\t{dihedral[0]:.5f}\t{dihedral[1]:.5f}\t{dihedral[2]:.5f}\t{dihedral[3]:.5f}\n")
            ltfile.write("     }\n")

        #Escritura de DESCRIPTORES DE ANGULO DIHEDRO IMPROPIO
            ltfile.write(f"     write('"'Data Impropers'"') { \n")

            for improper in compDic['improper']:
                Sy1 =improper[5]
                Sy2 =improper[6]
                Sy3 =improper[7]
                Sy4 =improper[8]
                ltfile.write(f"\t\t$improper:{Sy1}{Sy2}{Sy3}{Sy4}\t@improper:{Sy1}{Sy2}{Sy3}{Sy4}\t$atom:{Sy1}\t$atom:{Sy2}\t$atom:{Sy3}\t$atom:{Sy4}\n")
            ltfile.write("     }\n")

            #Escritura de PARAMETROS DE ANGULO DIHEDRO IMPROPIO
            ltfile.write(f"     write_once('"'In Settings'"') { \n")

            for improper in compDic['improper']:
                Sy1 =improper[5]
                Sy2 =improper[6]
                Sy3 =improper[7]
                Sy4 =improper[8]

                ltfile.write(f"\t\timproper_coeff\t@improper:{Sy1}{Sy2}{Sy3}{Sy4}\t\t{improper[0]:.5f}\t{int(improper[1])}\t{int(improper[2])}\n")
            ltfile.write("     }\n")
        ltfile.write("}\n")

def create_xyz(molecule, compDic):
    with open(f"{molecule}.xyz","w") as xyzfile:
        xyzfile.write(f" {len(compDic['atomic'])} \n")
        xyzfile.write(f" {molecule} \n")
        for atom in compDic['atomic']:
            xyzfile.write(f"{atom[10]}\t{atom[7]:2.5f}\t{atom[8]:2.5f}\t{atom[9]:2.5f}\n")

    with open(f"{molecule}_.xyz","w") as xyzfile:
        xyzfile.write(f" {len(compDic['atomic'])} \n")
        xyzfile.write(f" {molecule} \n")
        for atom in compDic['atomic']:
            sym = ''.join([i for i in atom[10] if not i.isdigit()])
            xyzfile.write(f"{sym}\t{atom[7]:2.5f}\t{atom[8]:2.5f}\t{atom[9]:2.5f}\n")

















