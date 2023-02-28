import numpy as np

mass2symb = {
    '1.008':[1,"H"],
    '12.011':[6,"C"],
    '15.999':[8,"O"],
    '14.007':[7,"N"],
    '30.974':[15,"P"],
    '35.450':[17,"Cl"]
    }


def read_system_after_run(file):
    with open (file,"r") as systemFile:
        data_system = systemFile.readlines()
    data = []
    for line in data_system:
        temp = line.replace("\n", "")
        data.append(temp)
        
    return data
    


def system_Natoms(data):
    iline = 0
    for line in data:
        if "atoms" in line:
            N_atoms = int(line.split(" ")[0])
            break
    return N_atoms

def system_data_xyz(data):
    iline = 0
    for line in data:
        if "atoms" in line:
            N_atoms = int(line.split(" ")[0])
        
        elif "atom types" in line:
            N_atomTypes = int(line.split(" ")[0])
        
        elif "Masses" in line:
            N_types = iline + 2
        
        elif "Atoms # full" in line:
            ini_line = iline + 2
            break
        iline += 1
    
    #se crea diccionario #atomo+Masa Atóica
    atomType = {}
    for i in range(N_atomTypes):
        atomType[i+1] = data[N_types+i].split(" ")[1]
        
    system_data = []
    
    for i in range(ini_line,ini_line+N_atoms):
        line = data[i].split(" ")
        line_ = line[4:7]
        simbolo = mass2symb[atomType[int(line[2])]][1]
        line_.insert(0, simbolo)
        system_data.append(line_)
    return system_data



    
if __name__ == "__main__":
    data = read_system_after_run("system_after_min.data")
    
    Natoms = system_Natoms(data)
    
    system = system_data_xyz(data)
    with open('sistemaFinal.xyz', 'w') as file_sys:
        file_sys.write(f"{Natoms}\n")
        file_sys.write(f"Sistema despues de minimizacion\n")


        for line in system:
            file_sys.write(f"{line[0]}\t {float(line[1]):0.5f} \t {float(line[2]):0.5f} \t {float(line[3]):0.5f}\t\n")
        
        
    
    
    
    
    
    