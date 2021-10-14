#!/usr/bin/python3
# coding: utf-8

# %% IMPORTS 
import numpy as np
import os
import sys
import pandas as pd

# %% CLASSES
def normal(n):
    a = '{:.9E}'.format(float(abs(n)))
    if n < 0:
        sign = "-"
    else:
        sign = " "
    e = a.find('E')
    return ' {}0.{}{}{}{:02d}'.format(sign,a[0],a[2:e],a[e:e+2],abs(int(a[e+1:])*1+1))

class MOLDEN:
    """
    This class handles MOLDEN files.
    """
    def __init__(self,file="MOLDEN",path=os.getcwd()):
        
        self.file = file
        self.path = path
        self.n_orb = 0 #number of orbitals
        #MOLDEN [ATOMS] section
        #         - pandas DataFrame {label , index, atom_number,coord}
        #         - contain [ATOMS] section
        #                 Name        |    Description        |   Type
        #             -   label       |   label of atom       |   str
        #             -   inex        |   index of atom       |   int
        #             -   atom_number |   number of atom      |   int
        #             -   coord       |   coordinates of atom |   np.array(3)
        self.atom_df =  pd.DataFrame({"label" : [], "index" : [], 
                                      "atom_number" : [], "coord" : []} )
        #MOLDEN [GTO] section
        self.gto_l = []
        
        #MOLDEN [MO] section
        #       - pandas DataFrame {Sym, Energy, Spin, Occup, Coeff}
        #       - contain [MO] section
        #               Name      |    Description        |   Type
        #             - Sym       |   symmetry of MO      |   str
        #             - Energy    |   energy of MO        |   float
        #             - Spin      |   alpha / beta        |   str
        #             - Occup     |   Occupitation of MO  |   float
        #             - Coeff     |   MO coefficients     |   np.array((n_orb,n_orb))
        self.mos_df = pd.DataFrame({"Sym" : [], "Energy" : [] , 
                                   "Spin" : [], "Occup" : [], "Coeff": [] })
        
        self.read() #read data from MOLDEN file
    
    def read(self):
        """
        Read molden file
        """
        path_fn = str(self.path) + '/' + str(self.file)    
            
        atom_data = {"label" : [], "index" : [], "atom_number" : [], "coord" : []}  
        mo_data = {"Sym" : [], "Energy" : [] , "Spin" : [], "Occup" : [], "Coeff": [] }
        
        first = 0
        
        with open(path_fn) as file:
            at_copy = False
            gt_copy = False
            mo_copy = False
            coeff_l = [] #np.zeros((num_orb,1))
            
            for line in file:
                if line.strip().lower() == "[atoms] au":
                    at_copy = True
                    continue
                elif line.strip().lower() == "[molden format]":
                    at_copy = False
                    continue
                elif line.strip().lower() == "[gto]":
                    at_copy = False
                    gt_copy = True
                    continue
                elif line.strip().lower() == "[mo]":
                    gt_copy = False
                    mo_copy = True
                    continue
                
                #Read [ATOMS] section
                if at_copy:
                    atom_data["label"].append(line.split()[0])
                    atom_data["index"].append(int(line.split()[1]))
                    atom_data["atom_number"].append(line.split()[2])
                    atom_data["coord"].append(np.array(list(np.float_(line.split()[3:6]))))
                #Read [GTO] section
                elif gt_copy:
                    self.gto_l.append(line)
                    if  line.split() != []:
                        if line.split()[0] == "s":
                            self.n_orb += 1
                        elif line.split()[0] == "p":
                            self.n_orb += 3
                        elif line.split()[0] == "d":
                            self.n_orb += 6
                        elif line.split()[0] == "f":
                            self.n_orb += 10
                        elif line.split()[0] == "g":
                            self.n_orb += 15     
                #Read [MO] section  
                elif mo_copy: 
                    if "=" in line:
                        ind = line.split("=")[0].replace(" ","").lower()
                        item = line.split("=")[1]
                        if ind == "sym":
                            mo_data["Sym"].append(item.split()[0])
                            i = 0
                        elif ind == "ene":
                            mo_data["Energy"].append(float(item))
                        elif ind == "spin":
                            mo_data["Spin"].append(item.split()[0]) 
                        elif ind == "occup":
                            mo_data["Occup"].append(float(item))
                    else:
                        tmp = int(line.split()[0])
                        coeff_l.append(float(line.split()[1]))
                        i += 1
                        if i == self.n_orb:
                            coeff_arr = np.array(coeff_l)
                            mo_data["Coeff"].append(coeff_arr)
                            coeff_l = []
        
        self.atom_df = self.atom_df.append(pd.DataFrame(atom_data))
        self.atom_df = self.atom_df.set_index("index")
        self.mos_df = self.mos_df.append(pd.DataFrame(mo_data))
        self.mos_df.index += 1                               #start index at 1
    
    def write(self,unit="AU"):
        """
        Write molden file
        
        Input
            atom_df     - pandas DataFrame {label, index, atom_number,coord}
                        - contain [ATOMS] section from MOLDEN file
            gto_l       - list 
                        - contain [GTO] section from MOLDEN file
            mos_df      - pandas DataFrame {Sym, Energy, Spin, Occup, Coeff}
                        - contain [MO] section from MODLEN file
        """
        if len(str(self.file).split(".")) > 1:
            path_fn = str(self.path) + '/' + str(self.file.split(".")[0]) + "_vmd." +\
                '.'.join([str(elem) for elem in self.file.split(".")[1:]])
        else:
            path_fn = str(self.path) + '/' + str(self.file) + "_vmd.mld"
            
        #if os.path.exists(path_fn): 
        #     os.rename(self.path + '/' + self.file,
        #               self.path + '/old_' + self.file)
        
        with open(path_fn,'w') as file:
            file.write("[Molden Format]\n")
            
            #Write ATOM section
            file.write("[Atoms] " + str(unit) + "\n")
            
            for atom_i, row in self.atom_df.iterrows():
                file.write('{0:>3}{1:>6.0f}{2:>3.0f}{3:>14.6f}{4:>14.6f}{5:>14.6f}\n'.format(
                    row["label"],atom_i,float(row["atom_number"]),row["coord"][0],row["coord"][1],row["coord"][2]))
     
            #file.write("[Molden Format]\n")
            
            #write GTO section
            file.write("[GTO]\n")
            for line in self.gto_l:
                #if line.replace(" ","") != "\n":
                if (len(line.split()) == 2 and
                 float(line.split()[1]) == 0 and
                 line.split()[1].find("E")==-1 and 
                 line.split()[1].find("D")==-1 and
                 line.split()[1].find("e")==-1 and
                 line.split()[1].find("d")==-1):
                    file.write("{0:>3}{1:>2}\n".format(line.split()[0],line.split()[1]))
                else:
                  if len(line.split()) == 2:
                    #print("second",line.replace("\n",""))
                    #print(normal(float(line.split()[0])),normal(float(line.split()[1])))
                    exps = str(normal(float(line.split()[0]))) + str(normal(float(line.split()[1]))) + "\n"
                    file.write(exps)
                #file.write(" ")
                if len(line.split()) == 3:
                    #print("last",line.replace("\n",""))
                    file.write("{0:>2}{1:>5.0f}{2:>5.2f}\n".format(
                        line.split()[0],float(line.split()[1]),float(line.split()[2])))
                if len(line.split()) < 2:
                    file.write("\n")
            #file.writelines(self.gto_l)
            
            #Write MO section
            file.write("\n")
            file.write("\n")
            file.write("[MO]\n")
            #mos_arr = MOs_matrix(mos_df)
            for index, row in self.mos_df.iterrows():
                #file.write('{0: >2}{1: >2}\n'.format("Sym=",row["Sym"]))
                file.write('{0: >2}  {1: 3.4f}\n'.format("Ene=",row["Energy"]))
                file.write('{0: >3} {1: >5}\n'.format("Spin=",row["Spin"]))
                file.write('{0: >4}{1: 1.1f}\n'.format("Occup=",abs(round(row["Occup"],0))))
                mo_sub_i = 1
                for coeff in row["Coeff"]:
                    file.write('{0: >4}{1:<2}{2: 3.6f}\n'.format(mo_sub_i,"",coeff))
                    mo_sub_i += 1
#%%Main
if __name__ == "__main__":
    fn = sys.argv[1]
    mol = MOLDEN(file=fn,path=os.getcwd())
    mol.write()
    