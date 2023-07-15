import numpy as np
from  matplotlib import pyplot as plt
import subprocess
import os
import shutil
import parameter_generator as pg


def write_parameters(p):
    global root
    strings = ["Cr","Cp","RHO_NM","ASS_NM"]

    with open(root + "/VSLAT/1_dft/src/hfbtho_v200d.f90", 'r') as fp:
        all_lines = fp.readlines()
        lines = all_lines[989:1011]
        
        data = map(list,lines)
        param_i = 0
        for i,line in enumerate(data):
            string = ""
            check = [st in lines[i] for st in strings]   
            if any(check):
                    for char in line:
                        string += char
                        if char == '=':
                            new_line = string
                            all_lines[i+989] = new_line + " " + str(p[param_i])+"_pr\n"
                            param_i += 1
                            break
        
        filehandle = open(root + "/VSLAT/1_dft/src/hfbtho_v200d.f90", "w")
        filehandle.writelines(all_lines)
        filehandle.close()

def data_1_dft(default = False,destination = None):

    lines = []
    
    if default:
        with open(destination+"dft_mesh_default.dat","r") as f:
            lines = f.readlines()
            f.close()
    else:
        with open(destination,"r") as f:
            lines = f.readlines()
            f.close()

    Beta = []
    Energy = []

    for line in lines[1:]:
        data = line.split()
        Beta.append(float(data[2]))
        Energy.append(float(data[6]))

    return Beta,Energy


def dft(param,index):
    global root,n,z,e
    
    params = list(pg.X_mean.T)
    params[index] = param
    l = []
    for i in range(10):
        l.append(params[i] - pg.X_mean[i])

    params = params[4:] + params[:4]
    write_parameters(params)
    os.chdir(root + '/VSLAT/1_dft')    
    subprocess.call("nice -n 10 ./run.sh -n "+n+" -z "+z+" -e "+e+" -t true ",shell=True)

    y_list = data_1_dft(False,root + '/VSLAT/1_dft/data/dft_mesh.dat')[1]
    y = min(y_list)
    return y


def derivative(index):
    global std_factor
    x = pg.X_mean[index]
    del_x = pg.del_xs[index]
    y1 = dft(x+del_x,index)
    y2 = dft(x-del_x,index)
    del_y = (y1 - y2)/(2*del_x)
    return del_y,y1,y2

def linearization(directory):
    global element,std_factor
    fl =  directory + "/log" + "_" +element+"_"+e+".txt"
    derivatives  = []
    std = 0
    if os.path.isfile(fl):
        with open(fl,"r") as fp:
            lines = fp.readlines()
            fp.close()
        for i in range(1,11):
            data = lines[i].split()
            data = (float(data[1]),float(data[2]),float(data[3]))
            derivatives.append(data)
    else:
        for i in range(10):
            with open(fl,"a") as fp:
                delta = derivative(i)
                derivatives.append(delta)
                fp.close()
                
    with open(fl,"w") as fp:
            string = "Parameter\tderivative\tBE_plus\t\tBE_minus\n"
            fp.write(string)
            for i,delta in enumerate(derivatives):
                string = str(pg.X_mean[i]) + "\t" + str(delta[0]) + "\t" +  str(delta[1]) + "\t" + str(delta[2]) + "\n"
                fp.write(string)
            fp.close()


    with open(fl,"a") as fp:
        string = "\nstd\tx_i\tBE_plus\t\tBE_minus\tx_j\tBE_plus\t\tBE_minus\n" 
        fp.write(string)
        fp.close()  

    for i in range(4,10):
        delta_xi,y1,y2 = derivatives[i]
        for j in range(4,10):
            with open(fl,"a") as fp:
                delta_xj,y3,y4 = derivatives[j]
                std += pg.R[i][j] * delta_xi * delta_xj * pg.C[i] * pg.C[j]
                string = str(std) + "\t" + str(pg.X_mean[i])+ "\t" + str(y1) + "\t" + str(y2)+ "\t" + str(pg.X_mean[j]) + "\t" + str(y3) + "\t" + str(y4)
                fp.write(string + "\n")
                fp.close()
    
    with open(fl,"a") as fp:
        fp.write("final std = " + str(std))
        print("final std = " + str(std))
        fp.close()

if __name__ == "__main__":
    global root,n,z,e,element,std_factor

    n = "12"
    z = "12"
    e = "4"
    element = "Mg24"    
    root = "/home/hammy"
    std_factor = 0.1
    directory = os.getcwd() + "/linearization_output" 
    linearization(directory)