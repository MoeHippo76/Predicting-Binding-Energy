#!/usr/bin/python3
"""
Created on Sun Apr 17 12:32:29 2022

@author: Ilham
"""

import numpy as np
from  matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import subprocess
import os
import shutil

#Correlation Matrix            
R = np.array([[  1.000000000000, -0.346456089960, -0.143602277420,  0.315457763048, -0.253113631056,  0.057551986284, -0.142555053005,  0.049793515737, -0.317267358672, -0.333077764936 ],
         [ -0.346456089960,  1.000000000000,  0.706703100972,  0.233695858004, -0.247166569223, -0.150693828263, -0.201000871763, -0.173751865385, -0.224620446870, -0.177159589117 ],
         [ -0.143602277420,  0.706703100972,  1.000000000000,  0.363377172265, -0.353336963669, -0.773888570041, -0.323747888208, -0.132193258832, -0.362528502960, -0.287229425362 ],
         [  0.315457763048,  0.233695858004,  0.363377172265,  1.000000000000, -0.991265140658, -0.216249061452, -0.855139894226, -0.096484542058, -0.992455543341, -0.971189044814 ],
         [ -0.253113631056, -0.247166569223, -0.353336963669, -0.991265140658,  1.000000000000,  0.190295331761,  0.908399390584,  0.065173748885,  0.984774451952,  0.965533888467 ],
         [  0.057551986284, -0.150693828263, -0.773888570041, -0.216249061452,  0.190295331761,  1.000000000000,  0.215647587362,  0.205781326002,  0.224223670029,  0.154010161981 ],
         [ -0.142555053005, -0.201000871763, -0.323747888208, -0.855139894226,  0.908399390584,  0.215647587362,  1.000000000000, -0.026672742907,  0.853817429008,  0.837871899707 ],
         [  0.049793515737, -0.173751865385, -0.132193258832, -0.096484542058,  0.065173748885,  0.205781326002, -0.026672742907,  1.000000000000,  0.104472994861,  0.067120558586 ],
         [ -0.317267358672, -0.224620446870, -0.362528502960, -0.992455543341,  0.984774451952,  0.224223670029,  0.853817429008,  0.104472994862,  1.000000000000,  0.964057333441 ],
         [ -0.333077764936, -0.177159589117, -0.287229425361, -0.971189044814,  0.965533888467,  0.154010161981,  0.837871899707,  0.067120558586,  0.964057333441,  1.000000000000 ]])

#Mean Vector
global X_mean
X_mean = np.array([0.158706769332587,28.986789057772100, 40.004790480413600, 0.992423332283364, -45.135131022237300,-145.382167908057000,-186.065399575124000,-206.579593890106000, -74.026333176459900, -35.658261114791700]).T
#List of variances
global C
C = np.array([0.000416242011032, 0.603949405916, 13.13633652871, 0.122729138391, 5.36118191809, 52.1685613068, 5.04824419651, 23.1467313535, 18.51551028812, 13.04856045886])
#Getting the covariance matrix S
D = np.diag(C)
S = D @ R @ D 

eVa,eVe = np.linalg.eig(S) #eigen vector and eigen values 

L = np.diag(eVa) 
L_2 = np.sqrt(L)
T = eVe @ L_2 #Calculates Transformation to apply to normal scatter plot.

'''Genegerates parameter sets'''
def random_params():
    #Applying transformation
    X = np.random.normal(0,1,10).T
    Y = (T @ X) + X_mean
    return Y


'''
Following method writes the parameters into hfbth_v200d under UNE1
The first 6 parameters in the input list must be the last 6 of the event from the generator.
The last 4 parameters in the input list me be the first 4 of the event from the generator.
'''
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

'''
Following method logs all the parameter sets.
'''
def save_params(params,destination,index = None):
    fout = open(destination+'/Parameters', 'a')
    
    if index is None:
        string = ""
    else:
        string = str(index) + " "

    for p in params:
        string = string + " " + str(p)
    string = string + '\n'
    fout.write(string)
    fout.close()

'''
renames data file
'''
def rename(index,destination):
    old = destination + "/dft_mesh.dat"
    if index == 0:
        new = destination + "/dft_mesh_default.dat"
    else: 
        new = destination + "/dft_mesh_" + str(index) + ".dat" 
    os.rename(old,new)

'''
Reads all dft data files to plot
'''
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

'''
Plots histogram of minimas (Maybe works?)
'''
def plot_ground_state(directory,directory65):
    global N,element,e
    x_plot,y_plot = data_1_dft(True,directory+"/")
    ground_states = [min(y_plot)]
    ground_states_i = [np.where(y_plot == ground_states[-1])]
    for i in range(N):
        x_list,y_list = np.array(data_1_dft(False,directory65+"/"+"dft_mesh_"+str(i+1)+".dat"))
        ground_states.append(min(y_list))
        ground_states_i.append(np.where(y_list == ground_states[-1]))

    plt.hist(ground_states,N)
    plt.title("Ground state histogram of "+ element+" N_shells = " + e)
    plt.xlabel("Binding Energy (MeV)")
    plt.savefig(directory+element+"_"+e+" ground state")
    plt.show()

'''
This box changes the hfbtho file and runs the box.
It then moves the file to desired location and renames them.
'''
def create_data(params,index,destination):
    global root,n,z,e
    #Writing parameters into 1_dft/hfbtho_v200
    write_parameters(params)
    #Running the boxes

    os.chdir(root + '/VSLAT/1_dft')    
    subprocess.call("nice -n 10 ./run.sh -n "+n+" -z "+z+" -e "+e+" -t true ",shell=True)
    
    #os.chdir(root + '/VSLAT/2_H')
    #subprocess.call("./run.sh -n "+n+" -z "+z+" -e "+e+" -t true",shell=True)

    #os.chdir(root+ '/VSLAT/3_chi')
    #subprocess.call("./run.sh",shell = True)

    #copying grid data into a seperate folder
    shutil.copy(root + "/VSLAT/1_dft/data/dft_mesh.dat",destination)
    #renaming the grid data.
    rename(index,destination)

'''
This method checks if all the parameters are within the 1,2,3 given standard deviations
'''
def eval_params_confidence(event,factor):
    flag = True
    mean = list(X_mean[4:]) + list(X_mean[:4])
    cov = list(C[4:]) + list(C[:4])
    for i,e in enumerate(event):
        if not (e <= mean[i] + (factor*cov[i]) and e >= mean[i] - (factor*cov[i])):
            flag = False
            break
    return flag

'''
This method reads parameters from an existing Parameter set file.
'''
def read_params(directory):
    lines = []
    with open(directory+"/Parameters", 'r') as fp:
        lines = fp.readlines()
    
    params = []
    for line in lines:
        data = line.split()
        if data == []:
            break
        p = [float(d) for d in data[1:]]
        params.append((int(data[0]) - 1,p))
    return params

'''
This method generates parameter sets. 
It returns a list of all the parameter sets with every parameter under 65% confidence band.
Logs all the parameter sets under this band.
'''
def generate_params(directory):
    global N,std_factor
    j = 0
    confidence_65 = []
    while len(confidence_65) < N:
        e = random_params()
        event = list(e[4:]) + list(e[:4])
        if eval_params_confidence(event,std_factor):
            save_params(event,directory,j+1)
            confidence_65.append((j,event))
            j += 1
    return confidence_65

'''
Saves Output files.
'''    
def save_all_data(index,destination):
    global root
    save_dest = destination + "/Output_files_"+str(index) 
    os.mkdir(save_dest)
    os.chdir(root + "/VSLAT/1_dft/")
    files = os.listdir()
    output = []
    for f in files:
        if "out" in f:
          output.append(f)
          shutil.copy(f,save_dest)

'''
Plots all the curves generated.
'''              
def plot_results(directory,directory65):
    global element,e
    
    os.chdir(directory65)
    files = os.listdir()
    for f in files:
        if not os.path.isdir(f):
            x_list,y_list = np.array(data_1_dft(False,directory65+'/'+f))
            plt.plot(x_list,y_list,color = 'red',label = "Within 1 standard deviation")
    
    x_plot,y_plot = data_1_dft(True,directory)
    plt.plot(x_plot,y_plot,color= 'black',label = "default")
        
    
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    plt.legend(handles, labels, loc='best')
    plt.xlabel("Beta_2")
    plt.ylabel("Binding Energy (MeV)")
    plt.title(element+ " " + "N_shells = " + e)
    plt.savefig(directory+element+"_"+e+".png")
    plt.show()

'''
This method plots slask.agr
'''
def plot_slask():
    z,x_1,x_2,y_1,y_2 = read_slask("/home/hammy/VSLAT/3_chi")
    plt.plot(x_1,y_1,label="fit")
    plt.plot(x_2,y_2,label = "dft")
    plt.legend()
    plt.show()

'''
Helper method used to read and return values from slask.agr
'''
def read_slask(path):
    with open(path+"/slask.agr",'r') as fp:
        lines = fp.readlines()
    
    beta_fit = []
    beta_dft = []
    fit = []
    flag = False
    dft = []
    for i,line in enumerate(lines):
        data = line.split()
        if data == []:
            flag = True
        if i == 0:
            int_strength = float(data[1])
        else:
            if (not flag):
                beta_fit.append(float(data[0]))
                fit.append(float(data[1]))
            elif data != []:
                beta_dft.append(float(data[0]))
                dft.append(float(data[1]))
    return int_strength,beta_fit,beta_dft,fit,dft

def dft(param,index):
    global root,n,z,e
    
    params = list(X_mean.T)
    params[index] = param
    l = []
    for i in range(10):
        l.append(params[i] - X_mean[i])

    params = params[4:] + params[:4]
    write_parameters(params)
    os.chdir(root + '/VSLAT/1_dft')    
    subprocess.call("nice -n 10 ./run.sh -n "+n+" -z "+z+" -e "+e+" -t true ",shell=True)

    y_list = data_1_dft(False,root + '/VSLAT/1_dft/data/dft_mesh.dat')[1]
    y = min(y_list)
    return y


def derivative(index):
    x = X_mean[index]
    del_x = C[index]/100
    if index in [0,1,2,3]:
        return 0,0,0
    y1 = dft(x+del_x,index)
    y2 = dft(x-del_x,index)
    del_y = (y1 - y2)/(2*del_x)
    return del_y,y1,y2

def linearization(directory):
    global element,e
    fl =  directory + "/log" + "_" +element+"_"+e+".txt" 
    std = 0
    with open(fl,"w") as fp:
        string = "Parameter\tderivative\tBE_plus\t\tBE_minus\n"
        fp.write(string)
        fp.close()
    
    derivatives  = []
    for i in range(10):
        with open(fl,"a") as fp:
            delta = derivative(i)
            derivatives.append(delta)
            string = str(X_mean[i]) + "\t" + str(delta[0]) + "\t" +  str(delta[1]) + "\t" + str(delta[2]) + "\n"
            fp.write(string)
            fp.close()
            
    with open(fl,"a") as fp:
        string = "\nstd\tx_i\tBE_plus\t\tBE_minus\tx_j\tBE_plus\t\tBE_minus\n" 
        fp.write(string)
        fp.close()        

    for i in range(10):
        if i in [0,1,2,3]:
            continue
        delta_xi,y1,y2 = derivatives[i]
        for j in range(10):
            if j in [0,1,2,3]:
                continue
            with open(fl,"a") as fp:
                delta_xj,y3,y4 = derivatives[j]
                std += R[i][j] * delta_xi * delta_xj
                string = str(std) + "\t" + str(X_mean[i])+ "\t" + str(y1) + "\t" + str(y2)+ "\t" + str(X_mean[j]) + "\t" + str(y3) + "\t" + str(y4)
                fp.write(string + "\n")
                fp.close()
    
    with open(fl,"a") as fp:
        fp.write("final std = " + str(std))
        fp.close()
    
if __name__ == "__main__":
    global root,n,z,e,element,N,std_factor

    '''
    Set n,z,e for the nucleus. Here e is n_max.
    Write the name of the element for title of directories and plot.
    Change root directory in the same format to specify where the data folder will be saved.
    N := number of parameter sets. Will be overwritten if reading an existing parameter set.
    Set run_boxes to False to not run any boxes in VSLAT.
    '''

    n = "92"
    z = "64"
    e = "8"
    element = "Gd156"    
    root = "/home/hammy"
    N = 20
    std_factor = 0.1
    
    run_VSLAT = False
    read_params_from_diff_dir = False
    visualize = False
    linearize = True

    
    direct_params = "" + "/Parameters"

    #default is saved in directory and varied data is saved in directory65
    directory65 = root +  "/dft_"+element+"_"+ e + "/data"
    directory = root + "/dft_"+element+"_"+ e
    
    
    #Directories created if needed.
    if not os.path.isdir(directory):
        os.mkdir(directory)
        os.mkdir(directory65)
    
    if linearize:
        lin = root+"/linearization"
        if not os.path.isdir(lin):
            os.mkdir(lin)
        linearization(lin)
    
    if read_params_from_diff_dir:
        shutil.copy(root + direct_params,directory)
        
    #Moving the first 4 parameters to the last for convenience of writing into hfbtho_v200d
    p = list(X_mean[4:]) + list(X_mean[:4])

        
    if run_VSLAT:
        #If Parameter set already exists it reads it.
        if os.path.isfile(directory+"/Parameters"):
            confidence_65 = read_params(directory)
            N = len(confidence_65)
        else:
            confidence_65 = generate_params(directory)

        #Running default parameter set.
        create_data(p,0,directory)
        save_all_data(0,directory)

        #Running varied parameter sets.
        for i in range(N):
            p = confidence_65[i][1]
            j = confidence_65[i][0]
            create_data(p,j+1,directory65)
            save_all_data(j+1,directory65)
    
    if visualize:     
        #Plotting results.
        plot_results(directory+"/",directory65)
        plot_ground_state(directory,directory65)
    
