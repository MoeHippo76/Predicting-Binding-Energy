import numpy as np
from  matplotlib import pyplot as plt
import subprocess
import os
import shutil
import parameter_generator as pg

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

def rename(index,destination):
    old = destination + "/slask.agr"
    if index == 0:
        new = destination + "/slask_default.agr"
    else: 
        new = destination + "/slask" + str(index) + ".agr" 
    os.rename(old,new)

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
    
    os.chdir(root + '/VSLAT/2_H')
    subprocess.call("./run.sh -n "+n+" -z "+z+" -e "+e+" -t true",shell=True)

    os.chdir(root+ '/VSLAT/3_chi')
    subprocess.call("./run.sh",shell = True)

    #copying grid data into a seperate folder
    shutil.copy(root + "/VSLAT/3_chi/slask.agr",destination)
    #renaming the grid data.
    rename(index,destination)

'''
This method checks if all the parameters are within the 1,2,3 given standard deviations
'''
def eval_params_confidence(event,factor):
    flag = True
    mean = list(pg.X_mean[4:]) + list(pg.X_mean[:4])
    cov = list(pg.C[4:]) + list(pg.C[:4])
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
This method generates parameter sets. 
It returns a list of all the parameter sets with every parameter under 65% confidence band.
Logs all the parameter sets under this band.
'''
def generate_params(directory):
    global N,std_factor
    j = 0
    confidence_65 = []
    while len(confidence_65) < N:
        e = pg.random_params()
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
    os.chdir(root + "/VSLAT/3_chi/")
    files = os.listdir()
    output = []
    for f in files:
        if ".png" in f:
          output.append(f)
          shutil.copy(f,save_dest)
    os.chdir(root+ "/VSLAT/1_dft")
    files = os.listdir()
    output = []
    for f in files:
        if ".out" in f:
          output.append(f)
          shutil.copy(f,save_dest)


def plot_all_slasks(directory):
    z_list = []
    for file in os.listdir(directory+"/data"):
        if file.endswith(".agr"):
            z,x_1,x_2,y_1,y_2 = read_slask(os.path.join(directory+"/data", file))
            plt.plot(x_1,y_1,label="fit within 1 std", color = 'hotpink')
            plt.plot(x_2,y_2,label = "dft within 1 std", color = 'red')
            z_list.append(z)

    z,x_1,x_2,y_1,y_2 = read_slask(directory+"/slask_default.agr")
    plt.plot(x_1,y_1,label="fit default", color = 'blue')
    plt.plot(x_2,y_2,label = "dft default", color = 'black')
    z_list.append(z)

    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    plt.legend(handles, labels, loc='best')
    plt.xlabel("Beta_2")
    plt.ylabel("Binding Energy (MeV)")
    plt.title(element+ " " + "N_shells = " + e)
    plt.savefig(directory+"/"+element+"_"+e)
    plt.show()
    
    plt.hist(z_list,N)
    plt.title("Histogram of "+ element+" N_shells = " + e + " int_strength")
    plt.xlabel("int_strength")
    plt.savefig(directory+"/"+element+"_"+e+"int_strength")
    plt.show()


'''
Helper method used to read and return values from slask.agr
'''
def read_slask(file):
    with open(file,'r') as fp:
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
    e = "10"
    element = "Gd156"    
    root = "/home/hammy"
    N = 20
    std_factor = 1
    
    run_VSLAT = False
    read_params_from_diff_dir = False
    visualize = True

    
    direct_params = "<target derectory>" + "/Parameters"

    #default is saved in directory and varied data is saved in directory65
    directory65 = os.getcwd() + "/box3_output/chi_"+element+"_"+ str(e) + "/data"
    directory = os.getcwd() + "/box3_output/chi_"+element+"_"+ str(e)
    
    
    #Directories created if needed.
    if not os.path.isdir(directory) and run_VSLAT:
        os.mkdir(directory)
        os.mkdir(directory65)

    
    if read_params_from_diff_dir:
        shutil.copy(root + direct_params,directory)
        
    #Moving the first 4 parameters to the last for convenience of writing into hfbtho_v200d
    p = list(pg.X_mean[4:]) + list(pg.X_mean[:4])

        
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
        plot_all_slasks(directory)
    