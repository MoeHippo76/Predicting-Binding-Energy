"""
Created on Sun Apr 17 12:32:29 2022

@author: Mulham
"""

from pprint import pprint
from unicodedata import decimal
import numpy as np
from  matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import subprocess
import os
import shutil
la = np.linalg

#Correlation Matrix
R = np.array([[    1,-0.35,-0.14, 0.32,-0.25,-0.06,-0.32,-0.33,-0.14, 0.05],
              [-0.35,    1, 0.71, 0.23,-0.25,-0.15,-0.22,-0.18,-0.20,-0.17],
              [-0.14, 0.71,    1, 0.36,-0.35,-0.77,-0.36,-0.29,-0.32,-0.13],
              [ 0.32, 0.23, 0.36,    1,-0.99,-0.22,-0.99,-0.97,-0.86,-0.10],
              [-0.25,-0.25,-0.35,-0.99,    1, 0.19, 0.98, 0.97, 0.91, 0.07],
              [-0.06,-0.15,-0.77,-0.22, 0.19,    1, 0.22, 0.15, 0.22, 0.21],
              [-0.32,-0.22,-0.36,-0.99, 0.98, 0.22,    1, 0.96, 0.85, 0.10],
              [-0.33,-0.18,-0.29,-0.97, 0.97, 0.15, 0.96,    1, 0.84, 0.07],
              [-0.14,-0.20,-0.32,-0.86, 0.91, 0.22, 0.85, 0.84,    1,-0.03],
              [ 0.05,-0.17,-0.13,-0.10, 0.07, 0.21, 0.10, 0.07,-0.03,    1]])

#Mean Vector
X_mean = np.array([0.15871,28.997, 40.005,0.992,-45.135,-145.382,
                   -186.065,-206.580,-74.026,-35.658]).T
#List of variances
C = np.array([0.00042,0.604,13.136,0.123,5.361,52.169,18.516,
              13.049,5.048,23.147])
#Getting the covariance matrix S
D = np.diag(C)
S = D @ R @ D 

eVa,eVe = la.eig(S) 

L = np.diag(eVa)
L_2 = np.sqrt(L)
T = eVe @ L_2

def random_params():
    X = np.random.normal(0,1,10).T

    Y = (T @ X) + X_mean
    return Y

def test_val(Y):
    t =  [(Y[0] >= 0.15 and Y[0] <= 0.17),
          (Y[1] >= 28 and Y[1] <= 36),
          (Y[2] >= 40 and Y[2] <= 100),
          (Y[3] >= 0.9 and Y[3] <= 1.5)]
    return t

def plot_distro(N,bins,Events,dimen = 1):
    for i in range(dimen):
        x = plot_distro_single(N,Events,i+1,bins)
        plt.hist(x,bins,label = str(dimen))
    plt.legend()
    plt.show()

def plot_distro_single(N,Events,dimen = 1,disp = False,bins = None):
    x = []
    for i in range(N):
        x.append(Events[i][dimen-1])

    if disp == True:
        plt.hist(x,bins,label = str(dimen))
        plt.legend()
        plt.show()

    return x

def plot_relation(par1,par2,all_draws,disp = False):
    P1 = []
    P2 = []
    for draw in all_draws:
        P1.append(draw[par1])
        P2.append(draw[par2])
    if disp == True:
        plt.scatter(P1,P2)
        plt.ylabel("Variable " + str(par2+1))
        plt.xlabel("Variable " + str(par1+1))   
        plt.show() 
    
    return P1,P2

def Monte_Carlo_events(N,bins = 0,vis = False,dimens = 10):
    List = []
    for i in range(N):
        Y = random_params()
        List.append(Y)   
    if vis == True:    
        plot_distro(N,bins,List,dimens)
    return List    


def Plot_correlation_graph():
    N = 100000
    bins = 100
    E = Monte_Carlo_events(N,bins,False,10)

    n = 10
    fig = plt.figure()
    gs = GridSpec(n, n)

    ax_scatters = np.array([])
    ax_hists = np.array([])
    for i in range(n):
        for j in range(n):
            if j > i:
                pass
            elif i == j:
                ax_hist_var = fig.add_subplot(gs[i,j])
                np.append(ax_hists,ax_hist_var)
                x = plot_distro_single(N,E,i+1)
                ax_hist_var.hist(x,bins)
            else:
                ax_scatter =  fig.add_subplot(gs[i, j])
                np.append(ax_scatters,ax_scatter)
                x,y = plot_relation(i,j,E)
                ax_scatter.scatter(x,y,s=0.001)
                

    plt.show()

def write_parameters(p):
     with open(r"/home/hammy/VSLAT/1_dft/src/hfbtho_v200d.f90", 'r') as fp:
        all_lines = fp.readlines()
        lines = all_lines[989:995]
        data = map(list,lines)
        for i,line in enumerate(data):
            string = ""
            for j,char in enumerate(line):
                string += char
                if char == '=':
                    new_line = string
                    all_lines[i+989] = new_line + " " + str(p[i])+"_pr\n"
                    break
        
        filehandle = open("/home/hammy/VSLAT/1_dft/src/hfbtho_v200d.f90", "w")
        filehandle.writelines(all_lines)
        filehandle.close()

def save_params(params,index = None):
    fout = open('/home/hammy/1_dft_8/Parameters', 'a')
    
    if index is None:
        string = ""
    else:
        string = str(index) + " "

    for p in params:
        string = string + " " + str(p)
    string = string + '\n'
    fout.write(string)
    fout.close()

def rename(index,destination):
    old = destination + "/dft_mesh.dat"
    new = destination + "/dft_mesh_" + str(index) + ".dat" 
    os.rename(old,new)

def data_1_dft(default = False,destination = None):

    lines = []
    
    if default:
        with open("/home/hammy/1_dft_8/dft_mesh_default.dat","r") as f:
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

def plot_ground_state():
    x_plot,y_plot = data_1_dft(default = True)
    ground_states = []
    for i in range(30):
        x_list,y_list = np.array(data_1_dft(False,i+1))
        ground_states.append(min(y_list))
    plt.hist(ground_states,10)
    plt.show()
        
    


def plot_bands(destination):
    x_plot,y_plot = data_1_dft(default = True)
    for i in range(30):
        x_list,y_list = np.array(data_1_dft(False,i+1,destination))
        plt.plot(x_list,y_list,color = 'red')
        
    
    ''' 

    all_y = np.array(all_y).T
    y_max = []
    y_min = []
    
    for y in all_y:
        y = np.sort(y)
        y_max.append(y[len(y)-2])
        y_min.append(y[1])

     
    plt.plot(x_plot,y_plot,'-k')
    #plt.fill_between(x_plot,y_min, y_max)
    plt.show()
    '''

def write_params_2_H(p):
    with open(r"/home/hammy/VSLAT/2_H/src/skyrme.f90", 'r') as fp:
        all_lines = fp.readlines()
        lines = all_lines[241:246]
        data = map(list,lines)
        for i,line in enumerate(data):
            for j,char in enumerate(line):
                if char == '=':
                    new_line = "".join(line[:(len(line)-j)])
            all_lines[i+989] = new_line + " " + str(p[i])+"_pr\n"
        
        filehandle = open("/home/hammy/VSLAT/1_dft/src/hfbtho_v200d.f90", "w")
        filehandle.writelines(all_lines)
        filehandle.close()



def create_data(params,index,destination):
    #Writing parameters into 1_dft/hfbtho_v200
    write_parameters(params)
    #Running the boxes
    os.chdir('/home/hammy/VSLAT/1_dft')
    subprocess.call("./run.sh -n 12 -z 12 -e 12 -t true",shell=True)
    #copying grid data into a seperate folder
    shutil.copy("/home/hammy/VSLAT/1_dft/data/dft_mesh.dat",destination)
    #renaming the grid data.
    rename(index,destination)

def evaluate_params():
    with open(r"/home/hammy/1_dft/Parameters","r") as f:
        lines = f.readlines()
    
    params = []
    for line in lines:
        data = line.split()
        data = [float(d) for d in data]
        params.append(np.array(data[1:]))
    all_events = np.array(params).T

    averages = list(map(np.mean,all_events))
    stds = list(map(np.std,all_events))
    
    counters = [0] * 6
    for i,event in enumerate(all_events):
        for e in event:
            if e <= X_mean[4+i] + (2*C[4+i]) and e >= X_mean[4+i] - (2*C[4+i]):
                counters[i] += 1
    confidence = [(c/32)*100 for c in counters]
    return confidence

def evaluate_event_generator():
    all_events = []
    N = 500
    for i in range(N):
        Event = np.array(random_params())
        all_events.append(Event)

    all_events = np.array(all_events).T
    averages = list(map(np.mean,all_events))
    stds = list(map(np.std,all_events))

    counters = [0] * 10

    for i,event in enumerate(all_events):
        for e in event:
            if e <= averages[i] + (2*stds[i]) and e >= averages[i] - (2*stds[i]):
                counters[i] += 1
           
    confidence = [(c/N)*100 for c in counters]
    print(confidence)

def eval_params_confidence(event,factor):
    flag = True
    for i,e in enumerate(event):
        if not (e <= X_mean[4+i] + (factor*C[4+i]) and e >= X_mean[4+i] - (factor*C[4+i])):
            flag = False
            break
    return flag

if __name__ == "__main__":
    
    directory65 = '/home/hammy/1_dft_12/conf_65'
    directory95 = '/home/hammy/1_dft_12/conf_95'
    directory = '/home/hammy/1_dft_12/'
    
    j = 0
    confidence_65 = []
    confidence_95 = []

    create_data(X_mean,0,directory)

    while j < 100:
        event = random_params()[4:]
        if eval_params_confidence(event,2):
            save_params(event,j+1)
            j += 1
            if eval_params_confidence(event,1):
                confidence_65.append((j,event))
            else:
                confidence_95.append((j,event))

    for i in range(2):
        p = confidence_65[i][1]
        j = confidence_65[i][0]
        print("65 percent",p,j)
        create_data(p,j,directory65)
        
    for i in range(2):
        p = confidence_95[i][1]
        j = confidence_95[i][0]
        print("95 percent",p,j)
        create_data(p,j,directory95)
    
    x_plot,y_plot = data_1_dft(default = True)
    plt.plot(x_plot,y_plot,color= 'black',label = "default")

    os.chdir(directory65)
    files = os.listdir()
    for f in files:
        x_list,y_list = np.array(data_1_dft(False,directory65+'/'+f))
        plt.plot(x_list,y_list,color = 'red',label = "65%")
    
    os.chdir(directory95)
    files = os.listdir()
    
    for f in files:
        x_list,y_list = np.array(data_1_dft(False,directory95+'/'+f))
        plt.plot(x_list,y_list,color = 'blue',label ="95%")
    
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    plt.legend(handles, labels, loc='best')
    plt.show()
    