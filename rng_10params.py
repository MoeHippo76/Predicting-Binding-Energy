"""
Created on Sun Apr 17 12:32:29 2022

@author: Mulham
"""

import numpy as np
from  matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
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
X_mean = np.array([0.15871,28.997, 40.005,0.992,-45.135,-148.382,
                   -186.065,206.580,-740.026,-35.658]).T
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

if __name__ == "__main__":
    Y = random_params()
    print(Y)
