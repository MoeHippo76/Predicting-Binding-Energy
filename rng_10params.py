# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 12:32:29 2022

@author: Mulham
"""

import numpy as np
from numpy import linalg as la
from  matplotlib import pyplot as plt

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
    
    
    #rerolling to make sure X[3] and X[2] are within allowed bounds.
    '''
    while(X[3]*C[3] + X_mean[3] < 0.9):
        X[3] = np.random.normal(0,1,1)[0]
    while(X[2]*C[2] + X_mean[2] < 40):
        X[2] = np.random.normal(0,1,1)[0] 
    '''
        
    
    Y = (T @ X) + X_mean
    return Y

def test_val(Y):
    t =  [(Y[0] >= 0.15 and Y[0] <= 0.17),
          (Y[1] >= 28 and Y[1] <= 36),
          (Y[2] >= 40 and Y[2] <= 100),
          (Y[3] >= 0.9 and Y[3] <= 1.5)]
    return t


if __name__ == "__main__":
    List = []
    for i in range(100):
         List.append(random_params())
    print(List)
        
        
    