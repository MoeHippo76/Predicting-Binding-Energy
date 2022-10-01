import numpy as np
from  matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

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

'''Generates and returns a list of Events. Also calls for plotting distribution'''
def Monte_Carlo_events(N,bins = 0,vis = False,dimens = 10):
    List = []
    ev = Event_Visualizer()
    for i in range(N):
        Y = random_params()
        List.append(Y)   
    if vis == True:    
        ev.plot_distro(N,bins,List,dimens)
    return List 

'''
Evaluates event generator.
Keeps track of how many events are within two standard deviations.
Prints the percentage.
'''

def evaluate_event_generator():
    all_events = []
    N = 30
    for i in range(N):
        Event = np.array(random_params())
        all_events.append(Event)

    all_events = np.array(all_events).T
    averages = list(map(np.mean,all_events))
    stds = list(map(np.std,all_events))
    corr = np.corrcoef(all_events)

    corr_error = []
    for i in range(10):
        b = []
        for j in range(10):
            b.append(round(abs(R[i][j] - corr[i][j]),5))
        corr_error.append(b)

    error = [abs(a[0] - a[1]) for a in zip(stds,C)]

    counters = [0] * 10

    for i,event in enumerate(all_events):
        for e in event:
            if e <= averages[i] + (2*stds[i]) and e >= averages[i] - (2*stds[i]):
                counters[i] += 1
           
    confidence = [(c/N)*100 for c in counters]
    print(corr_error)
    print(error)
    print(confidence)


'''
Following class is a collection of functions for plotting,
The plots are all different visualizations for the event generator.
'''
class Event_Visualizer:
    def __init__(self) -> None:
        super().__init__() 

    '''
    Used for Plotting the distribution for all parameters
    '''
    def plot_distro(self,N,bins,Events,dimen = 1):
        for i in range(dimen):
            x = self.plot_distro_single(N,Events,i+1,bins)
            plt.hist(x,bins,label = str(dimen))
        plt.legend()
        plt.show()

    '''
    Plots the distribution for single parameter. dimen = parameter index+1
    '''
    def plot_distro_single(self,N,Events,dimen = 1,disp = False,bins = None):
        x = []
        for i in range(N):
            x.append(Events[i][dimen-1])

        if disp == True:
            plt.hist(x,bins,label = str(dimen))
            plt.legend()
            plt.show()

        return x
    
    '''
    Plots the correlation. used a subfunction for plotting covariance table
    '''
    def plot_relation(self,par1,par2,all_draws,disp = False):
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

    '''
    Generates correlation matrix in the form of a graph. 
    The diagonal holds gaussian distribution for each variable.
    The rest is a scatter plot against on another.
    '''
    def Plot_correlation_graph(self):
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
                    x = self.plot_distro_single(N,E,i+1)
                    ax_hist_var.hist(x,bins)
                else:
                    ax_scatter =  fig.add_subplot(gs[i, j])
                    np.append(ax_scatters,ax_scatter)
                    x,y = self.plot_relation(i,j,E)
                    ax_scatter.scatter(x,y,s=0.001)
        plt.show()
