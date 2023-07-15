import numpy as np
import sys,os
sys.path.append('../Predicting-Binding-Energy')
import parameter_generator as pg


'''
Evaluates event generator.
Keeps track of how many events are within two standard deviations.
Prints the percentage.
'''

def evaluate_event_generator():
    all_events = []
    N = 100000
    for i in range(N):
        Event = np.array(pg.random_params())
        all_events.append(Event)

    all_events = np.array(all_events).T
    averages = list(map(np.mean,all_events))
    stds = list(map(np.std,all_events))
    corr = np.corrcoef(all_events)

    corr_error = []
    for i in range(10):
        b = []
        for j in range(10):
            b.append(round (abs ( (pg.R[i][j] - corr[i][j]) / pg.R[i][j] * 100 ), 2 ) )
        corr_error.append(b)

    error = [round(abs(a[0] - a[1])/a[1]*100 ,3) for a in zip(stds,pg.C)]

    counters = [0] * 10

    for i,event in enumerate(all_events):
        for e in event:
            if e <= averages[i] + (2*stds[i]) and e >= averages[i] - (2*stds[i]):
                counters[i] += 1
           
    confidence = [(c/N)*100 for c in counters]
    
    print(f"\tNumber of events = {N} \n")
    print("\tRelative Error Percentage of Standard Deviations")
    print(error)
    
    print("\n\tConfidence of all the parameters")
    print(confidence)

    print("\n\tRelative Error of correlation coefficients")
    for i in range(10):
        print(corr_error[i][:i+1])

dir = os.getcwd() + '/evaluating event generator'
evaluate_event_generator()

print("\n\t Building distribution plots!")
pg.Monte_Carlo_events(1000,100,True,10,dir)
EV = pg.Event_Visualizer()
EV.Plot_correlation_graph(dir)