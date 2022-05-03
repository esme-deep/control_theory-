import numpy as np
import matplotlib.pyplot as plt
import package_DBR
from IPython.display import display, clear_output

def LeadLag_RT(MV, Kp, Tlead, Tlag, Ts, PV, PVInit=0, method='EBD'):
    #slide 130
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K)) * PV[-1] + (Kp*K)/(1+K) * ((1 + (Tlead/Ts)) * MV[-1] - (Tlead/Ts) * MV[-2]))
            elif method == 'EFD':
                PV.append((1-K) * PV[-1] + Kp*K * ((Tlead/Ts) * MV[-1] + (1-Tlead/Ts) * MV[-2]))
            #elif method == 'TRAP':
                #PV.append()            
            else:
                PV.append((1/(1+K)) * PV[-1] + (Kp*K)/(1+K) * ((1 + (Tlead/Ts)) * MV[-1] - (Tlead/Ts) * MV[-2]))
    else:
        PV.append(Kp*MV[-1])