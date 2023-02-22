import numpy as np
import matplotlib.pyplot as plt
import random

# Parte 1 y 2
class Tanque:
    
    def __init__(self, Tamb: float):
        
        # Condiciones externas
        self.Tamb = Tamb
       
        # Parámetros del agua que ingresa
        self.Wr = 0.007848662004213181
        self.Tin = 30
        
        # Parámetros de los termostatos
        self.Ts = 60
        self.TdbH = 5
        self.TdbL = 10
        self.estado_H = random.choice(['ON', 'OFF'])
        self.estado_L = 'OFF' if self.estado_H == 'ON' else 'ON'

        # Consumo de potencia
        self.Pnom = random.choice([3e3, 4e3])#2e3, 2.5e3, 3e3, 3.5e3, 4e3])
        self.Qh1 = self.Pnom if self.estado_H == 'ON' else 0
        self.Qh2 = self.Pnom if self.estado_L == 'ON' else 0
        
        # Parámetros constructivos
        Cp1 = 2.531e5
        Cp2 = 2.531e5
        Cp3 = 2.531e5
        Gsc1 = 0.5265
        Gsc2 = 0.5263
        Gsc3 = 0.5265
        Ks1 = 0.4134
        Ks2 = 0.4106
        Ks3 = 0.4123
        
        self.Cp1 = 250235.83740734213
        self.Cp2 = self.Cp1#random.uniform(Cp2*0.9, Cp2*1.1)
        self.Cp3 = self.Cp1#random.uniform(Cp3*0.9, Cp3*1.1)
        self.Gsc1 = 0.5194449048715224
        self.Gsc2 = self.Gsc1#random.uniform(Gsc2*0.9, Gsc2*1.1)
        self.Gsc3 = self.Gsc1#random.uniform(Gsc3*0.9, Gsc3*1.1)
        self.Ks1 = 0.40759230446647454
        self.Ks2 = self.Ks1#random.uniform(Ks2*0.9, Ks2*1.1)
        self.Ks3 = self.Ks1#random.uniform(Ks3*0.9, Ks3*1.1)
        
        # Constante del agua
        self.p = 1
        self.cp = 4190    

        # Historial de resultados
        self.t = 0
        TH0 = self.Ts - random.random() * self.TdbH
        TL0 = self.Ts - random.random() * self.TdbL
        TM0 = (TL0 + TH0)/2
        self.y0 = np.array([TH0, TM0, TL0])
        self.historial_t = []
        self.historial_TH = []
        self.historial_TM = []
        self.historial_TL = []
        self.historial_P = []

    def y_dot(self, t: float, y: np.ndarray) -> np.ndarray:
        '''
        Devuelve derivada del vector de estados.
        
        y: array unidimensional
        '''
        
        TH, TM, TL = y

        Qh1 = self.Pnom if self.estado_H == 'ON' else 0
        Qh2 = self.Pnom if self.estado_L == 'ON' else 0

        Qf1 = self.Ks2*TM - self.Ks1*TH
        TH_dot = (1/self.Cp1) * (  self.p*self.Wr*self.cp*(TM - TH) \
                                 + self.Gsc1*(self.Tamb - TH) \
                                 + Qf1 + Qh1)
        
        Qf2 = self.Ks1*TH + self.Ks3*TL - 2*self.Ks2*TM
        TM_dot = (1/self.Cp2) * (  self.p*self.Wr*self.cp*(TL - TM) \
                                 + self.Gsc2*(self.Tamb - TM) \
                                 + Qf2)
            
        Qf3 = self.Ks2*TM - self.Ks3*TL
        TL_dot = (1/self.Cp3) * (  self.p*self.Wr*self.cp*(self.Tin - TL) \
                                 + self.Gsc3*(self.Tamb - TL) \
                                 + Qf3 + Qh2)
        
        return np.array([TH_dot,
                         TM_dot,
                         TL_dot])    

    
