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

    def solve_until(self, tf: float, h: float = 30) -> None:
        
        N = tf // h
        y_res = np.empty([3, N])
        y_res[:, 0] = self.y0

        for i in range(N-1):
            # Extraer temperaturas
            self.t += h
            yi = y_res[:, i]
            THi, TMi, TLi = yi
            self.historial_TH.append(THi)
            self.historial_TM.append(TMi)
            self.historial_TL.append(TLi)
            self.historial_t.append(self.t)
            
            # Actualizar termostato maestro
            if self.estado_H == 'ON' and THi > self.Ts:
                self.estado_H = 'OFF'

            elif self.estado_H == 'OFF' and THi < self.Ts - self.TdbH:
                self.estado_H = 'ON'

            # Actualizar termostato esclavo:
            # Si está encendido...
            if self.estado_L == 'ON':
                # y ya calentó lo suficiente...
                if TLi > self.Ts:
                    # apáguese
                    self.estado_L = 'OFF'
                
                # y no ha calentado lo suficiente pero el maestro se enciende...
                elif self.estado_H == 'ON':
                    # dele campo a él y espere
                    self.estado_L = 'WAITING'

            # Si está apagado...  
            elif self.estado_L == 'OFF':
                # Si debe encenderse y el maestro se lo permite
                if TLi < self.Ts - self.TdbL and self.estado_H == 'OFF':
                    # Enciéndase
                    self.estado_L = 'ON'

                # Si debe encenderse pero el maestro NO se lo permite
                elif TLi < self.Ts - self.TdbL and self.estado_H == 'ON':
                    # A esperar...
                    self.estado_L = 'WAITING'

            # Si está esperando...
            elif self.estado_L == 'WAITING':
                # y el maestro finalmente le permite encenderse...
                if self.estado_H == 'OFF':
                    # enciéndase
                    self.estado_L = 'ON'
                
            # Guardar la potencia
            if self.estado_H == 'ON' or self.estado_L == 'ON':
                self.historial_P.append(self.Pnom)
            else:
                self.historial_P.append(0)
            
            # Actualizar variables continuas
            k1 = self.y_dot(t=None, y=yi)
            k2 = self.y_dot(t=None, y=yi + 0.5*h*k1)
            k3 = self.y_dot(t=None, y=yi + 0.5*h*k2)
            k4 = self.y_dot(t=None, y=yi + 1.0*h*k3)
            
            y_res[:, i+1] = yi + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
            
# Definir las condiciones iniciales y el tiempo de integración
# Punto 1
horas = 12  
horizonte = horas * 60 * 60

'''
Creación del objeto tanque
'''
tanque = Tanque(Tamb=20)
sol_tanque = tanque.solve_until(horizonte, h=1)

plt.figure(0)
plt.plot(tanque.historial_t, tanque.historial_TH, label='th')
plt.plot(tanque.historial_t, tanque.historial_TM, label='tm')
plt.plot(tanque.historial_t, tanque.historial_TL, label='tl')
plt.legend()
plt.xlabel('Tiempo (h)')
plt.ylabel('Temperatura (grados C)')

# Punto 2
plt.figure(1)
plt.plot(tanque.historial_t, tanque.historial_P)
plt.xlabel('Tiempo (h)')
plt.ylabel('Potencia (kW)')
plt.show()


# Punto 3
calentadores = [Tanque(Tamb=random.choice([20, 40])) for i in range(5)]

for i, calentador in enumerate(calentadores):
    print(f'Resolviendo calentador número {i}')
    calentador.solve_until(horizonte, h=30)
    plt.figure(2)
    # plt.plot(calentador.historial_t, calentador.historial_TH, label='TH')
    # plt.plot(calentador.historial_t, calentador.historial_TM, label='TM')
    # plt.plot(calentador.historial_t, calentador.historial_TL, label='TL')
    plt.plot(calentador.historial_t,
             calentador.historial_P, label=f'Pot_cal {i}')
plt.xlabel('Tiempo (s)')
plt.ylabel('Potencia (kW)')
plt.legend()
plt.show()


# Punto 4
temperaturas = np.arange(50, 65, 5)
TH_resultado = []
TM_resultado = []
TL_resultado = []

for Tamb in temperaturas:
    tanque = Tanque(Tamb)
    tanque.solve_until(246060, h=60)
    TH_resultado.append(tanque.historial_TH)
    TM_resultado.append(tanque.historial_TM)
    TL_resultado.append(tanque.historial_TL)

for i in range(len(temperaturas)):
    plt.figure(3)
    plt.plot(tanque.historial_t,
             TH_resultado[i], label=f'Tamb={temperaturas[i]}°C')
plt.title('Temperaturas ambiente')
plt.xlabel('Tiempo (s)')
plt.ylabel('Temperatura (°C)')
plt.legend()
plt.show()
