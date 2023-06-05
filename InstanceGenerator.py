from Objects import Boat
import numpy as np
import pandas as pd
import json

# Abrir oja de parametros
with open('parameters.txt', 'r') as file:
    lines = file.readlines()

vars = []
info = []
for line in lines:
    data = line.strip().split(':')
    vars.append(data[0])
    info.append(float(data[1]))

# Cargar la información de todos los nodos
nodes = pd.read_csv('instance.csv', sep=';')
nodes.dist = nodes.dist.astype(float)
n = int(info[vars.index('n')])  # Numero de nodos a utilizar

# Filtrar el numero de nodos seleccionado por el usuario
pos = int(info[vars.index('pos')])
center = nodes.loc[nodes.node == 'Puerto Calvo'].index[0]   # Nodo del deposito
if pos == -1:
    limit = center - n
    N = nodes.iloc[limit:center+1].sort_values(by='pos', ascending=False).reset_index(drop=True)
else:
    limit = center + n
    N = nodes.iloc[center:limit+1].reset_index(drop=True)

# Estimar la información de distancias
dist_info = np.zeros((len(N), len(N)))

for i in range(len(N)):
    for j in range(len(N)):
        if i != j:
            N_aux = N.iloc[min(i, j)+1:max(i, j)+1]
            dist_info[i][j] = N_aux.dist.sum()
        else:
            dist_info[i][j] = 0


# Velocidades
v_min = float(info[vars.index('v_min')])  # km/h
v_max = float(info[vars.index('v_max')])  # km/h
v_step = float(info[vars.index('v_step')])  # km/h
V = np.linspace(v_min, v_max, num=int((v_max-v_min)/v_step)+1, endpoint=True, dtype=int)
v_river = float(info[vars.index('v_river')])  # km/h

# Capacidad de carga del bote en kg
Q = float(info[vars.index('Q')])

# Pesos
Q_min = float(info[vars.index('Q_min')])
Q_step = float(info[vars.index('Q_step')])
W = np.linspace(Q_min, Q, num=int((Q-Q_min)/Q_step)+1, endpoint=True, dtype=int) # Se debe convertir a m/s

# Aqui crear el diccionario con la información de las distancias
boat = Boat()
matrix = []

for v in V:
    for c in W:
        for i in range(len(dist_info)):
            for j in range(len(dist_info)):
                dist = dist_info[i][j] * 1000
                if N.iloc[i].pos > N.iloc[j].pos:
                    vel = (v - v_river) / 3.6  # pasarlo a metros por segundo, va a favor de la corriente
                    current = False
                else:
                    vel = (v + v_river) / 3.6   # Si va en contra corriente se le suma la velocidad del río
                    current = True
                w = c
                t = dist / (v / 3.6)  # Se calcula el tiempo de viaje con respecto a la velocidad determinada
                inst_pow = boat.drag_sav(vel, w)    # Se calcula la potencia instantanea con modelo del bote
                e = inst_pow * t / 3600     # Se calcula el consumo
                arc_info = {
                    "origin": int(i),
                    "destiny": int(j),
                    "v": float(v),
                    "w": float(c),
                    "dist": float(dist),
                    "t": float(t),
                    "e": float(e),
                    "flow": current
                }
                matrix.append(arc_info)

# Demanda en peso de cada nodo en orden
d = np.random.choice(W, size=n, replace=True)

# Capacidad de la batería
E = float(info[vars.index('E')])     # kWh

# Tasa de carga de la batería
f = float(info[vars.index('f')])  # kWh de carga

# Tasa de descarga de mercancía
s = float(info[vars.index('s')])    # horas/kg

# Tiempo minimo de permanencia del bote en el puerto entre viajes
T_min = float(info[vars.index('T_min')])

# Promesa de entrega para la demanda
service_upp_lim = float(info[vars.index('service_upp_lim')])
service_lower_lim = float(info[vars.index('service_lower_lim')])
mean_v = np.mean(V)     # Velocidad media
farthest_place = dist_info[0][-1]   # Lugar más lejano en km
mean_t = farthest_place * 1000 / (mean_v / 3.6)     # tiempo en segundos hacia el lugar más lejano
l = np.random.uniform(low=service_lower_lim, high=service_upp_lim, size=n) * mean_t

# Guarda en formato json
instance = {
    "matriz": matrix,
    "W": W.tolist(),
    "V": V.tolist(),
    "d": d.tolist(),
    "E": E,
    "Q": Q,
    "f": f,
    "s": s,
    "T_min": T_min,
    "l": l.tolist(),
}

with open('instance.json', 'w') as archivo:
    json.dump(instance, archivo)