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
center = nodes.loc[nodes.pos == 0].index[0]
upper_limit = min(center + n, len(nodes))
lower_limit = max(0, center - n)
N = nodes.iloc[lower_limit:upper_limit+1].reset_index(drop=True)    # Lista de nodos

# Estimar la información de distancias
dist_info = np.zeros((len(N), len(N)))
for i in range(len(N)):
    for j in range(len(N)):
        N_aux = N.iloc[min(i, j):max(i, j)]
        dist_info[i][j] = N_aux.dist.sum()

# Velocidades
v_min = float(info[vars.index('v_min')])  # km/h
v_max = float(info[vars.index('v_max')])  # km/h
v_step = float(info[vars.index('v_step')])  # km/h
V = np.linspace(v_min, v_max, num=int((v_max-v_min)/v_step)+1, endpoint=True, dtype=int)
v_river = float(info[vars.index('v_river')])  # metros por segundo

# Capacidad de carga del bote en kg
Q = float(info[vars.index('Q')])

# Pesos
Q_min = float(info[vars.index('Q_min')])
Q_step = float(info[vars.index('Q_step')])
W = np.linspace(Q_min, Q, num=int((Q-Q_min)/Q_step)+1, endpoint=True, dtype=int) # Se debe convertir a m/s

# Aqui crear el json con la información de las distancias
boat = Boat()
instance = []

# Crea los labels de la tabla de instancias
text = open("instance.txt", "a")
text.write(f"origin;destiny;v;w;d;t;e;flow\n")
text.close()

for v in V:
    for c in W:
        for i in range(len(dist_info)):
            for j in range(len(dist_info)):
                dist = dist_info[i][j] * 1000
                if N.iloc[i].pos > N.iloc[j].pos:
                    vel = (v / 3.6)   # pasarlo a metros por segundo
                    current = False
                else:
                    vel = (v / 3.6) + v_river   # Si va en contra corriente se le suma la velocidad del rio
                    current = True
                w = c
                t = dist / vel  # Se calcula el tiempo de viaje
                inst_pow = boat.drag_sav(vel, w)    # Se calcula la potencia instantanea con modelo del bote
                e = inst_pow * t / 3600     # Se calcula el consumo
                arc_info = {
                    "origin": int(i),
                    "destiny": int(j),
                    "v": float(v),
                    "w": float(c),
                    "d": float(dist),
                    "t": float(t),
                    "e": float(e),
                    "flow": current
                }
                instance.append(arc_info)
                text = open("instance.txt", "a")
                text.write(f"{int(i)};{int(j)};{float(v)};{float(c)};{float(dist)};{float(t)};{float(e)};{current}\n")
                text.close()

# Guarda en formato json
with open('instance.json', 'w') as archivo:
    json.dump(instance, archivo)

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

# Escribe el archivo con vectores adicionales
with open('instance2.txt', 'w') as archivo:
    archivo.write('W' + ':' + ','.join(map(str, W)) + '\n')
    archivo.write('V' + ':' + ','.join(map(str, V)) + '\n')
    archivo.write('d' + ':' + ','.join(map(str, d)) + '\n')
    archivo.write('E' + ':' + str(E) + '\n')
    archivo.write('Q' + ':' + str(Q) + '\n')
    archivo.write('f' + ':' + str(f) + '\n')
    archivo.write('s' + ':' + str(s) + '\n')
    archivo.write('T_min' + ':' + str(T_min) + '\n')
