import matplotlib.pyplot as plt
import numpy as np

t=np.linspace(0,2000000,100000)
premrna=np.zeros(100000)
mrna = np.zeros(100000)
mrnaout = np.zeros(100000)
car = np.zeros(100000)
k = np.zeros(100000)
a1 = 0.024
t1 = 151.35
d1 = 0.0033
t2 = 1800
v1 = 0.002
v2 = 0.00014
d2 = 0.000002778
i = 0
dt = 20
while i<99999:
    premrna_t = (a1 * 40 - (0.693/t1)*premrna[i] - (0.693/t2)*premrna[i])*dt + premrna[i]
    premrna[i+1] = premrna_t
    mrna_t = mrna[i]+ ((0.693/t1)*premrna[i] - v1*mrna[i])*dt
    mrna[i+1] = mrna_t
    mrnaout_t = mrnaout[i] + (v1*mrna[i]-d1*mrnaout[i])*dt
    mrnaout[i+1] = mrnaout_t
    car_t = car[i]+ (v2*mrnaout[i]-d2*car[i])*dt
    car[i+1] = car_t
    k[i] = (car[i+1]-car[i])/dt
    i+=1
plt.plot(t,premrna,'r')
plt.legend(["premRNA"],loc = 'lower right')
plt.xlabel('time/s')
plt.ylabel('amount')
plt.show()
plt.plot(t,mrna,'g')
plt.legend(["mRNA"],loc = 'lower right')
plt.xlabel('time/s')
plt.ylabel('amount')
plt.show()
plt.plot(t,mrnaout,'c')
plt.legend(["mRNAout"],loc = 'lower right')
plt.xlabel('time/s')
plt.ylabel('amount')
plt.show()
plt.plot(t,car,'violet')
plt.legend(["CAR"],loc = 'lower right')
plt.xlabel('time/s')
plt.ylabel('amount')
plt.show()
plt.plot(t,k,'violet')
plt.legend(["k"],loc = 'lower right')
plt.show()