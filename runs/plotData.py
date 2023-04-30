import numpy as np
import matplotlib.pyplot as plt
from momentSolver.PlottingUtility import PlottingUtility
# load data
# data should be in this order:
# [mom.z,mom.y,an_h,gamma_h,f_h,fp_h,df_h]
dt = np.load('runs/run0.npy', allow_pickle=True)

z = dt[0]
y = dt[1]

xr = y[0,:] + y[1,:] #x^2
yr= y[0,:] - y[1,:] #y^@
f_h= y[5,:]

print('This is Q', y[0,:])
plt.figure()

#plt.plot(z, xr)

plt.title('f_h')
plt.plot(f_h)# supposed to be plotting figure of merit versus iteration number in the optimization lop, is this the same as the x xis in the standard plots generated by this code
plt.ylabel('figure of merit')
#range(len(f_h)),
#plt.xlim([0,5e-10])
#plt.show()
#Plot of Error Plots
plt.figure()
Ratio=np.divide(xr,yr)
Difference=xr-yr
Error=np.divide(Difference**2,abs(xr))# my xr is negative.... so it's not x^2...
Error_num=np.sqrt(Error)
print(Error_num)
plt.figure()
plt.title('Ratio')
plt. xlabel('Z Postition')
plt.ylabel('X^2/Y^2')
plt.plot(z, Ratio, label='$<Ratio>$')
#plt.show()
plt.figure()
plt.title('Error')
plt.xlabel('Z Postition')
plt.ylabel('sqrt((X^2-Y^2)^2/(X^2))')
plt.plot(z, Error, label='$<Error>$')
# plt.show()
plt.show()