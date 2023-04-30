import numpy as np
import matplotlib.pyplot as plt

# load data
# data should be in this order:
# [mom.z,mom.y,an_h,gamma_h,f_h,fp_h,df_h]
dt = np.load('runs/run1.npy', allow_pickle=True)

z = dt[0]
y = dt[1]

xr = y[0,:] + y[1,:]

plt.figure()
plt.plot(z, xr)
plt.show()