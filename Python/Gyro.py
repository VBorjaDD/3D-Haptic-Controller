import time
from icm20948 import ICM20948
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


imu = ICM20948(0x69)

G = np.zeros((1,3))

Gangles = np.zeros((1,3))

timesample = 0.1

print('measurement start')
for i in range (40):

	x, y, z = imu.read_magnetometer_data()
	ax, ay, az, gx, gy, gz = imu.read_accelerometer_gyro_data()
	gx = (gx - 0.6328)*1.53
	gy = (gy - 0.8427)*1.53
	gz = (gz - 0.2779)*1.53

	gyro = np.array([gx, gy, gz])

	Gangles = Gangles + gyro * timesample

	G = np.append(G, Gangles, axis = 0)

	time.sleep(timesample)


plt.plot(G)
plt.legend(['X', 'Y', 'Z'])
plt.savefig('angles.png')