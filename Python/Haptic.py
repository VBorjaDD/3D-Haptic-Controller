import RPi.GPIO as GPIO  
import time
from icm20948 import ICM20948
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

GPIO.setmode(GPIO.BCM)

imu = ICM20948(0x69)


GPIO.setup(18,  GPIO.OUT)
GPIO.setup(23,  GPIO.OUT)
GPIO.setup(24, 	GPIO.OUT)
GPIO.setup(25, 	GPIO.OUT)

motor1 = GPIO.PWM(18, 48)
motor2 = GPIO.PWM(23, 48)
motor3 = GPIO.PWM(24, 48)
motor4 = GPIO.PWM(25, 48)


motor1.start(5)
motor2.start(5)
motor3.start(5)
motor4.start(5)

time.sleep(5)

for i in range(30):
	motor1.ChangeDutyCycle(4.4 + 0.1*i)
	motor2.ChangeDutyCycle(4.4 + 0.1*i)
	motor3.ChangeDutyCycle(4.4 + 0.1*i)
	motor4.ChangeDutyCycle(4.4 + 0.1*i)
	time.sleep(0.1)


a = 0.8
b = 0.15
c = 0.19
d = 0.03

k = np.array([[		0,	 	-a, 	 b,		0,		-c,		d],
			[		a,		0,		-b,		c,		0,		-d],
			[		0,		 a,	 	b,		0, 		c,		d],
			[		-a,		0,		-b,		-c,		0,		-d]])

R = np.array([	[-0.7071, -0.7071],
				[ 0.7071, -0.7071]])


percent = 0.06
Angle = np.zeros(3)
timesample = 0.015
err = np.zeros(6)
print('start of the controller')

for i in range(3000):

	start = time.time()

	ax, ay, az, gx, gy, gz = imu.read_accelerometer_gyro_data()

	a = np.array([ax, ay])
	a = np.matmul(R,a)

	AccY = -(np.arctan(a[0] / np.sqrt(np.power(a[1], 2) + np.power(az, 2))))
	AccX = (np.arctan(a[1] / np.sqrt(np.power(a[0], 2) + np.power(az, 2))))

	Acc = np.array([AccX, AccY, 0])
	
	gx = (gx - 0.6030)/37.44
	gy = (gy - 0.8427)/37.44
	gz = (gz - 0.1221)/37.44

	gyro = np.array([gx,gy])
	gyro = np.matmul(R,gyro)

	gyro = np.append(gyro, gz)

	Gangles = np.add(Angle, timesample * gyro)

	Angle = np.add(0.96 * Gangles, 0.04 * Acc) #in rad

	err = np.add(err,np.append(Angle,np.zeros(3)))
	state = np.append(Angle, gyro)
	
	t = np.add(-(1-percent)*np.matmul(k, state),- percent * np.matmul(k, err))*1600
	t = np.add(t, np.ones(4)*0.0)

	t[t < 0] = 0

	dc = np.sqrt(t/255) + 4.4

	motor1.ChangeDutyCycle(dc[0])
	motor2.ChangeDutyCycle(dc[1])
	motor3.ChangeDutyCycle(dc[3])
	motor4.ChangeDutyCycle(dc[2])
	
	end = time.time()
	elapsed = end - start
	time.sleep(timesample - elapsed)


motor1.stop()
motor2.stop()
motor3.stop()
motor4.stop()

GPIO.cleanup()

#####

