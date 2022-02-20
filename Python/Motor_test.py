import RPi.GPIO as GPIO  
import time
import numpy as np
from hx711 import HX711
import matplotlib
import matplotlib.pyplot as plt

GPIO.setmode(GPIO.BCM)
hx = HX711(dout_pin=21, pd_sck_pin=20)

GPIO.setup(18, GPIO.OUT)

soft_pwm = GPIO.PWM(18, 48)

soft_pwm.start(4.4)
n=4.4

T = np.array([])
D = np.array([])
P =	np.array([])
oldread = 0

time.sleep(5)

for i in range(40):
	soft_pwm.ChangeDutyCycle(n)
	n = n + 0.2
	read = (hx._read() + 101934)*100/20890
	if np.abs(oldread - read)>50:
	read = oldread
	T = np.append(T,read)
	D = np.append(D,n)
	time.sleep(0.5)
	oldread = read
		print('Amps: ')
		while True:
			try:
				current = float(input())
				break
			except ValueError:
				print("Oops!  That was no valid number.  Try again...")

		power = 30 * current
		P = np.append(P,power)
		print('DutyCyle of: ', n,'% has a thrust of :', read ,'and pulls a power of', power, 'W')
		print('DutyCyle of: ', n,'% has a thrust of :', read )

for i in range(40):
	soft_pwm.ChangeDutyCycle(n)
	n = n - 5
	time.sleep(0.5)
	if n<0:
		break

soft_pwm.stop()
GPIO.cleanup()


fig, ax1 = plt.subplots()
ax1.plot(D, T, 'g-')
ax2 = ax1.twinx()
ax1.plot(D, T, 'g-')
ax2.plot(D, P, 'b-')

ax1.set_xlabel('DutyCycle [%]')
ax1.set_ylabel('Thrust [g]', color='g')
ax2.set_ylabel('Power [W]', color='b')

plt.show()
plt.xlabel('DutyCycle [%]')
plt.savefig('Thrust_30V.png')