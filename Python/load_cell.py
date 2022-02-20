import RPi.GPIO as GPIO  
import time
import numpy as np
from hx711 import HX711

GPIO.setmode(GPIO.BCM)


hx = HX711(dout_pin=21, pd_sck_pin=20)

sum = 0

for i in range(100):
  read = (hx._read() + 101934)*100/20890
  print(read)
  time.sleep(0.1)
  sum = sum + read




GPIO.cleanup()

