import cv2
import numpy as np
import sys

print(sys.argv[1],sys.argv[2])
image = cv2.imread(sys.argv[1])
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
th3 = cv2.adaptiveThreshold(gray,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
                    cv2.THRESH_BINARY,45,15)

cv2.imwrite(sys.argv[2], th3)
