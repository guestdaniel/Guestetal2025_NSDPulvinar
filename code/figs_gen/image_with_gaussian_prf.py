from PIL.Image import ImageTransformHandler
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

img = plt.imread("example_image.jpg.png")

# define normalized 2D gaussian
def gaus2d(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

x = np.linspace(-5, 5, 937)
y = np.linspace(-5, 5, 930)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D
z = gaus2d(x, y, mx=1.75, my=-1.75, sx=1.25, sy=1.25)

#plt.subplot(1, 3, 1)
#plt.imshow(img)
#plt.subplot(1, 3, 2)
#plt.imshow(z)
#plt.subplot(1, 3, 3)
#plt.imshow(tensordot()

output = np.ones((930, 937, 4))
#output[:, :, 2] += 1
output[:, :, 3] = (1 - z/np.max(z))

plt.imsave("temp.png", output)
