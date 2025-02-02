import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# PARAMETER SETTING
file_loc = 'hypoDD.reloc'
file_sta = 'hypoDD.sta'
phiA = np.radians(99.29 - 90)
phiB = phiA - np.pi / 2
xshift, yshift, box_l, box_w, axlim, minmag = 0, 0, 2, 0.5, 0.1, 10

# Read events
mdat = pd.read_csv(file_loc, sep='\s+', header=None)
mag = mdat[16]
lon = mdat[2].values
lat = mdat[1].values
x = mdat[4] / 1000
y = mdat[5] / 1000
z = mdat[3]
ex = mdat[7] / 1000
ey = mdat[8] / 1000
ez = mdat[9] / 1000
print(f'# of events = {len(mdat)}')

fid1 = pd.read_csv(file_sta, sep = '\s+', header = None)
sta = fid1[0].values
slat = fid1[1].values
slon = fid1[2].values
print(f'# of stations = {len(sta)}')

# Plotting
plt.figure(figsize=(13, 13))

# STATION PLOT (map view)
plt.subplot(2, 2, 1)
plt.scatter(lon, lat, marker='o', s=2, color='r')
plt.scatter(slon, slat, marker='v', s=2)
for i, st in enumerate(sta):
    plt.text(slon[i], slat[i], st, fontsize=8)
plt.title('STATION MAP')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
#plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)


# PLOT EVENTS (map view)
plt.subplot(2, 2, 2)
plt.scatter(x, y, s=1, color='r')
for i in range(len(x)):
    plt.plot([x[i] - ex[i], x[i] + ex[i]], [y[i], y[i]], color='r')
    plt.plot([x[i], x[i]], [y[i] - ey[i], y[i] + ey[i]], color='r')

# Box Location on Map
plt.plot([(-box_l/2)*np.cos(phiA)-box_w*np.cos(phiB)+xshift,
          (box_l/2)*np.cos(phiA)-box_w*np.cos(phiB)+xshift,
          (box_l/2)*np.cos(phiA)+box_w*np.cos(phiB)+xshift,
          (-box_l/2)*np.cos(phiA)+box_w*np.cos(phiB)+xshift,
          (-box_l/2)*np.cos(phiA)-box_w*np.cos(phiB)+xshift],
         [(box_l/2)*np.sin(phiA)+box_w*np.sin(phiB)+yshift,
          (-box_l/2)*np.sin(phiA)+box_w*np.sin(phiB)+yshift,
          (-box_l/2)*np.sin(phiA)-box_w*np.sin(phiB)+yshift,
          (box_l/2)*np.sin(phiA)-box_w*np.sin(phiB)+yshift,
          (box_l/2)*np.sin(phiA)+box_w*np.sin(phiB)+yshift])
plt.plot([-box_l/2*np.cos(phiA)+xshift,box_l/2*np.cos(phiA)+xshift],
         [box_l/2*np.sin(phiA)+yshift-0.05,-box_l/2*np.sin(phiA)+yshift-0.05], color = 'b')

plt.scatter(x[mag > minmag], y[mag > minmag], marker='o', color='r')
plt.xlim([-axlim, axlim])
plt.ylim([-axlim, axlim])
plt.title('MAP VIEW')
plt.xlabel('Distance [km]')
plt.ylabel('Distance [km]')
plt.grid(True)


# Cross Section A-A'
plt.subplot(2, 2, 3)
i = np.abs((x-xshift)*np.cos(phiB)-(y-yshift)*np.sin(phiB)) < box_w
if len(i) > 0:
    x0 = x[i] - xshift
    y0 = y[i] - yshift
    z0 = z[i]
    mag0 = mag[i]
    ex0 = ex[i]
    ey0 = ey[i]
    ez0 = ez[i]
    cusp0 = cusp[i]
    
    plt.plot((x0*np.cos(phiA)-y0*np.sin(phiA)), z0, '.', markersize=1, color='r')
    if err == 1:
        for i in range(len(x0)):
            plt.plot([(x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)-ex0[i]),
                             (x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)+ex0[i])],
                            [z0[i], z0[i]], color='r')
            plt.plot([(x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)),
                             (x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA))],
                            [z0[i]-ez0[i], z0[i]+ez0[i]], color='r')

    if len(id) > 0:
        for i in range(len(x0)):
            for k in range(len(id)):
                if id[k] == cusp[i]:
                    plt.plot((x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)), z0[i], 'o', markersize=10, color='g')

    plt.plot((x0[mag0 > minmag]*np.cos(phiA)-y0[mag0 > minmag]*np.sin(phiA)),
                   z0[mag0 > minmag], 'o', color='r')

    plt.title("Cross Section: A-A'")
    plt.xlabel('Distance [km]')
    plt.ylabel('Depth [km]')
    plt.axis('equal')
    plt.axis([-box_l/2, box_l/2,
                    min(np.min(z0)-((np.max(z0)-np.min(z0)+0.01)/5), np.mean(z0)-box_l/2),
                    max(np.max(z0)+((np.max(z0)-np.min(z0)+0.01)/5), np.mean(z0)+box_l/2)])
    plt.xlim([-0.2, 0.2])
    plt.ylim([6.9, 6.5])
    plt.grid(True)

# Box Location on Map for B-B'
plt.subplot(2, 2, 4)

phiA = 99.29 + 90
tmp = box_w
box_w = box_l / 2
box_l = tmp * 2
phiB = np.radians(phiA - 90)
phiA = np.radians(phiA - 90)

# Cross Section B-B'
i = np.abs((x-xshift)*np.cos(phiB)-(y-yshift)*np.sin(phiB)) < box_w
if len(i) > 0:
    x0 = x[i] - xshift
    y0 = y[i] - yshift
    z0 = z[i]
    mag0 = mag[i]
    ex0 = ex[i]
    ey0 = ey[i]
    ez0 = ez[i]
    cusp0 = cusp[i]
    
    plt.plot((x0*np.cos(phiA)-y0*np.sin(phiA)), z0, '.', markersize=1, color='r')
    if err == 1:
        for i in range(len(x0)):
            plt.plot([(x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)-ex0[i]),
                             (x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)+ex0[i])],
                            [z0[i], z0[i]], color='r')
            plt.plot([(x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)),
                             (x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA))],
                            [z0[i]-ez0[i], z0[i]+ez0[i]], color='r')

    if len(id) > 0:
        for i in range(len(x0)):
            for k in range(len(id)):
                if id[k] == cusp[i]:
                    plt.plot((x0[i]*np.cos(phiA)-y0[i]*np.sin(phiA)), z0[i], 'o', markersize=10, color='g')

    plt.plot((x0[mag0 > minmag]*np.cos(phiA)-y0[mag0 > minmag]*np.sin(phiA)),
                   z0[mag0 > minmag], 'o', color='r')

    plt.title("Cross Section: B-B'")
    plt.xlabel('Distance [km]')
    plt.ylabel('Depth [km]')
    plt.axis('equal')
    plt.axis([-box_l/2, box_l/2,
                    min(np.min(z0)-((np.max(z0)-np.min(z0)+0.01)/5), np.mean(z0)-box_l/2),
                    max(np.max(z0)+((np.max(z0)-np.min(z0)+0.01)/5), np.mean(z0)+box_l/2)])
    plt.xlim([-0.2, 0.2])
    plt.ylim([6.9, 6.5])
    plt.grid(True)

plt.tight_layout()
plt.show()