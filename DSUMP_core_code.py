import math
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from math import radians, cos, sin, asin, sqrt


filepath,flight_path,strength,range_low,range_high,interval_x,interval_z=st.getsetting()
# get max value of bin
z_max_add=gmz.getmaxz()

with h5py.File(filepath, mode='r') as f:
    # read data
    latvar = f['/gt'+str(flight_path)+strength+'/heights/lat_ph']
    latitude = latvar[:]
    lat_vr = [latvar.attrs['valid_min'], latvar.attrs['valid_max']]

    # lonvar = f['/gt1l/geolocation/reference_photon_lon']
    lonvar = f['/gt'+str(flight_path)+strength+'/heights/lon_ph']
    longitude = lonvar[:]
    lon_vr = [lonvar.attrs['valid_min'], lonvar.attrs['valid_max']]

    # We'll plot h.
    h_var = f['/gt'+str(flight_path)+strength+'/heights/h_ph']
    temp = h_var[:]
    units = h_var.attrs['units']
    units = units.decode('ascii', 'replace')
    long_name = h_var.attrs['long_name']
    long_name = long_name.decode('ascii', 'replace')

    quality = f['/gt'+str(flight_path)+strength+'/heights/quality_ph']
    quality0 = quality[:]
    confidence = f['/gt'+str(flight_path)+strength+'/heights/signal_conf_ph']
    confidence4 = confidence[:]

    all_x,all_y,all_z=[],[],[]
    for i in range(len(temp)):
        if (latitude[i] > range_low and latitude[i] < range_high and temp[i]<z_max_add+0.2):
            all_x.append(latitude[i])
            all_y.append(longitude[i])
            all_z.append(temp[i])

    # DSUMP
    def count_twogs_point_area_clength(mu1,mu2,n1,gs_y1,n2,gs_y2):
        gs_y_max1=np.max(gs_y1)
        gs_y_max2 = np.max(gs_y2)
        sum_gs_rate1 = 0
        sum_gs_rate2 = 0
        # upper surface
        for i in range(len(n1)):
            length1 = abs((z_range_high-mu1)-(len(n1)-i)*0.1)
            if(n1[i]<(gs_y1[i]+gs_y1[i+1])*0.5):
                sum_gs_rate1=sum_gs_rate1+n1[i]*0.1/(math.e**(length1*1))
            else:
                if(n1[i]>gs_y_max1 and n1[i]<gs_y_max2):
                    sum_gs_rate1 = sum_gs_rate1 + gs_y1[i]*0.1/(math.e**(length1*1))
        for i in range(len(n2)):
            length2 = abs(( mu2-z_range_low) - i * 0.1)
            if(n2[i]<(gs_y2[i]+gs_y2[i+1])*0.5):
                sum_gs_rate2=sum_gs_rate2+n2[i]*0.1/(math.e**(length2*1))
            else:
                if(n2[i]>gs_y_max2 ):
                    sum_gs_rate2 = sum_gs_rate2 + gs_y2[i]*0.1/(math.e**(length2*1))

        # lower surface
        for i in range(len(n1)):
            length1 = abs((z_range_high-mu1)-(len(n1)-i)*0.1)
            if(n1[i]<(gs_y1[i]+gs_y1[i+1])*0.5):
                sum_gs_rate1=sum_gs_rate1+n1[i]*0.1/(math.e**(length1*1))
            else:
                if(n1[i]>gs_y_max1):
                    sum_gs_rate1 = sum_gs_rate1 + gs_y1[i]*0.1/(math.e**(length1*1))
        for i in range(len(n2)):
            length2 = abs(( mu2-z_range_low) - i * 0.1)
            if(n2[i]<(gs_y2[i]+gs_y2[i+1])*0.5):
                sum_gs_rate2=sum_gs_rate2+n2[i]*0.1/(math.e**(length2*1))
            else:
                if(n2[i]>gs_y_max2 and n2[i]<gs_y_max1):
                    sum_gs_rate2 = sum_gs_rate2 + gs_y2[i]*0.1/(math.e**(length2*1))
        return sum_gs_rate1,sum_gs_rate2,(sum_gs_rate1+sum_gs_rate2)*0.5

    # Gaussian
    def twopeakgs(g_break_min,z):
        # g_break_min = -1
        gauss1, gauss2 = [], []
        for i in range(len(z)):
            if (z[i] > g_break_min):
                gauss1.append(z[i])
        for i in range(len(z)):
            if (z[i] < g_break_min):
                gauss2.append(z[i])

        mu1 = np.mean(gauss1)
        sigma1 = np.std(gauss1)
        mu2 = np.mean(gauss2)
        sigma2 = np.std(gauss2)

        print("mu1:", mu1, "mu2:", mu2)

        g1bin = int((z_range_high - g_break_min) / interval_z)
        g2bin = int((g_break_min - z_range_low) / interval_z)
        n1, bins1, gs_y1 = create_gs(gauss1, mu1, sigma1, g1bin)
        n2, bins2, gs_y2 = create_gs(gauss2, mu2, sigma2, g2bin)

        sgr1_i, sgr2_i, sgr_i = count_twogs_point_area_clength(mu1, mu2, n1, gs_y1, n2, gs_y2)
        gslist_all.append(sgr_i)


