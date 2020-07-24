import numpy as np
import os
import random
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import numpy as np
import scipy as sci
import scipy.constants as scic
from Base_Station import *
from Room_Setting import *
import math as m
#from statsmodels.distributions.empirical_distribution import ECDF
def sum_arrays(array, rows, start_rows = 0):
    new_array = np.zeros((np.size(array, axis =1), np.size(array, axis =2)))
    for i in range(rows):
        new_array = new_array + array[i+start_rows,:,:]
    return new_array
#specifically just paths, not summing antenna elements
def antenna_reflections(reflex_array, antennas):
    paths = []
    ant_path = []
    for i in range(0,antennas*4):
        ant_path = reflex_array[i, :,:] + reflex_array[i +antennas*4, :,:] +reflex_array[i +antennas*8, :, :] + reflex_array[i +antennas*12,:,:]
        paths.append(ant_path)
    paths = np.array(paths)
    return paths

def CDF_generation(array, num_bins = 50):
    data = np.reshape(array, (1,-1))
    data = data[data !=0]
    print('without zeros = ',np.shape(data))   
    data = abs(data)
    #sorted_data = data[np.searchsorted(data, 0) :]
    counts, bin_edges =  np.histogram(data, bins = num_bins, normed = 'True')
    cdf =np.cumsum(counts)
    return cdf, bin_edges
def get_distances(list_points, xy):
    distances = []
    for i in range(len(list_points)):
        distances.append(np.linalg.norm(list_points[i] - xy))
    print(distances)
    return distances

def scatter_center(num_of_points, room_size, min_distance, room_edge):
    centroids = []
    x_range = np.array([room_edge, (room_size[0]-room_edge)])
    y_range = np.array([room_edge, (room_size[1]-room_edge)])
    i = 0 
    while i < num_of_points:
        x = round(random.uniform(x_range[0],x_range[1]), 3)
        y = round(random.uniform(y_range[0], y_range[1]),3)
        xy = np.array([x,y])
        if(i == 0):
            centroids.append(xy)

            i += 1
        else:
            dist = get_distances(centroids, xy)
            if all(q > min_distance for q in dist):

                centroids.append(xy)
                i += 1
    centroids = np.reshape(np.array(centroids), (num_of_points,2))
    return centroids
def scatter_clusters(points, num_in_cluster,  max_distance, min_distance):
    clusters = []
    radius = max_distance - min_distance
    for i in range(np.size(points, axis=0)):
        for j in range(num_in_cluster):
            alpha = 2*m.pi*random.random()
            r = min_distance + radius*m.sqrt(random.random())
            x = r*m.cos(alpha) + points[i,0]
            y = r*m.sin(alpha) + points[i,1]
            xy = (x,y)
            clusters.append(xy)
    clusters = np.reshape(np.array(clusters), (np.size(points, axis = 0)*num_in_cluster,2))
    return clusters    
def patchAntenna(angle_deg,wave_length, phi = 0):
    W = 0.5*wave_length
    L = W
    k = (2*m.pi)/wave_length
    phi = np.deg2rad(phi)
    angle = []
    for i in range(len(angle_deg)):
        angle.append(m.radians(angle_deg[i]))
    angle = np.array(angle)
    E_t = np.zeros(len(angle))
    for j in range(len(angle)):
        if phi == 0:
            E_t[j] = 10*m.log10(m.cos(((k*L)/2)*m.sin(angle[j])))
            
        else:
            E_t[j] = 10*np.log10((np.sin((k*W*np.sin(angle[j])*np.sin(phi))/2)/((k*W*np.sin(angle[j])*np.sin(phi))/2))
                        *np.cos((k*L/4)*(np.sin(angle[j]+phi)+np.sin(angle[j]-phi)))*np.cos(angle[j]))
            
    return E_t

    def calcAngle(BS_point, BS_ref, xxyy_points):
        lineA = np.array([[BS_point[0],BS_point[1]],[BS_ref[0],BS_ref[1]]] )
        lineB = np.array([[BS_point[0],BS_point[1]],[xxyy_points[0],xxyy_points[1]]])
        vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
        vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]

        angle1 = m.atan2(vA[1], vA[0])
        angle2 = m.atan2(vB[1], vB[0])
        ang_deg = (angle1-angle2) * 360 /(2*m.pi)
        return ang_deg
    #antenna patch eq
def roomSections(array, px, py):
    sect = 3
    divx = (len(px)+sect //2)//sect
    divy = (len(py)+sect//2)//sect
    corner = array[4:divx, 4:divy]
    center = array[divx+1:divx*2, divy+1:divy*2]
    yside = array[4:divx, divy+1:divy*2]
    xside = array[divx+1:divx*2, 4:divy]
    print('Corner average: ', np.average(corner))
    print('Center average: ', np.average(center))
    print('Y-side average: ', np.average(yside))
    print('X-side average: ', np.average(xside))
def main():

    room = np.array([30,20]) 
    antennas =8


    frequency = 5.9e9
    wave_len = (3e8/frequency)
    bs_spacing = wave_len/2 +0.005
    resolution = 0.5
    room_edge = 1
 
    #Room features set up
    nodes = Room_Setting(room[0], room[1], frequency, resolution, room_edge)
    BS = Base_Station(room, antennas, frequency, bs_spacing)
    BS.get_location()

    base_side =  BS.get_sidewallBS()
    base_corner = BS.get_cornerBS()
    sideBSreflect = BS.mirroring(base_side)
    cornerBSreflect = BS.mirroring(base_corner)
    sideRefRef = BS.mirroring(BS.BS_sidewall_ref)
    cornerRefRef = BS.mirroring(BS.BS_corner_ref)
    xx, yy = np.meshgrid(nodes.px, nodes.py, indexing = 'ij')
    xxyy = np.concatenate((np.reshape(xx,(-1,1)),np.reshape(yy,(-1,1))), axis = 1)


    centers = np.load(r"ScatteringPoints\30points.npy")
    #clusters = np.load("20ptclusters.npy")
    centers = np.array(centers)
    #clusters = np.array(clusters)
    ant8_30pt_corner = np.load('dBpaths\\8antennas30pts1Corner_dB32.npy')
    corner8antenna = np.load('dBpaths\\8antennaCornerdB.npy')
    cornerScattering = ant8_30pt_corner +corner8antenna
    cornerScattering = cornerScattering*np.conj(cornerScattering)
    cornerScattered_single = sum_arrays(cornerScattering, 8)
    cornerScattered_single2 = sum_arrays(cornerScattering, 8,8)
    cornerScattered_single3 = sum_arrays(cornerScattering, 8,16)
    cornerScattered_single4 = sum_arrays(cornerScattering, 8,24)
    cornerScattered = sum_arrays(cornerScattering, np.size(cornerScattering, axis = 0))

    ant8_30pt_side = np.load('dBpaths\\8antennas30pts1Side_dB32.npy')
    side8antenna = np.load('dBpaths\\8antennaSideddB.npy')
    sideScattering = ant8_30pt_side +side8antenna
    sideScattering = sideScattering*np.conj(sideScattering)
    sideScattered = sum_arrays(sideScattering, np.size(sideScattering, axis = 0))

    sideScattered_single = sum_arrays(sideScattering, 8)
    sideScattered_single2 = sum_arrays(sideScattering, 8,8)
    sideScattered_single3 = sum_arrays(sideScattering, 8,16)
    sideScattered_single4 = sum_arrays(sideScattering, 8,24)

    cdfs, bins = CDF_generation(10*np.log10(sideScattered))
    cdfs_bs2, bins_bs2 = CDF_generation(10*np.log10(sideScattered_single2))
    cdfs_bs1, bins_bs1 = CDF_generation(10*np.log10(sideScattered_single))
    cdfs_bs3, bins_bs3 = CDF_generation(10*np.log10(sideScattered_single3))
    cdfs_bs4, bins_bs4 = CDF_generation(10*np.log10(sideScattered_single4))


    cdfc_bs1, binc_bs1 = CDF_generation(10*np.log10(cornerScattered_single))
    cdfc_bs2, binc_bs2 = CDF_generation(10*np.log10(cornerScattered_single2))
    cdfc_bs3, binc_bs3 = CDF_generation(10*np.log10(cornerScattered_single3))
    cdfc_bs4, binc_bs4 = CDF_generation(10*np.log10(cornerScattered_single4))
    cdfc, binc = CDF_generation(10*np.log10(cornerScattered))
    
    fig_cdf, axCDF = plt.subplots(figsize = [8,4])
    axCDF.semilogy(bins[1:], cdfs/cdfs[-1], label = 'Side-all', linestyle = ':', color = 'tab:blue')
    axCDF.plot(binc[1:], cdfc/cdfc[-1], label = 'Corner-all',color = 'tab:blue')

    axCDF.semilogy(bins_bs1[1:], cdfs_bs1/cdfs_bs1[-1], label = 'Side1', linestyle = ':', color = 'tab:red')
    axCDF.plot(binc_bs1[1:], cdfc_bs1/cdfc_bs1[-1], label = 'Corner1',color = 'tab:red')

    axCDF.plot(bins_bs2[1:], cdfs_bs2/cdfs_bs2[-1], label = 'Side2',linestyle = ':', color = 'tab:orange')
    axCDF.plot(binc_bs2[1:], cdfc_bs2/cdfc_bs2[-1], label = 'Corner2',color = 'tab:orange')

    axCDF.plot(bins_bs3[1:], cdfs_bs3/cdfs_bs3[-1], label = 'Side3',linestyle = ':', color = 'tab:green')
    axCDF.plot(binc_bs3[1:], cdfc_bs3/cdfc_bs3[-1], label = 'Corner3',color = 'tab:green')

    axCDF.plot(bins_bs4[1:], cdfs_bs4/cdfs_bs4[-1], label = 'Side4',linestyle = ':', color = 'tab:purple')
    axCDF.plot(binc_bs4[1:], cdfc_bs4/cdfc_bs4[-1], label = 'Corner4',color = 'tab:purple')

    axCDF.set_title('CDF of Channel of Room for 30 Scattering Points, Single BS\'s vs All Channel')
    axCDF.set_xlabel('FSPL, dB scale')
    axCDF.grid()
    axCDF.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
    axCDF.legend()


    print(np.shape(cornerScattered))


    fig, ax = plt.subplots( figsize =[8,4])
    plotz = ax.pcolor(xx,yy,-10*np.log10(abs(sideScattered_single[:,:])), vmin = 0, vmax = 50)
    cbar = fig.colorbar(plotz)
    cbar.set_label('dB')
    ax.set_xlabel('meters')
    ax.set_ylabel('meters')
    x = [0,5,10, 15,20,25,30]
    y = [0,5,10,15,20]
    ax.xaxis.set_ticks(x)
    ax.yaxis.set_ticks(y)
    ax.scatter(base_side[0:7,0], base_side[0:7,1], s =2, color = 'r')
    ax.set_title('Coverage Map with 30 Scattering Points with FSPL For Sidewall Base Station')
    
    plt.show()
    #np.save('dBpaths\%dantennas%dpts9Side_dB32.npy' %(antennas,len(centers)), sideScattering)

    #np.save('dBpaths\%dantennas%dpts9Corner_dB32.npy' %(antennas,len(centers)), cornerScattering)

   
if __name__ == "__main__":
    main()


    







    # room = np.array([30,20]) 
    # antennas =64


    # frequency = 5.9e9
    # wave_len = (3e8/frequency)
    # bs_spacing = wave_len/2 +0.005
    # resolution = 0.5
    # room_edge = 1
 
    # #Room features set up
    # nodes = Room_Setting(room[0], room[1], frequency, resolution, room_edge)
    # BS = Base_Station(room, antennas, frequency, bs_spacing)
    # BS.get_location()

    # base_side =  BS.get_sidewallBS()
    # base_corner = BS.get_cornerBS()
    # sideBSreflect = BS.mirroring(base_side)
    # cornerBSreflect = BS.mirroring(base_corner)
    # sideRefRef = BS.mirroring(BS.BS_sidewall_ref)
    # cornerRefRef = BS.mirroring(BS.BS_corner_ref)
    # xx, yy = np.meshgrid(nodes.px, nodes.py, indexing = 'ij')
    # xxyy = np.concatenate((np.reshape(xx,(-1,1)),np.reshape(yy,(-1,1))), axis = 1)
    # cornerReflectRef = BS.mirroring(BS.BS_corner_ref)
    # sidewallReflectRef   = BS.mirroring(BS.BS_sidewall_ref)


    # slopes = (base_side[:,1]-BS.BS_sidewall_ref[:,1])/(base_side[:,0]-BS.BS_sidewall_ref[:,0])
    # print(slopes)
    
    # #testing
    # #angle_test = BS.calcAngle(base_side[3,:],BS.BS_sidewall_ref[3,:],[7.8,9] )    
    # los_paths_corner = BS.los_path(nodes.px, nodes.py, xxyy, location = 'corner')
    # ref_paths_corner = BS.wall_reflect(nodes.px, nodes.py, xxyy, location = 'corner')
    # attenuatedLOScorner = BS.attenuation(los_paths_corner, base_corner, BS.BS_corner_ref,wave_len, xxyy, nodes.px, nodes.py)
    # attenuatedreflectedcorner = BS.attenuation(ref_paths_corner, cornerBSreflect,  cornerReflectRef, wave_len, xxyy, nodes.px, nodes.py)
    # los_paths_side = BS.los_path(nodes.px, nodes.py, xxyy)
    # ref_paths_side = BS.wall_reflect(nodes.px, nodes.py, xxyy)
    # attenuatedLOSsidewall =BS.attenuation(los_paths_side, base_side, BS.BS_sidewall_ref, wave_len, xxyy, nodes.px, nodes.py)
    # attenuatedreflectedsidewall = BS.attenuation(ref_paths_side,sideBSreflect, sidewallReflectRef,wave_len,xxyy, nodes.px, nodes.py )
    # cornerRefPaths = antenna_reflections(attenuatedreflectedcorner, antennas)
    # sidewallRefPaths = antenna_reflections(attenuatedreflectedsidewall, antennas)

    # totalCornerPaths = cornerRefPaths +attenuatedLOScorner
    # totalSidewallPaths = sidewallRefPaths +attenuatedLOSsidewall
    # combinedCorner = sum_arrays((totalCornerPaths), antennas*4)
    # combinedSide = sum_arrays((totalSidewallPaths), antennas*4)
    # np.save('dBpaths/%dantennaCornerdB' % antennas, totalCornerPaths)
    # np.save('dBpaths/%dantennaSideddB'% antennas, totalSidewallPaths)
    # cdf_side, side_bins = CDF_generation(10*np.log10(abs(combinedSide)))
    # cdf_corner, corner_bins = CDF_generation(10*np.log10(abs(combinedCorner)))

    # fig, ax = plt.subplots()
    # ax.contourf(xx,yy,abs(combinedSide),levels = 50)
    # ax.legend()
    # cfig, cax = plt.subplots()
    # cax.semilogy(side_bins[1:], cdf_side/cdf_side[-1], label = 'Sidewall Coverage')
    # cax.plot(corner_bins[1:], cdf_corner/cdf_corner[-1], label = 'Corner Coverage')
    # cax.legend()
    # cax.set_title("CDFs of LOS and Reflected Components for %d Antennas" % antennas)

    # cdf_side, side_bins = CDF_generation(combinedSide)
    # cdf_corner, corner_bins = CDF_generation(combinedCorner)