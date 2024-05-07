# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:55:34 2024

@author: sgsreidg
"""

import numpy as np
import pandas as pd # pandas for opening and interacting with .csv files
import os #for interacting with operative system opening files from a list
import matplotlib.pyplot as plt


def data_RMSD (filename , datapos): # start function - datapos is dataposition in the describe list
    tdata = pd.read_csv(filename) # Open the .csv file
    tempdescribe = tdata["Reference.RMSD"].describe() # describe the reference RMSD, provides mean among other values
    mean = tempdescribe [datapos] # value from description list obtained
    return (filename , mean) # return the mean RMSD and filename as outputs

def makefilelist(dirname): #taking a directory (a string) as input
    a=os.listdir(dirname) #returning a list of files as output
    for k in range(0,len(a)):
        a[k]=dirname+a[k]
    return a

mydatafolder='M:/Documents/3rd Year/366 project/Gold files/Results/' #where the data is found

a=makefilelist(mydatafolder)

ligandlist = ['12a' , '18c' , '3rdd' ,'7PMT' , '7R2I' , '7R2J' , '7R2L' , '7TH7'] # List of all ligands that can be used in printing results

filenames = [] # list for storing filenames created

def rawdata (filename): # start function
    tdata = pd.read_csv(filename) # Open the .csv file
    rawdata = tdata["Reference.RMSD"] # collect raw RMSD values
    return (rawdata) # return the mean RMSD and filename as outputs

def quartiles(values, q1, q3):# Sample Code from https://matplotlib.org/stable/gallery/statistics/customized_violin.html 
    upper_adjacent_value = q3 + (q3 - q1) * 1.5 # 
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, values[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, values[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels): #Sample code from https://matplotlib.org/stable/gallery/statistics/customized_violin.html
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels) # remove the steps
    ax.set_xlim(0.25, len(labels) + 0.75)

spacing = len(a) / 4 # total number of files divided by number of different scoring systems
space = int(spacing) # converting to interger so it can be used as a limit



for c in range (0 , space): # one loop for each molecule/RMSD calculated (calculated by variable spacing and converted into the interger variable space)
    I = c*4 #In order to get files at the correct interval for it to be the same scoring system
    
    ASP_RMSDs = [] # create empty lists to be filled with the RMSD values
    CHEMPLP_RMSDs = []
    ChemScore_RMSDs = []
    GoldScore_RMSDs = []
    
    for x in range(0 , 9): # for each of 10 data points in each file
        ASP_RMSDs.append(rawdata(a[0+I])[x]) # creating lists of every RMSD value calculated fo every experiment in the respective scoring function's list
        CHEMPLP_RMSDs.append(rawdata(a[1+I])[x])
        ChemScore_RMSDs.append(rawdata(a[2+I])[x])
        GoldScore_RMSDs.append(rawdata(a[3+I])[x])
    
    ASP_array = np.array(ASP_RMSDs) # converting lists to arrays so that they can be plotted
    CHEMPLP_array = np.array(CHEMPLP_RMSDs)
    ChemScore_array = np.array(ChemScore_RMSDs)
    GoldScore_array = np.array(GoldScore_RMSDs)


    big_data = [sorted(ASP_array) , sorted(CHEMPLP_array) , sorted(ChemScore_array) , sorted(GoldScore_array)] # put all data into a big list

    fig1 = plt.figure(1, figsize=(14, 8)) # set size of the image created
    ax1 = fig1.add_subplot(111)

    ax1.set_title(ligandlist[c]) # setting labels for 
    ax1.set_ylabel('Root Mean Squared Deviation')  # label the x axis
    ax1.violinplot(big_data)


    ax2 = ax1.violinplot(big_data , showmeans = False , showmedians = False , showextrema = False) # create ax2 for overlay, without means, medians or extrema of ax2 visible

    for pc in ax2['bodies']:
        pc.set_facecolor('#d1419e') # set the olour of the fill of the violins
        pc.set_edgecolor('black') #black outline to define edges
        pc.set_alpha(1) # opacity of the filled in violin set to 0.9 (scale is 0 to 1)

    quartile1, medians, quartile3 = np.percentile(big_data, [25, 50, 75], axis=1)

    whiskers = np.array([
        quartiles(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(big_data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 0]

    inds = np.arange(1, len(medians) + 1)
    ax1.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax1.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5) 
    ax1.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    labels = ['ASP', 'CHEMPLP', 'ChemScore', 'GoldScore']
    set_axis_style(ax1, labels)
    ax1.set_yticks(range(0 , 18)) # interval of y axis points set to 1 and the maximum point on the plot to 17

    plt.show()

