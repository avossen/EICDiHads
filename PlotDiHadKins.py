#maybe we want to try root plotting...
import ROOT
import uproot
import awkward
from awkward import JaggedArray
import plotly.express as px
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from math import log,exp
from numpy import ma
import seaborn as sne

def loadRootFile(file):
    events = uproot.open(file)
    return events['tree'].lazyarrays(entrysteps=500)



#do regular 1d plot, takes some panda
def plotXVsEta(data):
    #let's do some sort of countour/2d plot here
    #define x, eta edges (x in log)
    #histogram the stuff and do matplotlib
   # print("in plot x vs eta")
#    xEdges=[];
#    xEdges=[10**(-3+i/3) for i in range(8)]
    xEdges=np.logspace(-3,-1,7)
    #xEdges.append(0.001)
    #xEdges.append(0.02)
    #xEdges.append(0.05)
    #xEdges.append(0.8)
    
    #xEdgesHigh=xEdgesLow[1:]
    etaEdges=[-3+i for i in range(6)]
    #this stores the eta for each hadron in the pair, so we have to flatten this structure
    arEta=data["hadEta"]
    arEta=arEta.flatten();
    #get x for the other dimension, now of course we have only have as many x values
    arX=data["x"]
    arX=np.repeat(arX,2)
    hist2d=np.histogram2d(arX,arEta,bins=(xEdges,etaEdges))
    #should be done better, namely at mean of log
    xMean=[exp((log(xEdges[i])-log(xEdges[i-1]))/2+log(xEdges[i-1]))  for i in range(len(xEdges)) if(i>0)]
    etaMean=[(etaEdges[i]-etaEdges[i-1])/2+etaEdges[i-1]  for i in range(len(etaEdges)) if(i>0)]
    fig, ax=plt.subplots(1,1)
    #histogram returns histogram plus axes...
    z=ma.masked_where(hist2d[0] <= 0, hist2d[0])
    cf=plt.contourf(etaMean, xMean, z, alpha=.75, cmap='jet',locator=ticker.LogLocator(numticks=20,subs=range(1,10)))
  #  cf=plt.contourf(etaMean, np.log10(xMean), z, 8, alpha=.75, cmap='jet',locator=ticker.LogLocator())
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(subs=range(1,10)))
    ax.yaxis.set_minor_locator(ticker.LogLocator(numticks=30))
  #  C = plt.contour(etaMean,xMean, z, 8, colors='black', linewidth=.5)
    fig.colorbar(cf)
    ax.set_title("$\eta$ vs $x$")
    plt.xlabel('$\eta$')
    plt.ylabel('$x$')
    return

#do regular 1d plot, takes some panda
def plotXVsEtaKDE(data):
    arEta=data["hadEta"]
    arEta=arEta.flatten();
    #get x for the other dimension, now of course we have only have as many x values
    arX=data["x"]
    arX=np.repeat(arX,2)
    fig, ax=plt.subplots(1,1)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(subs=range(1,10)))
    ax.yaxis.set_minor_locator(ticker.LogLocator(numticks=30))
    sns.kdeplot(arEta,arX,shade=True, clip=((-3,3),(0.001,0.1)))

    #  C = plt.contour(etaMean,xMean, z, 8, colors='black', linewidth=.5)
    fig.colorbar(cf)
    ax.set_title("$\eta$ vs $x$")
    plt.xlabel('$\eta$')
    plt.ylabel('$x$')
    return

#do some sort of 2d plot, takes some panda
def plotXVsEtaPt(data):
    return


    
#resolution in M, input are dataframes with true and smeared data
def resInM(data, smeared):
    return
