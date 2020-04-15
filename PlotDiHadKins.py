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
import awkward1 as ak

def loadRootFileArray(file):
    events = uproot.open(file)
    return events['tree'].arrays(namedecode="utf-8")
   # return events['tree'].lazyarrays(namedecode="utf-8",entrysteps=500)
def loadRootFile(file):
    events = uproot.open(file)
    return events['tree'].lazyarrays(entrysteps=500)

#takes in the arrays from two root files (gotten e.g. from loadRootFileA) ) and returns an awkward array event structure
def returnEventStructure(unsmearedEventsA, smearedEventsA):
    arrays = {name: ak.from_awkward0(array) for name, array in smearedEventsA.items()}
    arraysTruth = {name: ak.from_awkward0(array) for name, array in unsmearedEventsA.items()}
    events=ak.zip({"evnum":arrays["evnum"], "true":ak.zip({"x": arraysTruth["x"],"y": arraysTruth["y"],"Q2": arraysTruth["Q2"],"pair":ak.zip({"Z": arraysTruth["Z"],"hadP": arraysTruth["hadP"],
                                                                               "hadPt": arraysTruth["hadPt"],
                                                           "hadEta": arraysTruth["hadEta"],"hadPhi": arraysTruth["hadPhi"]},depthlimit=1),
                                                            "PhPerp": arraysTruth["PhPerp"],"PhEta": arraysTruth["PhEta"],"PhPhi": arraysTruth["PhPhi"],"Ph": arraysTruth["Ph"]},depthlimit=1),
                   "rec": ak.zip({"x": arrays["x"],"y": arrays["y"],"Q2": arrays["Q2"],"PhPerp": arrays["PhPerp"],"PhEta": arrays["PhEta"],"PhPhi": arrays["PhPhi"],"Ph": arrays["Ph"],
                                              "pair":ak.zip({"Z": arrays["Z"],"hadP": arrays["hadP"],
                                                            "hadPt": arrays["hadPt"],
                                                           "hadEta": arrays["hadEta"],"hadPhi": arrays["hadPhi"]},depthlimit=1)},depthlimit=1)})
    return events
    
def plotXVsEtaJ(data):
    return

def plotXVsEtaPtJ(events):
   
    return


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
    etaEdges=[-3+i for i in range(7)]
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
    #don't cut off half of the upper and lower bins
    xMean[0]=xEdges[0]
    xMean[-1]=xEdges[-1]
    etaMean[0]=etaEdges[0]
    etaMean[-1]=etaEdges[-1]

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
    ptEdges=[]
    ptEdges.append(0.1)
    ptEdges.append(0.5)
    ptEdges.append(1.0)
    ptEdges.append(2.0)
    ptEdges.append(5.0)
    etaEdges=[-3+i for i in range(7)]
    #this stores the eta for each hadron in the pair, so we have to flatten this structure
    arEta=data["hadEta"]
    arEta=arEta.flatten();
    arPt=data["hadPt"]
    arPt=arPt.flatten();
    #get x for the other dimension, now of course we have only have as many x values
    arX=data["x"]
    arX=np.repeat(arX,2)
    
    hist2d=np.histogram2d(arX,arPt,bins=(xEdges,ptEdges),weights=arX)
    hist2d_unweighted=np.histogram2d(arX,arPt,bins=(xEdges,ptEdges))
    
    z=ma.masked_where(hist2d[0] <= 0, hist2d[0])
    z_unweighted=ma.masked_where(hist2d_unweighted[0] <= 0, hist2d_unweighted[0])
    z=z/z_unweighted
    etaMean=[(etaEdges[i]-etaEdges[i-1])/2+etaEdges[i-1]  for i in range(len(etaEdges)) if(i>0)]
    ptMean=[(etaEdges[i]-etaEdges[i-1])/2+etaEdges[i-1]  for i in range(len(etaEdges)) if(i>0)]
    #don't cut off half of the upper and lower bins
    ptMean[0]=xEdges[0]
    ptMean[-1]=xEdges[-1]
    etaMean[0]=etaEdges[0]
    etaMean[-1]=etaEdges[-1]
    
    
    cf=plt.contourf(etaMean, ptMean, z, alpha=.75, cmap='jet',locator=ticker.LogLocator(numticks=20,subs=range(1,10)))
  #  cf=plt.contourf(etaMean, np.log10(xMean), z, 8, alpha=.75, cmap='jet',locator=ticker.LogLocator())
  #  C = plt.contour(etaMean,xMean, z, 8, colors='black', linewidth=.5)
    fig.colorbar(cf)
    ax.set_title("mean $x$, $\eta$ vs $p_T$")
    plt.xlabel('$\eta$')
    plt.ylabel('$p_T$')
    
    return


#resolution in M, input are dataframes with true and smeared data
def resInM(data, smeared):
    arM=data["Mh"]
    arMSmeared=smeared["Mh"]
    
    mEdges=[]
    mEdges.append(0.0)
    mEdges.append(0.5)
    mEdges.append(1.0)
    mEdges.append(1.5)
    mEdges.append(2.0)
    mEdges.append(5.0)

    dfM=pd.DataFrame()
    dfM["nonSmeared"]=data["Mh"]
    dfM["diffToSmeared"]=dfM["nonSmeared"]-smeared["Mh"]
    dfM['bin']=pd.cut(dfM["nonSmeared"],mEdges)
    
    s=dfM.groupby(dfM['bin'])
    sMean=s.mean()
    sStd=s.std()
    
    plt.errorbar(sMean['nonSmeared'],sMean['diffToSmeared'],yerr=sStd['diffToSmeared'],fmt='v')
    plt.xlabel('$M_{Inv}$')
    plt.ylabel('$M$-$M_{smeared}$')
    
    return
