from math import sqrt;
from numpy import arctanh;


def mytest2():
    print ("Hello world");
    return;

def mytest():
    print("second test");
    return;

def getEta(px,py,pz):
    pT=sqrt(px*px+py*py);
    pTotal=sqrt(px*px+py*py+pz*pz);
    eta=arctanh(pz/pTotal);

    return eta;


def getPhi(px,py,pz):
    return phi;
