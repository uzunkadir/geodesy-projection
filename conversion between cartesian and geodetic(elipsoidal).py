from math import sin,cos,atan2,degrees
from math import radians as rad

a=6378137.0000
b=6356752.3141
e=((a**2-b**2)/(a**2))**(0.5)

def ell2xyz(latitude,longtitude,ellHeight):
    
    N=a/((1-e**2*sin(rad(latitude))**2)**(0.5))
    
    x=(N+ellHeight)*cos(rad(latitude))*cos(rad(longtitude))
    y=(N+ellHeight)*cos(rad(latitude))*sin(rad(longtitude))
    z=((1-e**2)*N+ellHeight)*sin(rad(latitude))
    return  x,y,z


def xyz2ell(x,y,z):
    
    longtitude=atan2(y,x)

    p=(x**2+y**2)**(0.5)    
    
    N=a
    h=((x**2+y**2+z**2))**(0.5)-(a*b)**(0.5)
    latitude=atan2(z,((1-(e**2)*(N/(N+h)))*p))
    
    i=0
    while True:
        i=i+1
     
        eNi = N - a/((1-e**2*sin(latitude)**2)**(0.5))
        N=a/((1-e**2*sin(latitude)**2)**(0.5))
        
        ehi = h - ((p/(cos(latitude)))-N)
        h = (p/(cos(latitude)))-N
   
        elati = latitude - atan2(z,((1-(e**2)*(N/(N+h)))*p))
        latitude = atan2(z,((1-(e**2)*(N/(N+h)))*p))

        
        if abs(elati)<=10**(-8) and abs(eNi)<=10**(-4) and abs(ehi)<=10**(-4) :
            break

    ellHeight=h
    return degrees(latitude),degrees(longtitude),ellHeight


def degree2dms(decimal):
    d = int(decimal)
    m = int((decimal - d)*60)
    s = ((decimal - d)*60 - m)*60
    return "{} {} {:.4f}".format(d,m,s)

def dms2deg(dms):
    d, m, s      = list(map(float, dms.split(' ')))
    return d + m / 60 + s / 3600












