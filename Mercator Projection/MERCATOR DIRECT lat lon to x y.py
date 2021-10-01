from math import sin,cos,tan,radians,log,pi,atan
print("""================MERCATOR PROJECTION=================
      DIRECT PROBLEM
      Latitude,Longtitude ------> x,y""")

def deg2dms(deg):
    d            = int(deg)
    m            = int((deg - d)*60)
    s            = (deg - d - m/60)*3600
    return "{:d} {:02d} {:03f}".format(d, m, s)

def dms2deg(dms):
    d, m, s      = list(map(float, dms.split(' ')))
    return d + m / 60 + s / 3600


#GRS80
a=6378137
b=6356752.3141
f=(a-b)/a
print(f)
lat=28 #derece LATITUDE fi
lon=29 #derece LONGTITUDE lamda

e=((a**2-b**2)/(a**2))**(0.5)
eı=((a**2-b**2)/(b**2))**(0.5)
print("\ne:",e,"\ne':",eı,"\n")


lat=radians(lat)
lon=radians(lon)
q=log(tan(pi/4+lat/2)*((1-e*sin(lat))/(1+e*sin(lat)))**(e/2))

x=a*lon    
y=a*q
print("x:",x,"metre","\ny:",y,"metre\n")



"------------------SCALE FACTOR   terms of lat and lon------------------------"
N=a/((1-e**2*sin(lat)**2)**(0.5))

k=a/(N*cos(lat))
print("scale factor:",k)

"---------------MERIDIAN CONVERGENCE ANGLE terms of lat and lon---------------"
print("meridian convergence angle(gama) is always zero at Mercator Projection")




