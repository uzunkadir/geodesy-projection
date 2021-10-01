from math import sin,cos,tan,radians,degrees,atan,log,pi
print("""================MERCATOR PROJECTION=================
      INVERSE PROBLEM
      x,y ------> Latitude,Longtitude\n""")

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

x=3228265.2330049337
y=3228918.578910121


e=((a**2-b**2)/(a**2))**(0.5)
eı=((a**2-b**2)/(b**2))**(0.5)
print("e:",e,"\ne':",eı,"\n")


lon=x/a
q=y/a

print("isometric latitude (q):\n",q,"radian\n",degrees(q),"degree\n",deg2dms(degrees(q)),"dms\n")

iteration=0 
lat=2*atan(e**q)-pi/2

while True:
    iteration+=1
    latfunction=(1/2)*(log(1+sin(lat))-log(1-sin(lat))+
                e*log(1-e*sin(lat))-e*log(1+e*sin(lat)))-q
    latderivative=(1-e**2)/((1-e*e*sin(lat)*sin(lat))*cos(lat))
    
    
    lat1=lat-latfunction/latderivative 
    epsilon=abs(lat-lat1)
    print("iteration number:",iteration,"Epsilon:",epsilon,"\tlatitude:",lat1)
    lat=lat1 #update lat
        
    if (epsilon<=10**(-8)): # Control iteration step
        break

print("\niteration number:",iteration,"\nepsilon:",epsilon,"\n")

print("LATITUDE:\n",lat,"radian\n",degrees(lat),"degree\n",deg2dms(degrees(lat)),"dms\n")
print("LONGTITUDE:\n",lon,"radian\n",degrees(lon),"degree\n",deg2dms(degrees(lon)),"dms\n")




    
