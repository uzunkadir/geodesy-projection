from math import sin,cos,tan,radians,degrees,atan2
print("""========UNIVERSAL TRANSVERSE MERCATOR PROJECTION (UTM) =================
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

lonCM=27 #derece CENTRAL MERIDIAN
lat=dms2deg("57 43 19.23421") #derece LATITUDE fi
lon=dms2deg("27 30 30.2540") #derece LONGTITUDE lamda

k0=0.9996

e=((a**2-b**2)/(a**2))**(0.5)
eı=((a**2-b**2)/(b**2))**(0.5)
print("\ne:",e,"\ne':",eı,"\n")

lon=radians(lon) #convert to radian 
print("longtitude:",lon,"radian")

lonCM=radians(lonCM)
lon=lon-lonCM #DELTA lon
lat=radians(lat)

print("latitude:",lat,"radian","\nCentral Meridian (lon):",lonCM,"radian","\ndeltalon:",lon,"radian")


N=a/((1-e**2*sin(lat)**2)**(0.5))
t=tan(lat)
n=eı*cos(lat)

print("N:",N,"\nt",t,"\nn:",n,"\n")

A0=1-(1/4)*e**2-(3/64)*e**4-(5/256)*e**6-(175/16384)*e**8
A2=(3/8)*(e**2+(1/4)*e**4+(15/128)*e**6-(455/4096)*e**8)
A4=(15/256)*(e**4+(3/4)*e**6-(77/128)*e**8)
A6=(35/3072)*(e**6-(41/32)*e**8)
A8=-(315/131072)*(e**8)

print("A0:",A0,"\nA2",A2,"\nA4",A4,"\nA6",A6,"\nA8",A8,"\n")


Slat=a*(A0*lat-A2*sin(2*lat)+A4*sin(4*lat)-A6*sin(6*lat)+A8*sin(8*lat))
print("Meridian arc length (Sφ):\n",Slat,"metre\n")


x=(lon*cos(lat)+((lon**3*cos(lat)**3)/6)*(1-t**2+n**2)+
   ((lon**5*cos(lat)**5)/120)*(5-18*t**2+t**4+14*n**2-58*t**2*n**2+13*n**4+4*n**6-64*n**4*t**2-24*n**6*t**2)+
   ((lon**7*cos(lat)**7)/5040)*(61-479*t**2+179*t**4-t**6))*N

   
y=((Slat/N)+(lon**2*sin(lat)*cos(lat))/2+((lon**4*sin(lat)*cos(lat)**3)*(5-t**2+9*n**2+4*n**4))/24+
    ((lon**6*sin(lat)*cos(lat)**5)*(61-58*t**2+t**4+270*n**2-330*t**2*n**2+
    445*n**4+324*n**6-680*n**4*t**2+88*n**8-600*n**6*t**2-192*n**8*t**2))/720+
    ((lon**8*sin(lat)*cos(lat)**7)*(1385-311*t**2+543*t**4-t**6))/40320)*N
     
x=k0*x
y=k0*y

print("x:",x,"metre","\ny:",y,"metre")

x=x+500000

if y<0:
    y+=10_000_000
    
    

print("\nEASTING(x):",x,"metre","\nNORTHING(y):",y,"metre","\n")


"---------------MERIDIAN CONVERGENCE ANGLE terms of lat and lon----------------"
gama=(atan2((lon*sin(lat)*(1+(lon**2*cos(lat)*(1+t**2+3*n**2+2*n**4))/3+
                     ((lon**4*cos(lat)**4)/15)*(2+4*t**2+2*t**4+15*n**2+35*n**4-
                       40*t**2*n**4+33*n**6-60*t**2*n**6+18*n**8-24*t**2*n**8)+
                       (17/315)*(1+t**2)**3*lon**6*cos(lat)**6)),1))

print("meridian convergence angle:\n",gama,"radian\n",degrees(gama),"degree\n",deg2dms(degrees(gama)),"dms\n")



"------------------SCALE FACTOR   terms of lat and lon ------------------------------------"
k=1+((lon**2)/2)*cos(lat)**2*(1+n**2)+((lon**4*cos(lat)**4)/24)*(5-4*t**2+14*n**2+
    13*n**4-28*t**2*n**2+4*n**6-48*t**2*n**4-
    24*t**2*n**6)+((lon**6*cos(lat)**6)/720)*(61-148*t**2+16*t**4)
k=k*k0
print("Scale factor:",k)








