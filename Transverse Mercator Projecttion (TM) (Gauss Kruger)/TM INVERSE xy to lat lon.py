from math import sin,cos,tan,radians,degrees,atan2
print("""========TRANSVERSE MERCATOR PROJECTION (TM) (GAUSS KRUGER)===========
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

EASTING=500000
NORTHING=7704493.0944
lonCM=30  

x=EASTING-500_000
y=NORTHING  # EĞER nokta southern(güney) yarım kürede ise 
            # 10 000 000 çıkar buradan

e=((a**2-b**2)/(a**2))**(0.5)
eı=((a**2-b**2)/(b**2))**(0.5)
print("e:",e,"\ne':",eı,"\n")

A0=1-(1/4)*e**2-(3/64)*e**4-(5/256)*e**6-(175/16384)*e**8
A2=(3/8)*(e**2+(1/4)*e**4+(15/128)*e**6-(455/4096)*e**8)
A4=(15/256)*(e**4+(3/4)*e**6-(77/128)*e**8)
A6=(35/3072)*(e**6-(41/32)*e**8)
A8=-(315/131072)*(e**8)

print("A0:",A0,"\nA2",A2,"\nA4",A4,"\nA6",A6,"\nA8",A8,"\n")


def f(lat):
    return a*(A0*lat-A2*sin(2*lat)+A4*sin(4*lat)-
              A6*sin(6*lat)+A8*sin(8*lat))-y

def fı(lat):
    return a*(A0-2*A2*cos(2*lat)+4*A4*cos(4*lat)-
              6*A6*cos(6*lat)+8*A8*cos(8*lat))

lat0=y/a
lat1=lat0-f(lat0)/fı(lat0)
epsilon=abs(lat1-lat0)
limit=10**(-12)
i=1
print("iteration:",i,"epsilon:",epsilon,"\tlatitude:",lat1)
while epsilon>=limit:
    i+=1
    lat2=lat1-f(lat1)/fı(lat1)
    epsilon=abs(lat2-lat1)
    print("iteration:",i,"epsilon:",epsilon,"\tlatitude:",lat2)
    lat1=lat2
    
    
print("\niteration number:",i,"\nepsilon:",epsilon)
print("\nfootpoint latitude:\n",lat1,"radians\n",degrees(lat1),"degree","\n",deg2dms(degrees(lat1)),"dms\n")

N=a/((1-e**2*sin(lat1)**2)**0.5)
M=(a*(1-e**2))/((1-e**2*sin(lat1)**2)**1.5)
t=tan(lat1)
n=eı*cos(lat1)
print("M:",M,"metre","\nN:",N,"metre","\nt:",t,"\nn:",n,"\n")


lat=(lat1-(t*x**2)/(2*M*N)+((t*x**4)/(24*M*N**3))*(5+3*t**2+n**2-4*n**4-9*n**2*t**2)\
     -((t*x**6)/(720*M*N**5))*(61-90*t**2+46*n**2+45*t**4-252*t**2*n**2-3*n**4\
      +100*n**6-66*t**2*n**4-90*t**4*n**2+88*n**8+225*t**4*n**4+84*t**2*n**6-192*t**2*n**8)\
      +((t*x**8)/(40320*M*N**7))*(1385+3633*t**2+4095*t**4+1575*t**6))


lon=(1/cos(lat1))*(x/N-(1/6)*((x/N)**3)*(1+2*t**2+n**2)\
     +(1/120)*((x/N)**5)*(5+6*n**2+28*t**2-3*n**4+8*t**2*n**2+24*t**4
      -4*n**6+4*t**2*n**4+24*t**2*n**6)
     -(1/5040)*(x/N)**7*(61+662*t**2+1320*t**4+720*t**6))
     
lon=lon+radians(lonCM)

print("LATITUDE:\n",lat,"radians\n",degrees(lat),"degree","\n",deg2dms(degrees(lat)),"dms\n")
print("LONGTITUDE:\n",lon,"radian\n",degrees(lon),"degree","\n",deg2dms(degrees(lon)),"dms\n")


"---------------MERIDIAN CONVERGENCE ANGLE-----------------------------"
gama=atan2(((t/N)*x-(t/3)*(x/N)**3*(1-n**2-2*n**4)+
          (t/15)*(x/N)**5*(2+2*n**2+9*n**4+6*t**2*n**2+20*n**6+
          3*t**2*n**4-27*t**2*n**6+11*n**8-24*t**2*n**8)-((17*t)/315)*(x/N)**7),1)
          
print("meridian convergence angle (gama):\n",gama,"radian\n",degrees(gama),"degree\n",deg2dms(degrees(gama)),"dms\n")



"------------------SCALE FACTOR------------------------------------"
k=1+((1+n**2)/2)*(x/N)**2+((1+6*n**2+9*n**4+4*n**6-24*t**2*n**4-
    24*t**2*n**6)/24)*(x/N)**4+(1/720)*(x/N)**6

print("scale factor(k):",k)




