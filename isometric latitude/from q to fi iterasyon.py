from math import radians,atan,e,pi,log,sin,cos,degrees

#GRS80 Ellipsoid parameters
a=6378137
b=6356752.314140347
ek=((a**2-b**2)/(a**2))**(1/2) #eksentrisite

# WGS84 Ellipsoid Parameters
a=6378137
b=6356752.3142
ek=((a**2-b**2)/(a**2))**(1/2) #eksentrisite

#isometric latitude (q) 
qlist=[0.9934132219,43.4668126053,179.0306399306] # isometric latitude data

print("\t\t GRS80 ELLIPSOID\n","Geodetic Latitude(fi)\t\t","İsometric Latitude(q)\n")

for qd in qlist:
    iteration=0 
    q=radians(qd)
    fi=2*atan(e**q)-pi/2
    while True:
        iteration+=1
        fifunction=(1/2)*(log(1+sin(fi))-log(1-sin(fi))+
                    ek*log(1-ek*sin(fi))-ek*log(1+ek*sin(fi)))-q
        fiderivative=(1-ek**2)/((1-ek*ek*sin(fi)*sin(fi))*cos(fi))
        
        
        fi1=fi-fifunction/fiderivative 
        epsilon=abs(fi-fi1)
        
        fi=fi-fifunction/fiderivative #update fi
        
        if (epsilon<=10**(-8)): # Control iteration step
            print(degrees(fi1),"<---------------",degrees(q),)
            print("iteration number:",iteration)
            print("epsilon:",epsilon,"\n")
            break
