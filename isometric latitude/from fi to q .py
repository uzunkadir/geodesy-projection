from math import radians,degrees,log,sin,tan,pi

# EVEREST Ellipsoid parameters
a=6377276.345
b=6356075.41314024
e=((a*a-b*b)/(a*a))**(1/2) 


# GRS80 Ellipsoid parameters
a=6378137
b=6356752.314140347
e=((a**2-b**2)/(a**2))**(1/2) 


filist=[1,10,11,12,20,30,40,50,60,70,75,80,85,86,87,88,89]

print("\t\tEVEREST ELLIPSOID\n","Geodetic Latitude(fi)\t","İsometric Latitude(q)\n")

for fi in filist:
    
    fi=radians(fi)
    q=log(tan(pi/4+fi/2)*((1-e*sin(fi))/(1+e*sin(fi)))**(e/2))
    q=degrees(q)
    fi=degrees(fi)
    print("\t",format(fi,".0f"),"------------>",q,"\n")
    
