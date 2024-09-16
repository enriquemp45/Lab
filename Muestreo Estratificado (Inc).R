##############################################################################
###                                                                        ###
###    UNIDAD 3. MUESTREO ESTRATIFICADO                                    ###
###    FÓRMULAS BASADAS EN EL LIBRO "ELEMENTOS DE MUESTREO" DE MENDENHALL  ###
###                                                                        ###
##############################################################################


###########################################################################
"PROGRAMA #1. ESTIMACIÓN DE LA MEDIA POBLACIONAL EN MUESTREO ESTRATIFICADO"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(55,80,65)

"Introduce el vector de los tamaños de muestra ni por estrato"
n=c(14,20,16)

"Introduce el vector de las medias muestrales por estrato"
media=c(79.71,64.75,37.44)

"Introduce el vector de las varianzas muestrales por estrato"
s=c(105.14,158.20,186.13)

mstr.m=function(N,n,media,s)
{ 
Ntot=sum(N)
f=n/N
mestr=crossprod(N,media)/Ntot
varm=(s/n)*(1-f)
vstr=crossprod(N^2,varm)/Ntot^2
setr=sqrt(vstr)
limite=2*setr
a1=mestr-2*setr
b1=mestr+2*setr
cat("Media_Estr =",mestr,   "  Limite de error =",limite,   "  Intervalo =","(",a1, ",",b1,")","\n")
}
mstr.m(N,n,media,s)


#########################################################################
"PROGRAMA #2. ESTIMACIÓN DEL TOTAL POBLACIONAL EN MUESTREO ESTRATIFICADO"

Se deja como tareíta


################################################################################
"PROGRAMA #3. ESTIMACIÓN DE LA PROPORCIÓN POBLACIONAL EN MUESTREO ESTRATIFICADO"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(97,43,145,68)

"Introduce el vector de los tamaños de muestra ni por estrato"
n=c(39,17,69,33)

"Introduce el vector de las proporciones muestrales pi por estrato"
p=c(0.87,0.93,0.60,0.53)

mstr.p=function(N,n,p)
{ 
Ntot=sum(N)
f=n/N
q=1-p
pestr=crossprod(N,p)/Ntot
varp=(p*q/(n-1))*(1-f)
vstr=crossprod(N^2,varp)/Ntot^2
setr=sqrt(vstr)
limite=2*setr
a1=pestr-2*setr
b1=pestr+2*setr
cat("Proporción_Estr =",pestr,   "  Limite de error =",limite,   "  Intervalo =","(",a1, ",",b1,")","\n")
}
mstr.p(N,n,p)


################################################################################
"PROGRAMA #4. TAMAÑO DE MUESTRA PARA ESTIMAR LA MEDIA POBLACIONAL ESTRATIFICADA"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(155,62,93)

"Introduce el límite de error de estimación, B"
B=2

"Introduce el vector de las varianzas poblacionales por estrato"
s=c(25,225,100)

"Introduce el vector de los pesos wi por estrato"
w=c(1/3,1/3,1/3)

msa.nms=function(N,B,s,w)
{ 
D=(B^2)/4
Ntot=sum(N)
a=crossprod((N^2)*s,1/w)
b=crossprod(N,s)
n=a/(D*Ntot^2+b)
n1=ceiling(n)
cat("n =",n1)
}
msa.nms(N,B,s,w)


################################################################################
"PROGRAMA #5. TAMAÑO DE MUESTRA PARA ESTIMAR EL TOTAL POBLACIONAL ESTRATIFICADO"

Se deja como Taréta


#####################################################################################
"PROGRAMA #6. TAMAÑO DE MUESTRA PARA ESTIMAR LA PROPORCION POBLACIONAL ESTRATIFICADA"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(3000,5000,2000)

"Introduce el límite de error de estimación, B"
B=0.05

"Introduce el vector de las proporciones poblacionales pi por estrato"
p=c(0.3,0.5,0.2)

"Introduce el vector de los pesos wi por estrato"
w=c(1/3,1/3,1/3)

msa.nps=function(N,B,p,w)
{ 
Ntot=sum(N)
D=B^2/4
a=crossprod((N^2)*p*(1-p),1/w)
b=crossprod(N*p,1-p)
n=a/(D*Ntot^2+b)
n1=ceiling(n)
cat("n =",n1)
}
msa.nps(N,B,p,w)


################################################################################
"PROGRAMA #7. ASIGNACIÓN ÓPTIMA PARA ESTIMAR LA MEDIA POBLACIONAL ESTRATIFICADA"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(112,68,39)

"Introduce el límite de error de estimación, B"
B=0.6325

"Introduce el vector de las desviaciones poblacionales por estrato"
s=c(1.5,1.8,1.8)

"Introduce el vector de los costos unitarios ci por estrato"
cos=c(9,25,36)

msa.nmso=function(N,B,s,cos)
{ 
D=(B^2)/4
Ntot=sum(N)
a=crossprod(N*s,1/sqrt(cos))
b=crossprod(N*s,sqrt(cos))
n=(a*b)/(D*Ntot^2+crossprod(N,s^2))
n1=ceiling(n)
na=n1*(N*s*(1/sqrt(cos))/a)
cat("n =",n1, "   Tamaños por estrato =", na,"\n")
}
msa.nmso(N,B,s,cos)


################################################################################
"PROGRAMA #8. ASIGNACIÓN ÓPTIMA PARA ESTIMAR EL TOTAL POBLACIONAL ESTRATIFICADO"

Se deja como Taréita


###################################################################################
"PROGRAMA #9. ASIGNACIÓN DE NEYMAN PARA ESTIMAR LA MEDIA POBLACIONAL ESTRATIFICADA"

Se deja como Tareíta


####################################################################################
"PROGRAMA #10. ASIGNACIÓN DE NEYMAN PARA ESTIMAR EL TOTAL POBLACIONAL ESTRATIFICADO"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(112,68,39)

"Introduce el límite de error de estimación, B"
B=0.6325

"Introduce el vector de las desviaciones poblacionales por estrato"
s=c(1.5,1.8,1.8)

msa.ntn=function(N,B,s)
{ 
Ntot=sum(N)
D=(B^2)/(4*Ntot^2)
a=crossprod(N,s)
b=crossprod(N,s^2)
n=a^2/(D*Ntot^2+b)
n1=ceiling(n)
na=n1*((N*s)/a)
cat("n =",n1, "   Tamaños por estrato =", na,"\n")
}
msa.ntn(N,B,s)


#######################################################################################
"PROGRAMA #11. ASIGNACIÓN PROPORCIONAL PARA ESTIMAR LA MEDIA POBLACIONAL ESTRATIFICADA"

Se deja como Tareíta


#######################################################################################
"PROGRAMA #12. ASIGNACIÓN PROPORCIONAL PARA ESTIMAR EL TOTAL POBLACIONAL ESTRATIFICADO"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(112,68,39)

"Introduce el límite de error de estimación, B"
B=0.6325

"Introduce el vector de las varianzas poblacionales por estrato"
var=c(2.25,3.24,3.24)

msa.ntp=function(N,B,var)
{ 
Ntot=sum(N)
D=(B^2)/(4*Ntot^2)
s=var
a=crossprod(N,s)
n=a/(D*Ntot+a/Ntot)
n1=ceiling(n)
na=n1*(N/Ntot)
cat("n =",n1, "   Tamaños por estrato =", na,"\n")
}
msa.ntp(N,B,var)


######################################################################################
"PROGRAMA #13. ASIGNACIÓN ÓPTIMA PARA ESTIMAR LA PROPORCION POBLACIONAL ESTRATIFICADA"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(97,43,145,68)

"Introduce el límite de error de estimación, B"
B=0.05

"Introduce el vector de las proporciones poblacionales pi, por estrato"
p=c(0.9,0.9,0.5,0.5)

"Introduce el vector de los costos unitarios ci por estrato"
cos=c(4,4,8,8)

msa.npso=function(N,B,p,cos)
{ 
D=(B^2)/4
Ntot=sum(N)
q=1-p
a=crossprod(N*sqrt(p),sqrt(q)*(1/sqrt(cos)))
b=crossprod(N*sqrt(p),sqrt(q)*sqrt(cos))
n=(a*b)/(D*Ntot^2+crossprod(N*p,q))
n1=ceiling(n)
na=n1*(N*sqrt(p)*sqrt(q)*(1/sqrt(cos))/a)
cat("n =",n1, "   Tamaños por estrato =", na,"\n")
}
msa.npso(N,B,p,cos)


#########################################################################################
"PROGRAMA #14. ASIGNACIÓN DE NEYMAN PARA ESTIMAR LA PROPORCION POBLACIONAL ESTRATIFICADA"

Se deja como Tareíta


############################################################################################
"PROGRAMA #15. ASIGNACIÓN PROPORCIONAL PARA ESTIMAR LA PROPORCION POBLACIONAL ESTRATIFICADA"

"DATOS"
"Introduce el vector de los tamaños de población Ni por estrato"
N=c(97,43,145,68)

"Introduce el límite de error de estimación, B"
B=0.05

"Introduce el vector de las proporciones poblacionales pi, por estrato"
p=c(0.9,0.9,0.5,0.5)

msa.npp=function(N,B,p)
{ 
Ntot=sum(N)
D=(B^2)/4
q=1-p
a=crossprod(N*p,q)
n=a/(D*Ntot+a/Ntot)
n1=ceiling(n)
na=n1*(N/Ntot)
cat("n =",n1, "   Tamaños por estrato =", na,"\n")
}
msa.npp(N,B,p)

############################################################################################
