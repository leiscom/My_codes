
import numpy as np
import matplotlib.pyplot as plt


na = 32#numero de puntos angulares
bb = 0 
cc = 2*np.pi
m = 30
nr = 128#numero de puntos radiales
pmin = 10**-9#limite inferior radial
pmax = 1#limite superior radial
sep = np.log(10)/m #separacion para el cambio de variable
d = np.log(np.abs(1/(1-np.exp(-sep))))
z = np.log(np.abs(pmax/(pmin*(1-np.exp(-sep)))))
    
y = pmin

tole = 10**-3 # tolerancia para la convergencia
tole_error = 10**-6 # tolerancia para el error de la iteracion
mtsmax = 7 # numero de estados de Matsubara
    
    #np.polynomial.legendre.leggauss entrega entre -1 y 1 por lo tanto se 
    #efectua la correccion
    #pesos angulares     
xi, wi = np.polynomial.legendre.leggauss(na)
    
xi = -0.5*bb*(xi-1)+0.5*cc*(xi+1)
    
    #pesos radiales
ps, Rj = np.polynomial.legendre.leggauss(nr)
    
pi = -0.5*d*(ps-1)+0.5*z*(ps+1) 
    
pj = pmin*np.exp(pi)*(1-np.exp(-sep))
    
    #integrales angulares
angular = np.zeros((nr,nr,(mtsmax-(-mtsmax-1)+1),(mtsmax-(-mtsmax-1))+1),dtype=complex)
q=0




def mats(n,T):
    
    return (2*n+1)*np.pi*T-complex(0,1)*mu
def ff(k, p, x, m, n):
    
    return 1/(np.sqrt(k**2+p**2-2*k*p*np.cos(x)+(mats(n,T)-mats(m,T))**2))

def ker(mass):
    
    M = mass
    ker = np.zeros((nr, nr, (mtsmax-(-mtsmax-1)+1)),dtype=complex)
    kerint = np.zeros((nr, (mtsmax-(-mtsmax-1)+1)),dtype=complex)
    chsqrd = np.zeros((nr, (mtsmax-(-mtsmax-1)+1)),dtype=complex)
    s=0
    t=0
    
    for i in range(0,nr):
        for j in range(0,nr):
            for mf in range(0, (mtsmax-(-mtsmax-1)+1)):
                for nf in range(0, (mtsmax-(-mtsmax-1)+1)):
                    
                    ker[i][j][mf] = ker[i][j][mf]+alfa*(T/np.pi)*(((pj[i]*pj[i])*(M[i][mf])*angular[i][j][mf][nf])/(mats(nf+(-mtsmax-1),T)*mats(nf+(-mtsmax-1),T)+pj[i]*pj[i]+M[i][mf]*M[i][mf]))
    for j in range(0,nr):
        for mf in range(0, (mtsmax-(-mtsmax-1)+1)):
            for i in range(0, nr):
                
                kerint[j][mf] = kerint[j][mf]+ker[i][j][mf]*Rj[i]*0.5*(z-d)
    for j in range(0,nr):
        for mf in range(0, (mtsmax-(-mtsmax-1)+1)):   
            
            chsqrd[j][mf] = (kerint[j][mf]-M[j][mf])/(M[j][mf])
    error = np.amax(chsqrd)
    print(np.abs(error))
    return(error, kerint)




#listas de parametros que iremos iterando.
lista_alfa= [12,11,10,9,8,7,6,5,4,3,2,1]
lista_T   = [0.48,0.45,0.4]
lista_mu  = [0]


for j in range(0,len(lista_T)):
    valorT = str(lista_T[j])
    T = lista_T[j]
    for k in range(0,len(lista_alfa)):
        valoralfa = str(lista_alfa[k])
        alfa = lista_alfa[k]*(np.pi/8)
        intento=0
        for i in range(0,len(lista_mu)):
            valormu =str(lista_mu[i])
            mu = complex(0,lista_mu[i])
            if intento >= 20:
                continue
            errores = [0] #se da una lista de error, compara el error de la iteracion anterior con el actual.
    


    #comienza el loop de las integrales angulares
            for i in range(0,nr):
                for j in range(0,nr):
                    for mf in range(0, (mtsmax-(-mtsmax-1))+1):
                        for nf in range(0, (mtsmax-(-mtsmax-1))+1):
                            for l in range(0,na):
                    
                              angular[i][j][mf][nf] = angular[i][j][mf][nf]+complex((wi[l]*ff(pj[i],pj[j],xi[l],mf+(-mtsmax-1),nf+(-mtsmax-1)))*0.5*(cc-bb))
    #comienza el loop final
            M = np.zeros((nr,(mtsmax-(-mtsmax-1)+1)),dtype=complex)+complex(0.1)
            wildcard = 0
            it = 0
            while wildcard < 1:
                it += 1
                print(it)
                error, kerint = ker(M)
                errores.append(error)
                if np.abs(error) < tole:
                    wildcard += 1
                    Mfin = kerint
                elif np.abs(errores[it-1]-errores[it])<tole_error:
                    wildcard += 1
                    continue
                else:
                    M = kerint
            if np.abs(errores[it-1]-errores[it])<tole_error:
                continue
        
            cond = 0
            for i in range(0, nr):
                for mf in range(0, (mtsmax-(-mtsmax-1)+1)):
                    cond = cond + complex(((0.5*(z-d))*Rj[i])*(0.5*T/np.pi)*(((pj[i]**2)*M[i][mf])/((mats(mf + (-mtsmax-1),T))**2 + pj[i]*pj[i] + M[i][mf]*M[i][mf])))
            print("cond",cond)
    
    #seccion final
    #Se preparan los datos para entregar un txt con los datos de cond y Mm para cada T , alpfa y mu.
    #El orden es cond, Mm0, Mm1, Mm2, Mm3, Mm4, Mm5, Mm6, Mm7.
            Mm0 = np.zeros(nr,dtype=complex)
            Mm1 = np.zeros(nr,dtype=complex)
            Mm2 = np.zeros(nr,dtype=complex)
            Mm3 = np.zeros(nr,dtype=complex)
            Mm4 = np.zeros(nr,dtype=complex)
            Mm5 = np.zeros(nr,dtype=complex)
            Mm6 = np.zeros(nr,dtype=complex)
            Mm7 = np.zeros(nr,dtype=complex)
            
            for j in range(0,nr):
                Mm0[j] = Mfin[j][0]
                Mm1[j] = Mfin[j][1]
                Mm2[j] = Mfin[j][2]
                Mm3[j] = Mfin[j][3]
                Mm4[j] = Mfin[j][4]
                Mm5[j] = Mfin[j][5]
                Mm6[j] = Mfin[j][6]
                Mm7[j] = Mfin[j][7]
            file = open("/home/student01/resultados/T"+valorT+"alfa"+valoralfa+"mu"+valormu+"mu_negativo_imaginario_cond_corregido"+".txt","x")
            cond = repr(cond).replace("array","").replace("(","").replace(")","")
            Mm0 = repr(Mm0).replace("array","").replace("(","").replace(")","")
            Mm1 = repr(Mm1).replace("array","").replace("(","").replace(")","")
            Mm2 = repr(Mm2).replace("array","").replace("(","").replace(")","")
            Mm3 = repr(Mm3).replace("array","").replace("(","").replace(")","")
            Mm4 = repr(Mm4).replace("array","").replace("(","").replace(")","")
            Mm5 = repr(Mm5).replace("array","").replace("(","").replace(")","")
            Mm6 = repr(Mm6).replace("array","").replace("(","").replace(")","")
            Mm7 = repr(Mm7).replace("array","").replace("(","").replace(")","")
            file.write("cond-Mm0-Mm1-Mm2-Mm3-Mm4-Mm5-Mm6-Mm7"+"/"+cond+"/"+Mm0+"/"+Mm1+"/"+Mm2+"/"+Mm3+"/"+Mm4+"/"+Mm5+"/"+Mm6+"/"+Mm7)
            file.close()
            intento+=1

