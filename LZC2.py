import numpy as np
import matplotlib.pyplot as plt
#start = cputime; #Tempo Inicial da Simulação

#======================= ENTRADA DE DADOS PRINCIPAL========================

Espessura = 1.91e-3; # Espessura da Fita [m]
Largura = 1.074; # Largura da Fita [m]
Velo = 1.08; # Velocidade da Fita [m/s]
VazaoGN = 210; #Vazao por Zona de queima direta de gás natural [Nm3/h]
es = 0.2; # Emissividade da Lâmina
ew = 0.9; # Emissividade da Parede do Forno
sigma = 5.678e-8; # Constante de Stefan-Boltzmann
rhoLam = 7854; # Densidade do Lâmina kg/m3 
MassLam = Largura * Espessura* rhoLam ; # strip width x thickness x density = strip mass/m


#========================= DIMENSOES DO FORNO =============================


# Comprimento [m] de cada zona do forno
ComprimentoZonas = [2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 6.7, 6.7, 6.7, 6.7, 6.7, 6.7, 6.7, \
                    1, 1, 1, 1, 1, 1, 1]; 

# Numero de sessões em cada zona, de forma que cada sesão tenha
# aproximadamente 1m
SessoesZona = [3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1,1]


NS = sum(SessoesZona) # Numero Total de sessões
FimZonaAque = sum(SessoesZona[0:13])+1 # Sessões até finalizar zona de aquecimento
print('NS= ', NS)

#==========================================================================
# Matrizes de Calor, Temperatura e Malha de Solução 

TL = np.zeros((NS+1,1)) # Temperatura da lâmina em x,t
TLt1 = np.zeros((NS+1,1)) # Temperatura da lâmina em x,t+1
TP = np.zeros((NS+1,5)) # Temperatura da parede em x,t
TPt1 = np.zeros((NS+1,5)) #Temperatura da parede em x,t+1
Q = np.zeros((NS+1,1)) # Calor transferido para a lâmina em cada posição
QParede = np.zeros((NS+1,1)) # Fluxo de calor na parede

print('TP= ', TP)
print('np.shape(TP)= ', np.shape(TP))

# Malha de Solução
MalhaY = np.zeros((NS+1,5))
# Espaçamento na zona de aquecimento
MalhaY[0:FimZonaAque,1] = 0.0002 
MalhaY[0:FimZonaAque,2] = 0.01
MalhaY[0:FimZonaAque,3] = 0.06
MalhaY[0:FimZonaAque,4] = 0.2
# Espacamento do final da zona de aquecimento
MalhaY[FimZonaAque:NS+1,1] = 0.001/1.4 
MalhaY[FimZonaAque:NS+1,2] = 0.01/1.4
MalhaY[FimZonaAque:NS+1,3] = 0.04/1.4
MalhaY[FimZonaAque:NS+1,4] = 0.1/1.4

print('MalhaY= ', MalhaY)
#Pontos de referencia [s,n,w,e sul norte oeste leste)
Ps = np.zeros((NS+1,5))
Pn = np.zeros((NS+1,5))
print('Ps= ', Ps)


for y in range(1,5):
    Ps[:,y] = MalhaY[:,y] - MalhaY[:,y-1]

print('Ps= ', Ps)

for y in range(0,4):
    Pn[:,y] = MalhaY[:,y+1] - MalhaY[:,y]


print('Pn= ', Pn)

Ps[:,0] = Ps[:,1]
Pn[:,4] = Pn[:,3]

print('Ps= ', Ps)

Pe = np.zeros((NS+1,1))
Pw = np.zeros((NS+1,1))

MalhaX = np.zeros((NS+1,1)) # Posição dos pontos na malha-y de solução
auxiliar = 1

for i in range(0,20):

    for j in range(0,SessoesZona[i]):
        MalhaX[auxiliar] = MalhaX[auxiliar-1] + ComprimentoZonas[i] / SessoesZona[i]
        auxiliar = auxiliar+1
    
print('MalhaX= ', MalhaX)    


Pe[0:NS] = MalhaX[1:NS+1] - MalhaX[0:NS]
print('MalhaX[1:NS] - MalhaX[0:NS-1] = ', MalhaX[1:NS] - MalhaX[0:NS-1])
Pw[1:NS+1] = MalhaX[1:NS+1] - MalhaX[0:NS]
Pe[NS] = Pe[NS-1]
Pw[0] = Pw[1]



#=============== Critério de Estabilidade (Numero de Courant) =============
T = 100000 # Tamanho dominio do tempo (s)
N = 20*T # Numero de divisões do tempo
k = T/N # Duração de cada passo de tempo

print('Pe= ',  Pe)   
print('Pw= ',  Pw)   


Cont = np.min(Pe[1:NS])
co = k*Velo/Cont
if (k*Velo/Cont) > 1:
    print('ERRO - Criterio de Estabilidade Não Alcançado')

#==========================================================================

#=================================================================================
# Estimativa Inicial Para Temperatura das Camadas de Parede em Estado Estacionario
#=================================================================================

TP[:,0] = 800+273
TP[:,1] = 700+273
TP[:,2] = 500+273
TP[:,3] = 200+273
TP[:,4] = 50+273


#==========================================================================

#============== CALOR (Potencia Fornecida Por Zona) =======================

#PARA ZONAS 1 a 6
ALIMGN = VazaoGN # Alimentacao GN [Nm3/h]
PCIGN = 8600 #PCI do GN [KCal/Nm3]
POTGNKCAL = PCIGN*ALIMGN #Potencia em Kcal/h
POTGNW = POTGNKCAL*1.163

GC079 = 433 # Alimentacao GCO para zonas 7 a 9 [Nm3/h]

#PARA ZONA 7
ALIMGCO = GC079/25*10 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL = PCIGCO*ALIMGCO #Potencia em Kcal/h
POTGCOW7 = POTGCOKCAL*1.163 #Potencia W

#PARA ZONA 8
ALIMGCO = GC079/25*8 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL = PCIGCO*ALIMGCO #Potencia em Kcal/h
POTGCOW8 = POTGCOKCAL*1.163 #Potencia W

#PARA ZONA 9
ALIMGCO = GC079/25*7 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL = PCIGCO*ALIMGCO #Potencia em Kcal/h
POTGCOW9 = POTGCOKCAL*1.163 #Potencia W

GCO1013 = 400 # Alimentacao GCO para zonas 10 a 13 [Nm3/h]

#PARA ZONA 10
ALIMGCO2 = GCO1013/24*8 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL2 = PCIGCO*ALIMGCO2 #Potencia em Kcal/h
POTGCOW10 = POTGCOKCAL2*1.163

#PARA ZONA 11
ALIMGCO2 = GCO1013/24*6 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL2 = PCIGCO*ALIMGCO2 #Potencia em Kcal/h
POTGCOW11 = POTGCOKCAL2*1.163

#PARA ZONA 12
ALIMGCO2 = GCO1013/24*6 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL2 = PCIGCO*ALIMGCO2 #Potencia em Kcal/h
POTGCOW12 = POTGCOKCAL2*1.163

#PARA ZONA 13
ALIMGCO2 = GCO1013/24*4 # Alimentacao GCO [m3/h]
PCIGCO = 4050 #PCI do GCO [KCal/m3]
POTGCOKCAL2 = PCIGCO*ALIMGCO2 #Potencia em Kcal/h
POTGCOW13 = POTGCOKCAL2*1.163


#============= Calor (Potencia) de Cada Zona ==============================

PotZona = [POTGNW, POTGNW, POTGNW, POTGNW, POTGNW, POTGNW, POTGCOW7, POTGCOW8, \
    POTGCOW9, POTGCOW10, POTGCOW11, POTGCOW12, POTGCOW13, 0, 0, 0, 0, 0, 0, 0]

#====================== Eficiencia de Cada Zona ===========================

EficienciaZona = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0, 0, 0, 0, 0, 0, 0])*0.6
print('EficienciaZona= ', EficienciaZona)

#================ Contabilização do Calor por Sessão ======================


Calor = np.zeros((NS+1,1)) # calor
auxiliar = 1
for i in range(0,19):
    print('i= ', i)
    print('SessoesZona[i-1]= ', SessoesZona[i-1])
    for j in range(0,SessoesZona[i-1]):
        print('j= ', j)
        Calor[auxiliar] = PotZona[i]/SessoesZona[i] * EficienciaZona[i]
        auxiliar = auxiliar+1
    
    

# Parte Constante da Equação de (Q - Calor transferido para a lâmina em
# cada posição)
QConstante = (2 * Largura * es * sigma) / (1 + es*(1-ew)*Largura/(ew*2.6))


#==========================================================================
#====================== CALULO DAS TEMPERATURAS ===========================
#==========================================================================
Tent = 303 # Temperatura de Entrada da Fita [K]
Tinf = 303 # Temperatura Ambiente (K)


TL[:] = Tent

print('Pe= ',Pe)

for j in range(0,1): 
    #Posição de x = 1
    TLt1[0] = Tent
    TPt1[0,0] = TP[1,0]
    TPt1[0,1] = TP[1,1]
    TPt1[0,2] = TP[1,2]
    TPt1[0,3] = TP[1,3]
    TPt1[0,4] = TP[1,4]
    

    for x in range(1,FimZonaAque+1): 
        
        Cont = Pw[x]   
        # Calculo da transferencia de calor por radiação entre a parede do forno e a fita     
        Q[x] = QConstante * (TP[x,0]**4 - TL[x]**4)*Cont

        print('====================================')
        print('x= ', x)
        print('Pw= ', Pw[x])
        print('Q= ', Q[x])

        cp_steel = 0.0003275*TL[x]*TL[x] + 0.101*TL[x] + 396.4
        TLt1[x] = (1-Velo*k/Cont) * TL[x] + (Velo*k/Cont) * TL[x-1] + (Q[x]*k / (MassLam * cp_steel * Cont))
        # Calculo da transferencia de calor por radiação para a parede do forno
        QParede[x] = (Calor[x] - Q[x])/(2*2.6*Cont);
        TPt1[x,1] = (0.0040/0.25)*QParede[x] + TP[x,2]
        


        for y in range(1, 3):
            TPt1[x,y] = TP[x,y] + ((0.25 * 2 * k)/(497*500))*(TP[x+1,y]/(Pe[x]*(Pe[x]+Pw[x]))  \
            + TP[x-1,y]/(Pw[x]*(Pe[x]+Pw[x]))  \
            + TP[x,y+1]/(Pn[x,y]*(Pn[x,y]+Ps[x,y]))  \
            + TP[x,y-1]/(Ps[x,y]*(Pn[x,y]+Ps[x,y]))  \
            - TP[x,y]*(1/(Pw[x]*Pe[x]) + 1/(Pn[x,y]*Ps[x,y])));
    
           
        TPt1[x,4] = (TP[x,3] + (0.2579*5*8.5*Cont/0.25)*Tinf)/(1 + 0.2579*5*8.5*Cont/0.25)
    
   

    #======================= FINAL ULTIMA ZONA x = M+1 ========================

    x = NS
    Cont = Pw[x]
    # Calculo da transferencia de calor por radiação entre a parede do forno e a fita
    Q[x] = QConstante*(TP[x,1]**4 - TL[x]**4)*Cont
    cp_steel = 0.0003275*TL[x]*TL[x] + 0.101*TL[x] + 396.4
    TLt1[x] = (1-Velo*k/Cont) * TL[x] + (Velo*k/Cont) * TL[x-1] + (Q[x]*k / (MassLam * cp_steel * Cont));
    TPt1[NS,0] = TPt1[NS-1,0]
    TPt1[NS,1] = TPt1[NS-1,1]
    TPt1[NS,2] = TPt1[NS-1,2]
    TPt1[NS,3] = TPt1[NS-1,3]
    TPt1[NS,4] = TPt1[NS-1,4]

    #==================== Ajuste Para Proxima Iteração ========================

    TL[:,0] = TLt1[:,0]
    TP[:,:] = TPt1[:,:]

    print('TL[:,0]= ', TL[:,0])


#==========================================================================
#================================ FIM =====================================
#==========================================================================

DistaForno = MalhaX[1:61,1]
TempCelcius = TL[1:61,1]-273.15

plt.subplot(1,1,1); # Plotar Grafico Temperatura da Lâmina

plt.xlabel('Posição X [m]')
plt.ylabel('Temperatura [º C]')
plt.plot(MalhaX[1:61,1],TL[1:61,1]-273,'b')
plt.xlim([0, 65])
plt.ylim([0, 1000])
plt.legend('Temperatura da Lâmina em ºC')


#finish = cputime; #Tempo Final da Simulação
#duration = finish - start #Tempo Total de Simulação