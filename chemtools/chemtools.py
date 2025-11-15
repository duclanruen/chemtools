# Module d'aide au calcul en formation Chimie
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Importation des modules nécessaire
#----------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit, Bounds
from scipy.optimize import fsolve

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Programme Analyse Point d'Equivalence

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Définition des fonctions nécessaires à l'exécution du programme
#----------------------------------------------------------------
# Tangente hyperbolique modifiée
def modified_tanh(x, a, b, x0, c, d):
    return a * np.tanh(b * (x - x0)) + c * x + d

# Chargement du fichier Excel
def charge_xls(filename,delim):
    file_path = filename
    df = pd.read_csv(file_path,delimiter=delim)
    print("Aperçu des données")
    print(df.head())
    nom_colonnes = df.columns
    x_data = df[nom_colonnes[0]].to_numpy()
    y_data = df[nom_colonnes[1]].to_numpy()

    return x_data, y_data

# Sélection de la zone de travail
def borne_selection(x,y):
    plt.plot(x, y, label="Courbe à analyser")
    plt.title("Cliquez sur les bornes inférieures et supérieures")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)

    print(f"Cliquez les 2 points des bornes dans la fenêtre graphique")
    points = plt.ginput(2)
    print(f"Points sélectionnés :", points)

    for point in points:
        plt.plot(point[0], point[1], 'ro')

    return points

# Bornes de définition de la zone d'analyse
def zone_reduction(x,y,points):
    xmin = int(np.min(np.where(x>points[0][0])))
    xmax = int(np.max(np.where(x<points[1][0])))
    xr = x[xmin:xmax]
    yr = y[xmin:xmax]

    return xr,yr 

# Calcul de la pente de sortie de la courbe 
def calcul_pente(x,y,xr,yr):
    x_max_ini = np.max(x)
    y_max_ini = np.max(y)
    x_max_red = np.max(xr)
    y_max_red = np.max(yr)

    pente = ((y_max_ini-y_max_red)/(x_max_ini-x_max_red))

    return pente

# Interpolation par la fonction tangente hyperbolique affine
# Calcul du point d'équivalence
def interpolation_tanhf(xr,yr,pente):  
    initial_guess = [np.mean(yr), 0.1, np.mean(yr), pente, 2]
    params, _ = curve_fit(modified_tanh, xr, yr, p0=initial_guess)
    a, b, x0, c, d = params
    #print(f"Paramètres ajustés : a={a}, b={b}, x0={x0}, c={c}, d={d}")
    print("-------------------------------")
    print(f"Point d'équivalence : {x0}")

    return params

# Représentation graphique des points interpolés
def calcul_valeur_interp(x,y,params):
    a, b, x0, c, d = params
    x_interpolated = np.linspace(np.min(x), np.max(x), 100)
    y_interpolated = modified_tanh(x_interpolated, a, b, x0, c, d)

    plt.scatter(x, y, label="Données originales", color="red")  
    plt.plot(x_interpolated, y_interpolated, label="Courbe ajustée", color="blue")  

    plt.legend()
    plt.xlabel("volume [ml]")
    plt.ylabel("pH")
    plt.title("Interpolation par une tangente hyperbolique modifiée")
    plt.show()

def coeffR(xr,yr,params):
    a, b, x0, c, d = params
    y_fit =  modified_tanh(xr, a, b, x0, c, d)
    ss_res = np.sum((y_fit-yr)**2)
    ss_tot = np.sum((yr-np.mean(yr))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"Coefficient de régression : {r_squared}")

    return r_squared

def EquiPoint(filename,delim):
    x_data, y_data = charge_xls(filename,delim)
    points = borne_selection(x_data,y_data)
    plt.close('all')
    x_reduit,y_reduit = zone_reduction(x_data,y_data,points)
    pente = calcul_pente(x_data,y_data,x_reduit,y_reduit)
    params = interpolation_tanhf(x_reduit,y_reduit,pente)
    calcul_valeur_interp(x_reduit,y_reduit,params)
    rs = coeffR(x_reduit,y_reduit,params)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Calcul pH en solution acqueuse

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Définition des fonctions nécessaires à l'exécution du programme
#----------------------------------------------------------------
def pH_Sol_Acq(A,B,pKa,pKe,C):
    Ke = 10**(-pKe)
    Ka = 10**(-pKa)

    def equations ( vars ) :
        x , y , h , w = vars
        eq1 = x+y-C
        eq2 = x+h-w-C
        eq3 = x * Ka-h * y
        eq4 = w* h-Ke
        return [eq1,eq2,eq3,eq4]

    x,y,h,w = fsolve(equations,(1,1,1,1))

    print('  ')
    print('Résulats')
    print('----------------------------------')
    print("C[{}] = {:.3e} [mol/L]".format(B,x))
    print("C[{}] = {:.3e} [mol/L]".format(A,y))    
    print("C[{}] = {:.3e} [mol/L]".format("H3O+",h))
    print("C[{}] = {:.3e} [mol/L]".format("OH-",w))
    print("Somme de concentration : {:.3e}".format(x+y+h+w))
    print('----------------------------------')
    print("pH = ",round(-np.log10(h),3))
    print('----------------------------------')

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

# Diagramme de distribution en fonction des concentration et du pKa

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
def diag_Dist_C(A,B,pKa):

    pHmax=14

    fig,ax = plt.subplots()
    pH = np.arange(0,pHmax,0.1)

    alpha0 = 100/(1+10**(pH-pKa))
    alpha1 = 100-alpha0

    a = plt.plot(pH,alpha0,label = A)
    b = plt.plot(pH,alpha1,label = B)
    plt.legend(loc = 'lower left')
    plt.xlabel('$pH$')
    plt.title("Diagramme de distribution de {}".format(A))

    plt.yticks(np.arange(0,100,10))

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(symbol= ' %' ))

    # réglage de la grille secondaire :
    minor_xticks = np.arange(0,14,1)    # grille secondaire en x
    minor_yticks = np.arange(0,101,10)  # grille secondaire en y
    ax.set_xticks(minor_xticks,minor = True)
    ax.set_yticks(minor_yticks,minor = True)
    plt.axis([0,pHmax,0,105])           # limites des axes
    plt.ylabel("proportion des éspèces (%) ")
    plt.xlabel("$pH$" )

    #Légende des courbes et points p a r t i c u l i e r s
    plt.text(2,100,r'$\alpha_o$' , backgroundcolor = 'white',fontsize = 12)
    plt.text(12,100,r'$\alpha_1$ ' , backgroundcolor= 'white',fontsize = 12)

    # Point d 'intersection des deux courbes
    plt.plot([0,10],[50,50],'g--')      # ligne alpha = 50%
    plt.plot([pKa,pKa],[0,70],'g--')     # ligne pH qui correspont à pKa

    plt.grid()
    plt.show()

