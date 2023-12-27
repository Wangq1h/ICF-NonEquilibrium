import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import scienceplots
import sys


# Define the constants
# All in SI units
A_alpha = 8e39 
m_e = 9.1*10**(-31) 
m_p= 1.637*10**(-27) 
m_i = 5*m_p
A_r = 8.95336*10**12 
A_e = 1.78455*10**(-10) 
A_omegaei = 1.90766*10**21 
kB = 1.381*10**(-23)
GammaB = 6600.87 
Rhoc = 1*10**6 
Cvi = 3/2*kB/m_i 
Cve = 3/2*kB/m_e 
cut_off = 2
hbar = 1.0545718*10**(-34)

# Define the functions
def Th(Te,Ti):
    return (Te+Ti)/2 

def SigmaV(Ti):
    return 9.1*10**(-16)*np.exp(-0.572*np.abs(np.log(Ti/(64.2*11604525.0062)))**(2.4))*10**(-6)

def Rhoe(rho):
    return 2*m_e*rho/(5*m_p+2*m_e) # g/cm^3

def Rhoi(rho):
    return 5*m_p*rho/(5*m_p+2*m_e) # g/cm^3

def Ne(rho):
    return Rhoe(rho)/m_e # 10^21 cm^-3

def Ni(rho):
    return Rhoi(rho)/m_i # 10^21 cm^-3

def Bmax(Te,rho):
    ne = Ne(rho)
    return (kB*Te/(4*np.pi*ne*(1.6e-19)**2))**(1/2)

def Bmin(Te,Ti):
    v = np.sqrt(3*Te*kB/m_e)
    wavelength = hbar/(2*m_e*v)
    r0 = (2*1.6e-19)**2/(3*kB*Ti)
    return np.max([wavelength,r0])

def LnLambda(Te,Ti,rho):
    bmin = Bmin(Te,Ti)
    bmax = Bmax(Te,rho)
    lnlambda = 0.5*np.log(1+(bmax/bmin)**2)
    return lnlambda

def Taualpha(Te,Ti,rho,Rh):
    th = Th(Te,Ti)
    lnLambda = LnLambda(Te,Ti,rho)
    return 45*lnLambda/5*rho*1000*Rh/(th/11604525.0062)**(3/2)

def f_alpha(Te,Ti,rho,Rh):
    tau = Taualpha(Te,Ti,rho,Rh)
    if tau<=1/2:
        return 3/2*tau-4/5*tau**2
    else:
        return 1-1/(4*tau)+1/(160*tau**3)

def W_alpha(Te,Ti,rho,Rh):
    sigmav = SigmaV(Ti)
    falpha = f_alpha(Te,Ti,rho,Rh)
    # print(sigmav)
    return A_alpha*Rhoi(rho)**2*sigmav*falpha

def f_alphai(Te):
    return Te/(Te+32*11604525.0062)

def f_alphae(Te):
    return 32*11604525.0062/(Te+32*11604525.0062)

def W_r(rho, Te):
    return A_r*rho**2*np.sqrt(Te)

def W_e(Te, Ti,rho, Rh):
    lnlambda = LnLambda(Te,Ti, rho)
    # th = Th(Te,Ti)
    return A_e*3*Te**(7/2)/(lnlambda*Rh**2)

def W_ie(Te, Ti, rho):
    lnlambda = LnLambda(Te,Ti,rho)
    return A_omegaei*lnlambda*rho**2*(-Te+Ti)/Te**(3/2)

def P_e(Te, rho):
    ne = Ne(rho)
    return ne*kB*Te

def P_i(Ti, rho):
    ni = Ni(rho)
    return ni*kB*Ti

def Uh(Te, Ti, rho):
    th = Th(Te,Ti)
    return (3/4*GammaB*th*rho/Rhoc)**(1/2)*4*10**(-15)

def W_mi(Te, Ti, rho, Rh):
    uh = Uh(Te, Ti, rho)
    pi = P_i(Ti, rho)
    return 3*pi*uh/Rh

def W_me(Te, Ti, rho, Rh):
    uh = Uh(Te, Ti, rho)
    pe = P_e(Te, rho)
    return 3*pe*uh/Rh

def Vh(Rh):
    return 4*np.pi/3*Rh**3


# Debug
def printall(Y,t):
    y = [Y[t,0],Y[t,1],Y[t,2],Y[t,3]]
    Te = y[0]
    Ti = y[1]
    rho = y[2]
    Rh = y[3]
    Walpha = W_alpha(y[0],y[1],y[2],y[3])
    falpha = f_alpha(y[0],y[1],y[2],y[3])
    Wme = W_me(y[0],y[1],y[2],y[3])
    Wmi = W_mi(y[0],y[1],y[2],y[3])
    We = W_e(y[0],y[1],y[2],y[3])
    Wr = W_r(y[2],y[0])
    Wie = W_ie(y[0],y[1],y[2])
    falphai = f_alphai(y[0])
    falphae = f_alphae(y[0])
    vh = Vh(y[3])
    uh = Uh(y[0],y[1],y[2])
    return [Walpha, Wie, Wr, We, Wme, Wmi]
    


# Define the differential equations
# y = [Te,Ti,rho,Rh]
def model(y, t):
    # dydt2 = 1/(Cvi*Rhoi(y[2]))*(W_alpha(y[0],y[1],y[2],y[3])*f_alphai(y[0])-W_ie(y[0],y[1],y[2])-W_mi(y[0],y[1],y[2],y[3]))
    # dydt1 = 1/(Cve*Rhoe(y[2]))*(W_alpha(y[0],y[1],y[2],y[3])*f_alphae(y[0])+W_ie(y[0],y[1],y[2])-W_me(y[0],y[1],y[2],y[3])-W_r(y[2],y[0])-W_e(y[0],y[1],y[2],y[3]))
    Walpha = W_alpha(y[0],y[1],y[2],y[3])
    falpha = f_alpha(y[0],y[1],y[2],y[3])
    Wme = W_me(y[0],y[1],y[2],y[3])
    Wmi = W_mi(y[0],y[1],y[2],y[3])
    We = W_e(y[0],y[1],y[2],y[3])
    Wr = W_r(y[2],y[0])
    Wie = W_ie(y[0],y[1],y[2])
    falphai = f_alphai(y[0])
    falphae = f_alphae(y[0])
    vh = Vh(y[3])
    uh = Uh(y[0],y[1],y[2])
    # print(Walpha, falpha, Wme, Wmi, We)#, Wr, Wie, falphai, falphae, vh, uh)
    dydt2 = 1/(Cvi*Rhoi(y[2]))*(Walpha*falphai-Wie-Wmi)
    dydt1 = 1/(Cve*Rhoe(y[2]))*(Walpha*falphae+Wie-Wme-Wr-We)
    # dydt2 = 1/(Cvi*y[2])*(Walpha*falphai-Wie-Wmi)
    # dydt1 = 1/(Cve*y[2])*(Walpha*falphae+Wie-Wme-Wr-We)
    dydt3 = (((Wr+We+Walpha*(1-falpha)/falpha)*vh)/(Cvi*y[1]+Cve*y[0])-4*np.pi*y[3]**2*y[2]*uh)/vh
    dydt4 = uh
    if np.isnan(y[0]):
        print("Te is nan")
        sys.exit()
    if np.isnan(y[1]):
        print("Ti is nan")
        sys.exit()
    if np.isnan(y[2]):
        print("rho is nan")
        sys.exit()
    if np.isnan(y[3]):
        print("Rh is nan")
        sys.exit()
    gamma = -dydt2/dydt1
    if t==0:
        print("gamma is", gamma)
    # print(dydt1, dydt2, dydt3, dydt4)
    return [dydt1, dydt2, dydt3, dydt4]


# Simulate the model
def Simulate(f,c, rho0=120*10**3, Rh0=30*10**(-6), Th0=8*11604525.0062, plot= True):
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange', 'purple']
    # Set up the initial conditions
    Th0 = 8*11604525.0062 
    Ti0 = Th0*f # keV
    Te0 = Th0*(2-f) # keV
    rho0 = 120*10**3 # g/cm^3
    Rh0 = 30*10**(-6) # cm rho0*rh0=3.6
    y0 = [Te0, Ti0, rho0, Rh0]
    t = np.linspace(0, 50*10**(-12), 10000)
    # Y0 = [Te0, Ti0, rho0, Rh0]
    # tt = t[1]-t[0]
    # Y =[]
    # Y.append(y0)
    # for time in t[0:-1]:
    #     y = Y[-1]
    #     dydt = model(y, tt)
    #     y2 = [y[0]+dydt[0]*tt, y[1]+dydt[1]*tt, y[2]+dydt[2]*tt, y[3]+dydt[3]*tt]
    #     Y.append(y2)
    #     # print(y2)

    # Ti = []
    # Te = []
    # for y in Y:
    #     Ti.append(y[1]/11604525.0062)
    #     Te.append(y[0]/11604525.0062)
    # Call the odeint function
    y = odeint(model, y0, t)
    # for y in Y:
    #     printall(y)
    # print(y[:, 0])
    # Print the solution
    Y = y.copy()
    for i in range(2):
        num = len(y[:, i])
        for j in range(num):
            y[j,i] = y[j,i]/11604525.0062
    t = [tt*10**12 for tt in t]   
    Walpha = []
    Wie = []
    Wr = []
    We = []
    Wme = []
    Wmi = []
    for tt in range(len(t)):
        info = printall(Y,tt)
        Walpha.append(np.abs(info[0]))
        Wie.append(np.abs(info[1]))
        Wr.append(np.abs(info[2]))
        We.append(np.abs(info[3]))
        Wme.append(np.abs(info[4]))
        Wmi.append(np.abs(info[5]))
        # print(Walpha)
    if not plot:
        if y[-1,0]>=y[0,0] and y[-1,1]>=y[0,1]:
            return True
    
    # 创建第一个子图，绘制 Walpha 和 Wie
    if plot:
        ax1.plot(t, Walpha, label='Walpha,f='+str(f),linestyle=c,color=colors[0])
        ax1.plot(t, Wie, label='Wie,f='+str(f),linestyle=c,color=colors[1])
        ax1.plot(t, Wr, label='Wr,f='+str(f),linestyle=c,color=colors[2])
        ax1.plot(t, We, label='We,f='+str(f),linestyle=c,color=colors[3])
        ax1.plot(t, Wme, label='Wme,f='+str(f),linestyle=c,color=colors[4])
        ax1.plot(t, Wmi, label='Wmi,f='+str(f),linestyle=c,color=colors[5])
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Walpha, Wie')
        ax1.set_yscale('log')

        # 创建第二个子图，共享 x 轴，绘制 Ti/Te
        
        # ax2.plot(t, y[:,0], label='Te,f='+str(f),linestyle=c,color=colors[6])
        # ax2.plot(t, y[:,1], label='Ti,f='+str(f),linestyle=c,color=colors[7])
        # ax2.plot(t, y[:,2], label='rho,f='+str(f),linestyle=c,color=colors[8])
        ax2.plot(t, y[:,3], label='Rh,f='+str(f),linestyle=c,color=colors[6])
        # ax2.set_ylabel('Ti/Te')
        # ax2.set_yscale('log')

    # 显示图例
    
    # plt.plot(t, Wr, label='Wr,f='+str(f),linestyle=':',color=c)
    # plt.plot(t, We, label='We,f='+str(f),linestyle='-.',color=c)
    # plt.plot(t, Wme, label='Wme,f='+str(f),linestyle='-',color=c)
    # plt.plot(t, Wmi, label='Wmi,f='+str(f),linestyle='--',color=c)
    # plt.xlabel('t (ps)')
    # plt.ylabel('W (J/m^3)')
    # plt.legend(loc='best', ncol=3)
    # plt.rc('legend', fontsize=6)
    # plt.plot(t, Te, label='Te,f='+str(f))
    # plt.plot(t, Ti, label='Ti,f='+str(f))
    # plt.plot(t, y[:,0], label='Te,f='+str(f),linestyle='--',color=c)
    # plt.plot(t, y[:,1], label='Ti,f='+str(f),linestyle='-',color=c)
    # plt.xlabel('t (ps)')
    # plt.ylabel('T (keV)')
    # plt.legend(loc='best', ncol=3)
    # plt.rc('legend', fontsize=6) 
    # plt.yscale('log')
    # plt.show()

# Plot the results
def Plot(f1=0.8,f2=1.0,f3=1.2):
    with plt.style.context('science'):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        Simulate(f1,'-')
        Simulate(f2,'--')
        Simulate(f3,':')
        fig.legend(loc='upper right')
        # if cut_off>0:
        #         plt.annotate( 'coulomb logarithm is cut off at '+str(cut_off),xy=(0, 1), xycoords='axes fraction',xytext=(0.03, 0.8),fontsize=7)
        # else:
        #         plt.annotate( 'coulomb logarithm is not cut off ',xy=(0, 1), xycoords='axes fraction',xytext=(0.03, 0.8),fontsize=7)
        plt.show()

def Scan(RhohRh, Th, f):
    Ignit_success = []
    for rhohrh in RhohRh:
        for th in Th:
            for rh in np.linspace(1/100*rhohrh, rhohrh, 100):
                if Simulate(f, 'k', rhohrh/rh, rh, th, plot=False):
                    Ignit_success.append([rhohrh, th])
                    break
            