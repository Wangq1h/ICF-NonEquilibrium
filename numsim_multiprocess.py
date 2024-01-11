import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import scienceplots
import sys
# import multiprocessing
# from multiprocessing import Pool
# import threading
from tqdm import tqdm
import time
import matplotlib.ticker as ticker
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
# import warnings

# # 将运行时警告转换为错误
# warnings.filterwarnings('error', category=RuntimeWarning)

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
    if isobaric:
        return 0
    if not isobaric:
        th = Th(Te,Ti)
        return (3/4*GammaB*th*rho/Rhoc)**(1/2)*4*10**(-15)
    # th = Th(Te,Ti)
    # return (3/4*GammaB*th*rho/Rhoc)**(1/2)*4*10**(-15)

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
    dydt2 = 1/(2*Cvi*Rhoi(y[2]))*(Walpha*falphai-Wie-Wmi)
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
    # if t==0:
        # print("gamma is", gamma)
    # print(dydt1, dydt2, dydt3, dydt4)
    return [dydt1, dydt2, dydt3, dydt4]


# Simulate the model
def Simulate(ax1, ax2, f,c, rho0=120*10**3, Rh0=30*10**(-6), Th0=8*11604525.0062, plot= True):
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange', 'purple']
    # Set up the initial conditions
    # Th0 = 8*11604525.0062 
    Ti0 = Th0*f # keV
    Te0 = Th0*(2-f) # keV
    if not isobaric:
        rho0 = Rhoc
    # rho0 = 120*10**3 # g/cm^3
    # Rh0 = 30*10**(-6) # cm rho0*rh0=3.6
    y0 = [Te0, Ti0, rho0, Rh0]
    t = np.linspace(0, 50*10**(-12), 1000)
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
        return Is_ignite(Y)

    
    # 创建第一个子图，绘制 Walpha 和 Wie
    if plot:
        # print(max(max_derivative_change(Y[:,0]), max_derivative_change(Y[:,1]), max_derivative_change(Y[:,2]), max_derivative_change(Y[:,3])))
        # ax1.plot(t, Walpha, label='Walpha,f='+str(f),linestyle=c,color=colors[0])
        # # ax1.plot(t, Wie, label='Wie,f='+str(f),linestyle=c,color=colors[1])
        # # ax1.plot(t, Wr, label='Wr,f='+str(f),linestyle=c,color=colors[2])
        # # ax1.plot(t, We, label='We,f='+str(f),linestyle=c,color=colors[3])
        # # ax1.plot(t, Wme, label='Wme,f='+str(f),linestyle=c,color=colors[4])
        # # ax1.plot(t, Wmi, label='Wmi,f='+str(f),linestyle=c,color=colors[5])
        # ax1.set_xlabel('Time (ps)')
        # # ax1.set_ylabel('Walpha, Wie')
        # ax1.set_yscale('log')
        # # ax1.xaxis.set_major_locator(ticker.MaxNLocator(4))  # 限制x轴刻度的数量
        # # ax1.yaxis.set_major_locator(ticker.MaxNLocator(3))  # 限制y轴刻度的数量
        # # ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x)))  # 设置x轴刻度标签的格式
        # # ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x)))  # 设置y轴刻度标签的格式
        # ax1.yaxis.set_ticklabels([])  # 隐藏y轴刻度标签

        # 创建第二个子图，共享 x 轴，绘制 Ti/Te
        
        ax2.plot(t, y[:,0], label='Te,f='+str(f),linestyle=c,color=colors[6])
        ax2.plot(t, y[:,1], label='Ti,f='+str(f),linestyle=c,color=colors[7])
        # ax2.plot(t, y[:,2], label='rho,f='+str(f),linestyle=c,color=colors[8])
        # ax2.plot(t, y[:,3], label='Rh,f='+str(f),linestyle=c,color=colors[6])
        ax2.set_ylabel('Ti/Te')
        # ax2.set_yscale('log')
        # ax2.xaxis.set_major_locator(ticker.MaxNLocator(4))  # 限制x轴刻度的数量
        # ax2.yaxis.set_major_locator(ticker.MaxNLocator(3))  # 限制y轴刻度的数量
        # ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x)))  # 设置x轴刻度标签的格式
        # ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x)))  # 设置y轴刻度标签的格式

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

#Lawson Criterion
def Lawson(th, rhoh, rh):
    # print(th, rhoh, rh)
    lhs = rhoh/(10**3)*rh*100*th/11604525.0062
    if lhs>10:
        rhs = (1.1*np.sqrt(th/11604525.0062))/(1-3.47*(th/11604525.0062)**(-3/2))*np.sqrt(rhoh/Rhoc)
    else:
        rhs = 6*np.sqrt(rhoh/Rhoc)
    # print(lhs, rhs)
    # input("Press Enter to continue...")
    # rhs = 6*np.sqrt(rhoh/Rhoc)
    if lhs>rhs:
        return True
    return False

def Precheck(th, rhoh, rh):
    lhs = rhoh/(10**3)*rh*100*th/11604525.0062
    if th>28*11604525.0062:
        rhs = (1.1*np.sqrt(th/11604525.0062))/(1-3.47*(th/11604525.0062)**(-3/2))*np.sqrt(rhoh/Rhoc)
    else:
        rhs = 6*np.sqrt(rhoh/Rhoc)
    if lhs>1/5*rhs:
        return True
    return False

def Is_smooth(points, threshold=1):
    # 计算点的导数
    derivatives = np.gradient(points)
    A = np.trapz(derivatives)
    derivatives = derivatives/A
    # 计算导数的变化
    changes = np.abs(np.diff(derivatives))
    # print(max(changes))
    # 如果导数的任何连续变化大于阈值，则认为曲线不平滑
    if np.isnan(changes).any():
        return False
    else:
        if np.any(changes > threshold):
            return False
        else:
            return True

def max_derivative_change(points):
    # 计算点的导数
    derivatives = np.gradient(points)
    # 归一化
    A = np.trapz(derivatives)
    derivatives = derivatives/A

    # 计算导数的变化
    changes = np.abs(np.diff(derivatives))
    
    # 返回导数变化的最大值
    return np.max(changes)

def Is_ignite(Y):
    Temp_up = False
    Smooth = True
    for i in range(2):
        if not Is_smooth(Y[:,i]):
            Smooth = False
        if Y[-1,i]>5*11604525.0062:
            if Y[-1,i]>Y[-2,i] or np.abs(Y[-1,i]-Y[-2,i])/Y[-1,i]<0.01:
                Temp_up = True
    # print(Temp_up, Smooth)
    if Temp_up and Smooth:
        delta_T = [np.abs(Y[i,0]-Y[i,1]) for i in range(len(Y[:,0]))]
        min_value = min(delta_T)
        equibllrium = delta_T.index(min_value)
        sample_Te = [Y[i+equibllrium,0] for i in range(len(Y[:,0])-equibllrium)]
        sample_Ti = [Y[i+equibllrium,1] for i in range(len(Y[:,1])-equibllrium)]
        if len(sample_Te)<5 or len(sample_Ti)<5:
            return False
        first_derivative_e = np.gradient(sample_Te)
        second_derivative_e = np.gradient(first_derivative_e)
        first_derivative_i = np.gradient(sample_Ti)
        second_derivative_i = np.gradient(first_derivative_i)
        count = 0
        for value in second_derivative_i:
            if value > 0:
                count += 1
                if count >= 5:
                    return True
            else:
                count = 0
        return False
    else:
        return False
    

# Plot the results
def Plot(f1=0.8,f2=1.0,f3=1.2):
    with plt.style.context('science'):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        Simulate(ax1, ax2, f1,'-')
        Simulate(ax1, ax2, f2,'--')
        Simulate(ax1, ax2, f3,':')
        fig.legend(loc='upper right')
        # if cut_off>0:
        #         plt.annotate( 'coulomb logarithm is cut off at '+str(cut_off),xy=(0, 1), xycoords='axes fraction',xytext=(0.03, 0.8),fontsize=7)
        # else:
        #         plt.annotate( 'coulomb logarithm is not cut off ',xy=(0, 1), xycoords='axes fraction',xytext=(0.03, 0.8),fontsize=7)
        plt.show()

def simulate_wrapper(args):
    try:
        ax1, ax2, f, rhohrh, rh, th = args
        return Simulate(ax1, ax2, f, 'k', rhohrh/rh, rh, th, plot=False)
    except Exception as e:
        print(f"Error in process: {e}")
        return None

# def Scan(RhohRh, Th, f):
#     ax1 = 1
#     ax2 = 2
#     Ignit_success = []

#     total = len(RhohRh)
#     pbar = tqdm(total=total, desc="Processing", ncols=80)

#     with ProcessPoolExecutor() as executor:
#         for th in Th:
#             futures = []
#             for rhohrh in RhohRh:
#                 pbar.update(1)
#                 if isobaric:
#                     for rh in np.linspace(10**(-5), 10**(-4), 10):
#                         if rhohrh/10*th/11604525.0062>0.5:
#                             future = executor.submit(simulate_wrapper, (ax1, ax2, f, rhohrh, rh, th))
#                             futures.append(future)
#                 else:
#                     if rhohrh/10*th/11604525.0062>0.5:
#                         future = executor.submit(simulate_wrapper, (ax1, ax2, f, Rhoc, rhohrh/Rhoc, th))
#                         futures.append(future)
#                 for future in as_completed(futures):
#                     result = future.result()
#                     if result is not None and result:
#                         Ignit_success.append([rhohrh/10, th/11604525.0062])
#                         break
#             else:
#                 continue
#             break
#         print("Simulation finished")
#         print("Waiting for all processes to finish...")

#     plt.plot([i[0] for i in Ignit_success], [i[1] for i in Ignit_success], 'o')
#     plt.show()

def Scan(RhohRh, Th, f):
    ax1 = 1
    ax2 = 2
    Ignit_success = []

    total = len(RhohRh)
    pbar = tqdm(total=total, desc="Processing", ncols=80)

    with ProcessPoolExecutor() as executor:
        for rhohrh in RhohRh:
            pbar.update(1)
            futures = []
            is_find = False
            # num_process = 0
            for th in Th:
                # num_process += 1
                if isobaric:
                    for rh in np.linspace(10**(-6), 10**(-4), 20):
                        if rhohrh/10*th/11604525.0062>0.5:
                            future = executor.submit(simulate_wrapper, (ax1, ax2, f, rhohrh, rh, th))
                            futures.append(future)
                else:
                    if rhohrh/10*th/11604525.0062>0.5:
                        future = executor.submit(simulate_wrapper, (ax1, ax2, f, rhohrh, rhohrh/Rhoc, th))
                        futures.append(future)
                # if num_process%10==0:
                for future in as_completed(futures):
                    if not future.cancelled():
                        result = future.result()
                        if result is not None and result:
                            Ignit_success.append([rhohrh/10, th/11604525.0062])
                            for future in futures:
                                future.cancel()
                            is_find = True
                            break
                if is_find:
                    break
        print("Simulation finished")
    with open(str(isobaric)+'Ignit_success_f='+str(f)+'.txt', 'w') as file:
        for item in Ignit_success:
            file.write(f'{item[0]}, {item[1]}\n')
    print("Waiting for all processes to finish...")

    plt.plot([i[0] for i in Ignit_success], [i[1] for i in Ignit_success], 'o')
    plt.show()

def Insight(th, rhohrh,f):
    fig, axs = plt.subplots(2, 5, figsize=(15, 6))  # 创建一个2行5列的子图网格
    axs = axs.flatten()  # 将子图网格转换为一维数组，以便我们可以在循环中使用它
    total = 10
    pbar = tqdm(total=total, desc="Processing", ncols=80)
    for i, rh in enumerate(np.linspace(10**(-6), 10**(-4), total)):
        pbar.update(1)
        ax1 = axs[i]
        ax2 = ax1.twinx()
        if isobaric:
            Simulate(ax1, ax2,f,'-', rhohrh/rh, rh, th, plot=True)
        else:
            Simulate(ax1, ax2,f,'-', Rhoc, rhohrh/Rhoc, th, plot=True)
        ax1.set_title(f'rho: {rhohrh/rh/1000:.1f}g/cm$^3$, rh: {rh*1e6:.1f}$\mu$ m\n th: {th/11604525.0062:.1f}keV, rhohrh: {rhohrh/10:.1f}g/cm$^2$')  # 设置子图的标题
    fig.legend(bbox_to_anchor=(0.5, -0.05), loc='upper center', ncol=5)  # 创建图例
    plt.tight_layout()  # 调整子图的位置，以确保它们不会重叠
    plt.show()
            

isobaric = True
# Plot()
# np.seterr(all='ignore')
def main():
    # 这里是你的主程序
    global isobaric
    isobaric = True
    Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 100), 0.8)
    # Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 200), 1)
    # Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 200), 1.2)
    isobaric = False
    # Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 200), 0.8)
    # Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 200), 1)
    # Scan(np.linspace(0.01*10**(1), 1.5*10**(1), 20), np.linspace(1*11604525.0062, 20*11604525.0062, 200), 1.2)

if __name__ == '__main__':
    main()