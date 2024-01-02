# Non-equbllirum Hot Spot
This programme investigate the dynamic model of ICF, with seperated ion temperature and electron temperature. 
## 参考图
![Alt text](image-4.png)
## Lawson判据
来自Atezni:
$$
\begin{aligned}
\rho_{\mathbf{h}}R_{\mathbf{h}}T_{\mathbf{h}}& >\frac{9\sqrt{3}}{4}\frac{\Gamma_{\mathrm{B}}^{3/2}T_{\mathrm{h}}^{5/2}}{A_{\alpha}\langle\sigma v\rangle-A_{\mathrm{b}}T_{\mathrm{h}}^{1/2}}\left(\frac{\rho_{\mathrm{h}}}{\rho_{\mathrm{c}}}\right)^{1/2}  \\
&=\frac{1.1T_\mathbf{h}^{1/2}}{1-3.47T_\mathbf{h}^{-3/2}}\left(\frac{\rho_\mathbf{h}}{\rho_\mathbf{c}}\right)^{1/2}\text{g/cm}^2,
\end{aligned}
$$
```
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
```

来自张老师：
$$
\rho_{\mathbf{h}}R_{\mathbf{h}}T_{\mathbf{h}}>6\left(\frac{\rho_\mathbf{h}}{\rho_\mathbf{c}}\right)^{1/2}\text{g/cm}^2
$$
两者在$T_h=28$keV左右衔接。
## 采取的判据
升温：
```
if Y[-1,1]>Y[-2,1]:
    if Y[-1,1]>Y[0,1]:
        return True
    # return True
    return False
```
使用Lawson判据：
![Alt text](image.png)
在$0.7$g/cm$^2$处：
![Alt text](image-1.png)
![Alt text](image-2.png)

使用升温判据：
![Alt text](image-3.png)

有趣的现象，在某些初值，$f=1$情形出现分叉
![Alt text](8cb75127084fa2df5191053cd3a3247.png)

## 修改了比热
![Alt text](image-5.png)