import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.special as sp
import torch
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# 定义新的函数
def new_func(k, n):
    return sp.factorial(n) * 2**(-n) / (sp.factorial((n+k)/2) * sp.factorial((n-k)/2))

# 定义拟合函数
def fit_func(k, A, B):
    return 1/A * np.exp(-k**2 / B)

# 随机游走模拟
def random_walk(n):
    steps = np.random.choice([-1, 1], size=n)
    return np.sum(steps)

# 随机游走模拟 GPU 版本
def random_walk_gpu(n, num_simulations):
    # 创建一个在 GPU 上的随机步骤数组
    steps = torch.randint(0, 2, (num_simulations, n), device='cuda') * 2 - 1
    # 计算每次模拟的总步数
    return steps.sum(dim=1)

# 模拟次数
num_simulations = 1000000

# 步数
n = 256

# # 存储结果
# results = np.zeros(num_simulations)

# # 运行模拟
# for i in range(num_simulations):
#     results[i] = random_walk(n)

# 在 GPU 上运行模拟
results = random_walk_gpu(n, num_simulations)

# 将结果移动到 CPU 上，以便后续处理
results = results.cpu().numpy()


# # 计算概率分布
# unique, counts = np.unique(results, return_counts=True)
# probabilities = counts / num_simulations

# # 计算方差
# variance = np.var(results)

# # 打印方差
# print("Variance:", variance)

# # 绘制概率分布
# plt.bar(unique, probabilities)
# plt.xlabel('Final Position')
# plt.ylabel('Probability')
# plt.title('Probability Distribution of Final Positions')
# plt.show()

unique, counts = np.unique(results, return_counts=True)
probabilities = counts / num_simulations

# 拟合数据
# 提供初始参数猜测
initial_guess = [5.0, 10.0]
popt, pcov = curve_fit(fit_func, unique, probabilities, p0=initial_guess)
Popt = popt.copy()
Popt[0] = popt[0]/np.sqrt(2*np.pi)

# 打印拟合参数
print("A:", Popt[0])
print("B:", Popt[1])

# 绘制原始数据
plt.bar(unique, probabilities, label='Original Data')

# 绘制拟合曲线
k_values = np.linspace(min(unique), max(unique), 1000)
plt.plot(k_values, fit_func(k_values, *popt), 'r', label='Fit: A=%5.3f, B=%5.3f' % tuple(Popt))
plt.plot(k_values, new_func(k_values, n), 'g', label='Theoretical Function')

# 计算卡方值
residuals = probabilities - fit_func(unique, *popt)
chi_square = np.sum((residuals ** 2) / fit_func(unique, *popt))

print("Chi-square:", chi_square)

plt.xlabel('Final Position')
plt.ylabel('Probability')
plt.title('Probability Distribution of Final Positions')
plt.legend()
plt.show()