import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 设置随机种子以保证结果可重复
np.random.seed(42)

# 1. 设定物理参数 (保证理论恒等式 K = Omega + A * exp(-kappa) 严格成立)
Omega_true = 6000.0
A_true = 10000.0
kappa_true = 0.223
beta_true = 0.26

# 精确计算 k=1 处的理论交点基数 K
K_true = Omega_true + A_true * np.exp(-kappa_true * 1.0)  # 约 13998.6

# 2. 模拟样本梯度 k (1 到 34 个基因组)
k_data = np.arange(1, 35)
num_replicates = 15  # 每个梯度 15 次重抽样

raw_k = []
raw_V_pan = []
raw_V_core = []

# 3. 模拟数据（核心修改：k=1 处共享完全一致的单样本抽样点）
for k in k_data:
    if k == 1:
        # 在 k=1 时，两者代表完全相同的单基因组容量抽样
        shared_noise = np.random.normal(0, 300, num_replicates)
        v_k1 = K_true + shared_noise
        
        raw_k.extend([1] * num_replicates)
        raw_V_pan.extend(v_k1)
        raw_V_core.extend(v_k1)  # 强制点完全一致
    else:
        # k > 1 时，按各自的膨胀/衰减规律生成散点
        base_pan = K_true * (k ** beta_true)
        noise_pan = np.random.normal(0, 300 + 80 * np.sqrt(k), num_replicates)
        
        base_core = Omega_true + A_true * np.exp(-kappa_true * k)
        noise_core = np.random.normal(0, 300 * np.exp(-0.05 * k), num_replicates)
        
        raw_k.extend([k] * num_replicates)
        raw_V_pan.extend(base_pan + noise_pan)
        raw_V_core.extend(base_core + noise_core)

raw_k = np.array(raw_k)
raw_V_pan = np.array(raw_V_pan)
raw_V_core = np.array(raw_V_core)

# 4. 计算各 k 梯度的均值用于拟合
pan_means = np.array([np.mean(raw_V_pan[raw_k == k]) for k in k_data])
core_means = np.array([np.mean(raw_V_core[raw_k == k]) for k in k_data])

# 5. 定义模型与拟合
def heaps_law(k, K, beta):
    return K * (k ** beta)

def exp_decay(k, Omega, A, kappa):
    return Omega + A * np.exp(-kappa * k)

popt_pan, _ = curve_fit(heaps_law, k_data, pan_means, p0=[10000, 0.3])
popt_core, _ = curve_fit(exp_decay, k_data, core_means, p0=[5000, 8000, 0.2])

K_fit, beta_fit = popt_pan
Omega_fit, A_fit, kappa_fit = popt_core

# 6. 绘图
plt.figure(figsize=(9, 5.5), dpi=120)

# A. 绘制原始抽样小叉叉
plt.scatter(raw_k, raw_V_pan, marker='x', color='#2b5c8f', s=10, 
            linewidths=0.8, alpha=0.7, label='Pangenome iterations (x)')

plt.scatter(raw_k, raw_V_core, marker='x', color="#ff6f00", s=10, 
            linewidths=0.8, alpha=0.7, label='Core genome iterations (x)')

# B. 绘制拟合平滑曲线
k_smooth = np.linspace(1, 34, 200)
V_pan_fit = heaps_law(k_smooth, K_fit, beta_fit)
V_core_fit = exp_decay(k_smooth, Omega_fit, A_fit, kappa_fit)

plt.plot(k_smooth, V_pan_fit, color="#000000", linewidth=1, 
         label=f'Pangenome fit: $V(n) = {K_fit:.0f} \\cdot k^{{{beta_fit:.3f}}}$')

plt.plot(k_smooth, V_core_fit, color="#000000", linewidth=1, 
         label=f'Core fit: $Y(n) = {Omega_fit:.0f} + {A_fit:.0f} \\cdot e^{{-{kappa_fit:.3f} k}}$')

# 图形美化
plt.xlabel("Number of Genomes (n)", fontsize=11)
plt.ylabel("Number of Orthogroups", fontsize=11)
plt.xticks(np.arange(2, 35, 2))
plt.grid(True, linestyle='--', alpha=0.4)
plt.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5, loc='center right')

plt.tight_layout()
plt.show()