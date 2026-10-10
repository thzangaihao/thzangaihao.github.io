import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 设置随机种子以确保可重复性
np.random.seed(42)

# 1. 模拟样本梯度 k (1 到 34 个基因组)
k_data = np.arange(1, 35)

# 2. 模拟 C(n, k) 重抽样的原始数据点 (每个 k 梯度模拟 15 次重抽样)
num_replicates = 15
raw_k = []
raw_V = []

for k in k_data:
    # 模拟该 k 下的基准均值 (Heap's Law 真实轨迹)
    base_mean = 14000 * (k ** 0.26)
    # 添加抽样方差 (k 较小时方差稍大)
    noise = np.random.normal(0, 400 + 100 * np.sqrt(k), num_replicates)
    V_samples = base_mean + noise
    
    raw_k.extend([k] * num_replicates)
    raw_V.extend(V_samples)

raw_k = np.array(raw_k)
raw_V = np.array(raw_V)

# 3. 计算每个 k 梯度的均值用于拟合
V_means = np.array([np.mean(raw_V[raw_k == k]) for k in k_data])

# 4. 定义希普夫定律模型并进行曲线拟合
def heaps_law(k, K, beta):
    return K * (k ** beta)

popt, pcov = curve_fit(heaps_law, k_data, V_means, p0=[10000, 0.3])
K_fit, beta_fit = popt

# 5. 绘图
plt.figure(figsize=(8, 5), dpi=120)

# A. 绘制原始抽样散点 (用小叉叉 'x' 表示)
plt.scatter(raw_k, raw_V, marker='x', color='#3B6EA8', s=10, alpha=0.7, linewidths=0.8,
            label='Subsampling iterations (x)')

# B. 绘制希普夫拟合黑线
k_smooth = np.linspace(1, 34, 200)
V_fit = heaps_law(k_smooth, K_fit, beta_fit)
plt.plot(k_smooth, V_fit, color='black', linewidth=1.5, 
         label=f'Fitted Heap\'s Curve: $Y = {K_fit:.0f} \\cdot N^{{{beta_fit:.3f}}}$')

# 图形美化
plt.title("Pangenome Accumulation Curve (Rarefaction)", fontsize=12, fontweight='bold')
plt.xlabel("Genome number (N)", fontsize=11)
plt.ylabel("Orthogroup number (Y)", fontsize=11)
plt.xticks(np.arange(2, 35, 2))
plt.grid(True, linestyle='--', alpha=0.4)
plt.legend(frameon=True, facecolor='white', framealpha=0.9)

plt.tight_layout()
plt.show()