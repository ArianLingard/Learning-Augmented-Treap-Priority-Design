# -*- coding: utf-8 -*-
"""
plot_all.py —— 读取 C++ 程序 (LA_BST_CSV.cpp) 输出的 7 个 CSV，为每个实验生成一张 PDF
（fig1~fig7），并在控制台打印表 1-表 6。

用法：在项目根目录执行  python notebooks/plot_all.py
输出：figures/fig1_e1_basic.pdf ... figures/fig7_e7_bucketing.pdf
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------
# 设计系统参数（dataviz 参考调色板，浅色模式）
# ---------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

CAT = {
    "Random": "#2a78d6",   # blue
    "Chen":   "#eb3434",   # red
    "f1":     "#d2c400",   # yellow
    "f2":     "#00b812bd",   # green
    "f3":     "#eb38e5",   # pink
    "f4":     "#9007C6",   # purple
}
METHOD_ORDER = ["Random", "Chen", "f1", "f2", "f3", "f4"]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "grid.color": GRID,
    "grid.linestyle": "--",
    "grid.linewidth": 0.6,
    "axes.grid": True,
    "axes.axisbelow": True,
    "font.size": 11,
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def read(name):
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        print(f"[warn] missing {path} — 先运行 C++ 程序")
        return None
    return pd.read_csv(path)


def method_color(m):
    return CAT.get(m, "#888888")


def save(fig, fname):
    fig.savefig(os.path.join(FIG, fname), format="pdf")
    plt.close(fig)
    print(f"  wrote figures/{fname}")


# ---------------------------------------------------------------
# 图 1：E1 常规分布基础性能（Zipf α=1 分组柱状图）
# ---------------------------------------------------------------
def fig1(e1):
    if e1 is None:
        return
    df = e1[(e1["distribution"] == "Zipf") & (e1["param"] == 1.0)]
    methods = [m for m in METHOD_ORDER if m in set(df["method"])]
    ns = sorted(df["n"].unique())
    n_colors = {ns[0]: "#86b6ef", ns[1]: "#1c5cab"}

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(methods))
    width = 0.38
    for i, n in enumerate(ns):
        sub = df[df["n"] == n].set_index("method")
        vals = [sub.loc[m, "mean_comparisons"] for m in methods]
        errs = [sub.loc[m, "std_comparisons"] for m in methods]
        bx = x + (i - 0.5) * width
        ax.bar(bx, vals, width, color=n_colors[n], label=f"n = {n}", zorder=3)
        ax.errorbar(bx, vals, yerr=errs, fmt="none",
                    ecolor=INK2, elinewidth=0.8, capsize=3, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_xlabel("方法")
    ax.set_ylabel("总比较次数（m = 10⁵）")
    ax.set_title("图 1  E1  Zipf α = 1 下各方法总比较次数")
    ax.legend(title="规模")
    ax.margins(x=0.05)
    fig.tight_layout()
    save(fig, "fig1_e1_basic.pdf")


# ---------------------------------------------------------------
# 图 2：E2 绝对量纲漂移鲁棒性（比较次数 + Kendall τ）
# ---------------------------------------------------------------
def fig2(e2):
    if e2 is None:
        return
    df = e2.copy()
    df["label"] = df.apply(lambda r: f"{r['transform']}={r['param']}", axis=1)
    labels = sorted(df["label"].unique())
    methods = [m for m in METHOD_ORDER if m in set(df["method"])]
    x = np.arange(len(labels))
    w = 0.18

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # 左：比较次数
    for i, m in enumerate(methods):
        sub = df[df["method"] == m].set_index("label")
        vals = [sub.loc[l, "mean_comparisons"] for l in labels]
        errs = [sub.loc[l, "std_comparisons"] for l in labels]
        bx = x + (i - (len(methods) - 1) / 2) * w
        axes[0].bar(bx, vals, w, color=method_color(m), label=m, zorder=3)
        axes[0].errorbar(bx, vals, yerr=errs, fmt="none",
                         ecolor=INK2, elinewidth=0.7, capsize=2, zorder=4)
    axes[0].set_xticks(x); axes[0].set_xticklabels(labels)
    axes[0].set_xlabel("量纲变换"); axes[0].set_ylabel("总比较次数")
    axes[0].set_title("比较次数")
    # 右：Kendall τ
    for i, m in enumerate(methods):
        sub = df[df["method"] == m].set_index("label")
        vals = [sub.loc[l, "mean_kendall_tau"] for l in labels]
        bx = x + (i - (len(methods) - 1) / 2) * w
        axes[1].bar(bx, vals, w, color=method_color(m), label=m, zorder=3)
    axes[1].axhline(1.0, color=AXIS, lw=0.8, ls=":")
    axes[1].set_xticks(x); axes[1].set_xticklabels(labels)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel("量纲变换"); axes[1].set_ylabel("Kendall τ")
    axes[1].set_title("拓扑不变性（τ = 1 表示拓扑不变）")

    axes[0].legend(ncol=2)
    fig.suptitle("图 2  E2  绝对量纲漂移的鲁棒性", fontsize=13)
    fig.tight_layout()
    save(fig, "fig2_e2_drift.pdf")


# ---------------------------------------------------------------
# 图 3：E3 右端奇点与数值溢出压力测试
# ---------------------------------------------------------------
def fig3(e3):
    if e3 is None:
        return
    df = e3.copy()
    df["max_priority"] = df["max_priority"].replace([np.inf, -np.inf], np.nan)
    ws = sorted(df["w"].unique())
    wlabels = [f"{w:g}" for w in ws]
    methods = [m for m in METHOD_ORDER if m in set(df["method"])]
    x = np.arange(len(ws))

    fig, ax = plt.subplots(figsize=(8, 5))
    for m in methods:
        sub = df[df["method"] == m].set_index("w")
        vals = [sub.loc[w, "max_priority"] for w in ws]
        ax.plot(x, vals, marker="o", markersize=5, lw=1.6, color=method_color(m), label=m, zorder=3)
    # 标注 Chen 在 w→1 的溢出
    ax.annotate("Chen → ∞（溢出）", xy=(len(ws) - 1, 32), xytext=(len(ws) - 2.2, 30),
                color=CAT["Chen"], arrowprops=dict(arrowstyle="->", color=CAT["Chen"]))
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(wlabels)
    ax.set_xlabel("预测权重 w（右端奇点）")
    ax.set_ylabel("最大 priority")
    ax.set_ylim(-6, 34)
    ax.set_title("图 3  E3  右端奇点压力测试：Chen 发散，其余有界")
    ax.legend()
    fig.tight_layout()
    save(fig, "fig3_e3_singularity.pdf")


# ---------------------------------------------------------------
# 图 4：E4 理论静态最优代价分析（ρ）
# ---------------------------------------------------------------
def fig4(e4):
    if e4 is None:
        return
    df = e4[e4["n"] == 5000].copy()
    df["dist"] = df.apply(
        lambda r: f"{r['distribution']}" + (f" α={r['param']}" if r["distribution"] == "Zipf" else ""),
        axis=1)
    dists = sorted(df["dist"].unique())
    methods = [m for m in METHOD_ORDER if m in set(df["method"])]
    x = np.arange(len(dists))
    w = 0.14

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, m in enumerate(methods):
        sub = df[df["method"] == m].set_index("dist")
        vals = [sub.loc[d, "mean_rho"] for d in dists]
        errs = [sub.loc[d, "std_rho"] for d in dists]
        bx = x + (i - (len(methods) - 1) / 2) * w
        ax.bar(bx, vals, w, color=method_color(m), label=m, zorder=3)
        ax.errorbar(bx, vals, yerr=errs, fmt="none",
                    ecolor=INK2, elinewidth=0.7, capsize=2, zorder=4)
    ax.axhline(1.0, color=AXIS, lw=0.8, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(dists)
    ax.set_xlabel("分布")
    ax.set_ylabel("静态最优比 ρ = 比较次数 / 香农熵下界")
    ax.set_title("图 4  E4  静态最优比（n = 5000，ρ 越接近常数越优）")
    ax.legend(ncol=3)
    fig.tight_layout()
    save(fig, "fig4_e4_static_optimal.pdf")


# ---------------------------------------------------------------
# 图 5：E5 非完美预测器（KL 噪声）鲁棒性
# ---------------------------------------------------------------
def fig5(e5):
    if e5 is None:
        return
    df = e5.copy()
    methods = [m for m in METHOD_ORDER if m in set(df["method"])]
    fig, ax = plt.subplots(figsize=(8, 5))
    for m in methods:
        sub = df[df["method"] == m].sort_values("kl")
        ax.errorbar(sub["kl"], sub["mean_comparisons"], yerr=sub["std_comparisons"],
                    color=method_color(m), marker="o", markersize=5, lw=1.6,
                    capsize=3, label=m, zorder=3,
                    elinewidth=0.8, ecolor=method_color(m))
    ax.set_xlabel("KL 散度  D_KL(P ‖ Q)")
    ax.set_ylabel("总比较次数（m = 10⁵）")
    ax.set_title("图 5  E5  预测器 KL 噪声下的鲁棒性")
    ax.set_xticks(sorted(df["kl"].unique()))
    ax.legend(ncol=2)
    fig.tight_layout()
    save(fig, "fig5_e5_kl_noise.pdf")


# ---------------------------------------------------------------
# 图 6：E6 自适应扰动 vs 固定幅度扫描（跨桶逆序率 + 查询代价）
# ---------------------------------------------------------------
def fig6(e6):
    if e6 is None:
        return
    df = e6.copy()
    methods = ["f3", "f4"]
    fixed_order = ["0.01", "0.25", "0.49", "0.75", "0.99", "1.25", "1.75", "2.0"]
    x = np.arange(len(fixed_order))
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8))

    for m in methods:
        sub = df[df["method"] == m].set_index("perturb")
        fixed = sub.loc[fixed_order]
        # 左：跨桶逆序率（机制指标）
        axes[0].plot(x, fixed["mean_inversion_ratio"], marker="o", markersize=5, lw=1.8,
                     color=CAT[m], label=m, zorder=3)
        axes[0].fill_between(x,
            fixed["mean_inversion_ratio"] - fixed["std_inversion_ratio"],
            fixed["mean_inversion_ratio"] + fixed["std_inversion_ratio"],
            color=CAT[m], alpha=0.15, zorder=2)
        # 右：查询代价（后果指标）
        axes[1].errorbar(x, fixed["mean_comparisons"], yerr=fixed["std_comparisons"],
                         color=CAT[m], marker="o", markersize=5, lw=1.8,
                         capsize=3, elinewidth=0.8, ecolor=CAT[m], label=m, zorder=3)

    # 自适应参考（虚线）：恒在安全区
    for m in methods:
        sub = df[(df["method"] == m) & (df["perturb"] == "adaptive")]
        if len(sub) == 0:
            continue
        axes[0].axhline(sub["mean_inversion_ratio"].iloc[0], color=CAT[m], ls="--", lw=1.1, alpha=0.7)
        axes[1].axhline(sub["mean_comparisons"].iloc[0], color=CAT[m], ls="--", lw=1.1, alpha=0.7)

    for ax in axes:
        ax.axvline(4.5, color=AXIS, lw=0.9, ls=":")   # φ = 1（桶宽）分界
        ax.set_xticks(x)
        ax.set_xticklabels(fixed_order)
        ax.set_xlabel("固定扰动幅度 φ（桶宽 = 1）")

    axes[0].set_ylabel("跨桶逆序率")
    axes[0].set_ylim(-0.005, 0.14)
    axes[0].set_title("机制指标：跨桶逆序率（φ ≥ 1 从 0 跳升）")
    axes[1].set_ylabel("总比较次数（m = 10⁵）")
    axes[1].set_title("后果指标：查询代价")

    dash = plt.Line2D([], [], color=INK2, ls="--", lw=1.2, label="自适应（虚线，恒安全）")
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles + [dash], labels=labels + [dash.get_label()], ncol=3, fontsize=9)

    fig.suptitle("图 6  E6  自适应扰动：固定幅度跨过桶宽（φ ≥ 1）即破坏学习排序", fontsize=13)
    fig.tight_layout()
    save(fig, "fig6_e6_adaptive_perturb.pdf")


# ---------------------------------------------------------------
# 图 7：E7 分桶 + 扰动必要性消融（对抗键序）
# ---------------------------------------------------------------
def fig7(e7):
    if e7 is None:
        return
    df = e7.copy()
    labels = {"no_bucket_no_perturb": "无分桶\n无扰动",
              "bucket_no_perturb": "分桶\n无扰动",
              "bucket_perturb": "分桶\n+固定扰动",
              "bucket_adaptive": "分桶\n+自适应扰动"}
    settings = list(df["setting"])
    x = np.arange(len(settings))
    colors = ["#d03b3b", "#d03b3b", "#0ca30c", "#0ca30c"]   # 退化=红，平衡=绿
    metrics = [
        ("mean_height", "std_height", "树高（n = 5000）"),
        ("mean_comparisons", "std_comparisons", "总比较次数（m = 2×10⁴）"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for ax, (mc, _, ylab) in zip(axes, metrics):
        vals = df[mc].values
        ax.bar(x, vals, 0.5, bottom=1, color=colors, zorder=3)
        ax.set_yscale("log")
        ax.set_ylim(bottom=1)
        for xi, v in zip(x, vals):
            ax.annotate(f"{v:,.0f}", (xi, v), textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=8.5, color=INK2)
        ax.set_xticks(x)
        ax.set_xticklabels([labels[s] for s in settings])
        ax.set_ylabel(ylab)
    axes[0].set_title("树高（退化 = 链，≈ n）")
    axes[1].set_title("查询代价")
    fig.suptitle("图 7  E7  分桶 + 扰动必要性消融（对抗键序，Zipf α=1, n=5000）", fontsize=13)
    fig.tight_layout()
    save(fig, "fig7_e7_bucketing.pdf")


# ---------------------------------------------------------------
# 控制台表格（表 1 - 表 6）
# ---------------------------------------------------------------
def tables(e1, e2, e3, e4, e5, e6, e7):
    print("\n" + "=" * 70)
    print("表 1  E1 常规分布基础性能（总比较次数均值）")
    if e1 is not None:
        t = e1.pivot_table(index=["distribution", "param", "n"], columns="method",
                           values="mean_comparisons")
        print(t.round(0).to_string())
    print("\n表 2  E2 量纲漂移（Kendall τ）")
    if e2 is not None:
        t = e2.pivot_table(index=["transform", "param"], columns="method",
                           values="mean_kendall_tau")
        print(t.round(3).to_string())
    print("\n表 3  E3 右端奇点（最大 priority / 是否溢出）")
    if e3 is not None:
        print(e3.to_string(index=False))
    print("\n表 4  E4 静态最优比 ρ")
    if e4 is not None:
        t = e4.pivot_table(index=["distribution", "param", "n"], columns="method",
                           values="mean_rho")
        print(t.round(3).to_string())
    print("\n表 5  E6 自适应扰动（跨桶逆序率 / 比较次数 / 最大深度）")
    if e6 is not None:
        print(e6.to_string(index=False))
    print("\n表 6  E7 分桶消融（树高 / 比较次数）")
    if e7 is not None:
        print(e7.to_string(index=False))


def main():
    print("读取 data/*.csv ...")
    e1 = read("e1_basic.csv")
    e2 = read("e2_drift.csv")
    e3 = read("e3_singularity.csv")
    e4 = read("e4_static_optimal.csv")
    e5 = read("e5_kl_noise.csv")
    e6 = read("e6_adaptive_perturb.csv")
    e7 = read("e7_bucketing.csv")

    fig1(e1)
    fig2(e2)
    fig3(e3)
    fig4(e4)
    fig5(e5)
    fig6(e6)
    fig7(e7)
    tables(e1, e2, e3, e4, e5, e6, e7)
    print("\ndone.")


if __name__ == "__main__":
    main()
