"""
plot_env_composition_gui.py

用途：
对 SKOR 蛋白中 312 位点邻域的“残基环境组成”进行统计与可视化，
用于比较 WT 与突变体（如 D312N）在局部理化环境上的差异。

本脚本只负责：
- 读取已整理好的环境统计 CSV
- 按残基理化性质分类
- 计算比例 / 数量
- 生成对比图（WT vs Mut）

不负责：
- 结构解析
- 距离计算
- 残基邻域判定
"""

import csv
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 残基分类规则（理化性质）
# 这是整个分析的“生物学假设核心”
# ------------------------------------------------------------
RESIDUE_CATEGORIES = {
    "hydrophobic": {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"},
    "polar": {"SER", "THR", "ASN", "GLN", "TYR", "CYS"},
    "positive": {"LYS", "ARG", "HIS"},
    "negative": {"ASP", "GLU"},
}

def classify_residue(res_name: str) -> str:
    """
    根据三字母残基名，将其归类为某一理化性质类别。

    参数：
        res_name: 三字母残基名（如 'ASP', 'ASN'）
    返回：
        类别字符串：
        'hydrophobic' / 'polar' / 'positive' / 'negative'
        若不在定义中，返回 'other'
    """
    for category, members in RESIDUE_CATEGORIES.items():
        if res_name in members:
            return category
    return "other"

def pick_file(title: str) -> Path:
    """
    GUI 弹窗选择 CSV 文件。

    参数：
        title: 对话框标题
    返回：
        Path 对象
    """
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()

        p = filedialog.askopenfilename(
            title=title,
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not p:
            raise SystemExit("取消选择，已退出。")
        return Path(p)

    except Exception:
        import subprocess
        script = f'POSIX path of (choose file with prompt "{title}")'
        out = subprocess.check_output(["osascript", "-e", script], text=True).strip()
        if not out:
            raise SystemExit("取消选择，已退出。")
        return Path(out)

def load_env_table(path: Path) -> pd.DataFrame:
    """
    读取环境统计 CSV。

    CSV 至少应包含以下字段：
    - State：WT / Mut
    - ResidueName：三字母残基名
    - Is312Neighbor：是否属于 312 位点邻域（True/False 或 1/0）

    参数：
        path: CSV 文件路径
    返回：
        pandas DataFrame
    """
    df = pd.read_csv(path)

    required_cols = {"State", "ResidueName", "Is312Neighbor"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV 缺少必要字段：{required_cols}")

    return df

def summarize_composition(df: pd.DataFrame) -> pd.DataFrame:
    """
    对 312 位点邻域残基进行分类统计。

    逻辑：
    - 只保留 Is312Neighbor == True 的行
    - 按 State（WT / Mut）分组
    - 统计每种理化类别的数量
    - 同时计算比例

    返回：
        每一行一个 State，每一列一个类别（数量或比例）
    """
    df = df[df["Is312Neighbor"] == True].copy()

    # 新增一列：ResidueCategory
    df["Category"] = df["ResidueName"].apply(classify_residue)

    # 按 State + Category 统计数量
    count_table = (
        df.groupby(["State", "Category"])
          .size()
          .unstack(fill_value=0)
    )

    # 转成比例
    ratio_table = count_table.div(count_table.sum(axis=1), axis=0)

    # 给列名加后缀，避免歧义
    ratio_table = ratio_table.add_suffix("_ratio")
    count_table = count_table.add_suffix("_count")

    # 合并数量与比例
    summary = pd.concat([count_table, ratio_table], axis=1)

    return summary

def plot_composition(summary: pd.DataFrame, out_dir: Path):
    """
    绘制 WT vs Mut 的环境组成对比图（堆叠柱状图）。

    参数：
        summary: summarize_composition() 输出的 DataFrame
        out_dir: 图像保存目录
    """
    categories = [c for c in summary.columns if c.endswith("_ratio")]

    ax = summary[categories].plot(
        kind="bar",
        stacked=True,
        figsize=(6, 4)
    )

    ax.set_ylabel("Residue composition ratio")
    ax.set_xlabel("State")
    ax.set_title("Environment composition around residue 312")

    plt.tight_layout()

    out_path = out_dir / "env_312_composition_stacked.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

def main():
    """
    主流程：
    1) 选择环境统计 CSV
    2) 汇总 312 位点邻域残基组成
    3) 输出汇总表 + 可视化图像
    """
    csv_path = pick_file("选择 312 位点环境统计 CSV")
    out_dir = csv_path.parent

    df = load_env_table(csv_path)
    summary = summarize_composition(df)

    # 输出汇总表
    out_csv = out_dir / "env_312_composition_summary.csv"
    summary.to_csv(out_csv)

    # 画图
    plot_composition(summary, out_dir)

    print("Done ✅")
    print("输出：", out_csv)

if __name__ == "__main__":
    main()