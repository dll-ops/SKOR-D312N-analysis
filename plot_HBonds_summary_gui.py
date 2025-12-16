from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def pick_file(title: str) -> Path:
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

def find_col(cols, must_contain_all=(), must_contain_any=()):
    cols_l = {c: c.lower() for c in cols}
    for c, cl in cols_l.items():
        if all(x in cl for x in must_contain_all) and (not must_contain_any or any(x in cl for x in must_contain_any)):
            return c
    return None

csv_path = pick_file("选择 hbonds_312_summary.csv")
out_dir = csv_path.parent

df = pd.read_csv(csv_path)

# 列名自动匹配（不依赖 Å 符号）
state_col = find_col(df.columns, must_contain_any=("state",)) or df.columns[0]
count_col = find_col(df.columns, must_contain_all=("hbonds",)) or find_col(df.columns, must_contain_any=("count", "number"))
mind_col  = find_col(df.columns, must_contain_all=("min", "d"))

if count_col is None or mind_col is None:
    print("❌ 列名没匹配上。当前 CSV 列名是：")
    for c in df.columns:
        print(" -", repr(c))
    raise SystemExit(1)

# 图1：氢键数量
plt.figure()
plt.bar(df[state_col], df[count_col])
plt.ylabel("HBonds involving 312")
plt.xticks(rotation=15, ha="right")
plt.tight_layout()
plt.savefig(out_dir / "plot_hbonds_count.png", dpi=300)
plt.close()

# 图2：最短距离
plt.figure()
plt.bar(df[state_col], df[mind_col])
plt.ylabel("Min D..A (Å)")
plt.xticks(rotation=15, ha="right")
plt.tight_layout()
plt.savefig(out_dir / "plot_hbonds_min_dist.png", dpi=300)
plt.close()

print("Done ✅ 输出：")
print(out_dir / "plot_hbonds_count.png")
print(out_dir / "plot_hbonds_min_dist.png")
print("（如果你想确认脚本识别到的列名：）")
print("state_col =", state_col)
print("count_col =", count_col)
print("mind_col  =", mind_col)