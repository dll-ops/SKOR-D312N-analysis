import re
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt


def pick_file(title: str) -> Path:
    """弹窗选文件（优先 tkinter；不行就用 macOS osascript）"""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        p = filedialog.askopenfilename(
            title=title,
            filetypes=[("Text/RTF files", "*.txt *.rtf"), ("All files", "*.*")]
        )
        if not p:
            raise SystemExit("取消选择，退出。")
        return Path(p)
    except Exception:
        import subprocess
        script = f'POSIX path of (choose file with prompt "{title}")'
        out = subprocess.check_output(["osascript", "-e", script], text=True).strip()
        if not out:
            raise SystemExit("取消选择，退出。")
        return Path(out)


# 三字母 -> 单字母
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 分类
GROUPS = {
    "Positive": set("KRH"),
    "Negative": set("DE"),
    "Polar": set("NQSTY"),
    "Hydrophobic": set("AVLIFMWPCG"),
}

# 从任意文本（含 RTF）里抓 “链ID  三字母残基  编号”
# 例如：A ASP 312
PAT = re.compile(r"\b([A-Za-z0-9])\s+([A-Za-z]{3})\s+(\d+)\b")


def read_unique_residues(file_path: Path) -> dict:
    """读取文件，返回去重后的残基字典：{(chain, resnum): aa1}"""
    text = file_path.read_text(errors="ignore")
    uniq = {}
    for chain, aa3, num in PAT.findall(text):
        aa3 = aa3.upper()
        aa1 = AA3_TO_1.get(aa3)
        if aa1:
            uniq[(chain, int(num))] = aa1
    return uniq


def count_groups(uniq_dict: dict) -> Counter:
    """把残基字典按分类计数"""
    c = Counter({k: 0 for k in GROUPS})
    for aa1 in uniq_dict.values():
        for g, s in GROUPS.items():
            if aa1 in s:
                c[g] += 1
                break
    return c


def identify_state(uniq_dict: dict) -> str:
    """根据 312 位点是 D 还是 N 自动标注 WT/Mut"""
    aa312 = None
    for (chain, num), aa1 in uniq_dict.items():
        if num == 312:
            aa312 = aa1
            break
    if aa312 == "D":
        return "WT(ASP312)"
    if aa312 == "N":
        return "Mut(ASN312)"
    return f"Unknown312({aa312})"


def main():
    f1 = pick_file("选择第一个 env 文件（txt/rtf 都行）")
    f2 = pick_file("选择第二个 env 文件（txt/rtf 都行）")

    d1 = read_unique_residues(f1)
    d2 = read_unique_residues(f2)

    label1, c1 = identify_state(d1), count_groups(d1)
    label2, c2 = identify_state(d2), count_groups(d2)

    labels = list(GROUPS.keys())
    vals1 = [c1[l] for l in labels]
    vals2 = [c2[l] for l in labels]

    # 画图（只画一次、只保存一次）
    x = list(range(len(labels)))
    w = 0.35

    plt.figure()
    plt.bar(x, vals1, w, label=label1)
    plt.bar([i + w for i in x], vals2, w, label=label2)
    plt.xticks([i + w / 2 for i in x], labels)
    plt.ylabel("Residue count (within 8 Å of 312)")
    plt.legend()
    plt.tight_layout()

    out = f1.parent / "plot_env_composition_312.png"
    plt.savefig(out, dpi=300)  # ✅ 只保存这一次
    plt.close()

    print("Done ✅ 输出：", out)
    print(label1, c1)
    print(label2, c2)


if __name__ == "__main__":
    main()