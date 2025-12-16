import re
import csv
from pathlib import Path
from collections import defaultdict

# 兼容你这种 ChimeraX H-bonds 输出行
pat = re.compile(
    r"#(?P<m1>\d+)/(?P<c1>\S)\s+(?P<r1>[A-Z]{3})\s+(?P<n1>\d+)\s+(?P<a1>\S+)\s+.*?"
    r"#(?P<m2>\d+)/(?P<c2>\S)\s+(?P<r2>[A-Z]{3})\s+(?P<n2>\d+)\s+(?P<a2>\S+)\s+.*?"
    r"\s(?P<da>\d+\.\d+)\s"
)

def is_312(x: str) -> bool:
    return int(x) == 312

def pick_file(title: str) -> Path:
    # 优先 tkinter 弹窗；如果你的 python 没 tkinter，就退回到 mac 的 osascript 选文件框
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        p = filedialog.askopenfilename(
            title=title,
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
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

def parse_hbonds(path: Path, label: str):
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        m = pat.search(line)
        if not m:
            continue
        d = m.groupdict()
        if not (is_312(d["n1"]) or is_312(d["n2"])):
            continue

        # 统一方向：让 312 永远在左边（center）
        if is_312(d["n2"]):
            d = {
                "m1": d["m2"], "c1": d["c2"], "r1": d["r2"], "n1": d["n2"], "a1": d["a2"],
                "m2": d["m1"], "c2": d["c1"], "r2": d["r1"], "n2": d["n1"], "a2": d["a1"],
                "da": d["da"],
            }

        center = f"{d['r1']}{d['n1']}{d['a1']}"      # e.g. ASN312ND2 / ASP312N
        partner = f"{d['r2']}{d['n2']}{d['a2']}"     # e.g. TYR308O
        da = float(d["da"])
        rows.append((label, d["c1"], center, d["c2"], partner, da))
    return rows

def main():
    mut = pick_file("选择突变体 H-bonds txt（ASN312 / model #1）")
    wt  = pick_file("选择 WT H-bonds txt（ASP312 / model #2）")

    # 输出目录：默认放到突变体 txt 所在文件夹（两份一般同目录）
    out_dir = mut.parent

    rows = []
    rows += parse_hbonds(mut, "Mut_ASN312(#1)")
    rows += parse_hbonds(wt,  "WT_ASP312(#2)")

    # 去重：同一状态下同一对(中心原子-伙伴原子)保留最短距离
    best = {}
    for state, c1, center, c2, partner, da in rows:
        key = (state, c1, center, c2, partner)
        best[key] = min(best.get(key, 1e9), da)

    detail = sorted(
        [(k[0], k[1], k[2], k[3], k[4], v) for k, v in best.items()],
        key=lambda x: (x[0], x[2], x[4], x[5])
    )

    out_detail = out_dir / "hbonds_312_detail.csv"
    with out_detail.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["State", "CenterChain", "CenterAtom(312)", "PartnerChain", "PartnerAtom", "Min_D..A(Å)"])
        for r in detail:
            w.writerow([r[0], r[1], r[2], r[3], r[4], f"{r[5]:.3f}"])

    summary = defaultdict(lambda: {"count": 0, "min_da": 999.0})
    for state, *_ , da in detail:
        summary[state]["count"] += 1
        summary[state]["min_da"] = min(summary[state]["min_da"], da)

    out_sum = out_dir / "hbonds_312_summary.csv"
    with out_sum.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["State", "HBonds_involving_312", "Min_D..A(Å)"])
        for state in summary:
            w.writerow([state, summary[state]["count"], f"{summary[state]['min_da']:.3f}"])

    print("Done ✅")
    print("输出：", out_detail)
    print("输出：", out_sum)

if __name__ == "__main__":
    main()