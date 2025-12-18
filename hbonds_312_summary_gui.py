"""
hbonds_312_summary_gui.py

用途：
1) 让用户选择两份 ChimeraX H-bonds 输出 txt（突变体 vs WT）。
2) 解析每行氢键信息，筛选出“涉及 312 位点”的记录。
3) 统一记录方向：确保 312 位点永远作为“中心（center）”放在左侧。
4) 去重：同一状态下同一对（中心原子-伙伴原子）可能出现多次，保留最短 D..A 距离。
5) 输出两份 CSV：
   - hbonds_312_detail.csv：明细（每个配对一行）
   - hbonds_312_summary.csv：汇总（每个状态一行）

输入：
- ChimeraX 的 H-bonds 输出文本（txt），要求行格式能被本脚本的正则 pat 匹配。

输出：
- 与输入文件同目录下的 hbonds_312_detail.csv 和 hbonds_312_summary.csv
"""

import re
import csv
from pathlib import Path
from collections import defaultdict

# ------------------------------------------------------------
# 正则：兼容你这种 ChimeraX H-bonds 输出行
# 说明：
# - 抓取两端参与氢键的对象：model/chain/residue_name/residue_number/atom_name
# - 并抓取 D..A 距离（da）
# 注意：
# - 这里假定 chain 是单字符（\S），residue name 是三字母大写（[A-Z]{3}）
# - 若你的输出格式变化（例如链名更长、残基名不是三字母、字段间空格变化），需要改 pat
# ------------------------------------------------------------
pat = re.compile(
    r"#(?P<m1>\d+)/(?P<c1>\S)\s+(?P<r1>[A-Z]{3})\s+(?P<n1>\d+)\s+(?P<a1>\S+)\s+.*?"
    r"#(?P<m2>\d+)/(?P<c2>\S)\s+(?P<r2>[A-Z]{3})\s+(?P<n2>\d+)\s+(?P<a2>\S+)\s+.*?"
    r"\s(?P<da>\d+\.\d+)\s"
)

def is_312(x: str) -> bool:
    """
    判断一个 residue number 是否为目标位点 312。

    参数：
        x: 从正则中抓出来的 residue number（字符串形式）
    返回：
        True 表示该 residue number == 312
    """
    return int(x) == 312

def pick_file(title: str) -> Path:
    """
    弹窗选择文件，返回 Path。

    逻辑：
    1) 优先使用 tkinter 的文件选择对话框（跨平台常用）
    2) 若 tkinter 不可用（或异常），在 macOS 上退回使用 osascript 调用系统 choose file

    参数：
        title: 文件选择对话框的提示标题
    返回：
        用户选择的文件路径（Path）

    异常/退出：
        用户取消选择时直接 SystemExit 退出（脚本终止）
    """
    # 优先 tkinter 弹窗；如果你的 python 没 tkinter，就退回到 mac 的 osascript 选文件框
    try:
        import tkinter as tk
        from tkinter import filedialog

        # 创建一个 Tk 根窗口，但隐藏它（只用文件对话框）
        root = tk.Tk()
        root.withdraw()

        # 弹出打开文件对话框
        p = filedialog.askopenfilename(
            title=title,
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        # 如果用户取消，askopenfilename 返回空字符串
        if not p:
            raise SystemExit("取消选择，已退出。")

        return Path(p)

    except Exception:
        # macOS fallback：osascript 调用 Finder 的 choose file
        import subprocess

        # 返回 POSIX 路径（/Users/...）
        script = f'POSIX path of (choose file with prompt "{title}")'
        out = subprocess.check_output(["osascript", "-e", script], text=True).strip()

        if not out:
            raise SystemExit("取消选择，已退出。")

        return Path(out)

def parse_hbonds(path: Path, label: str):
    """
    解析单个 ChimeraX H-bonds 输出 txt，返回“涉及 312 位点”的记录列表。

    参数：
        path: 输入 txt 的路径
        label: 状态标签，用于区分 WT / Mut（会写入输出 CSV 的 State 列）
    返回：
        rows: list[tuple]，每个元素为：
              (state_label, center_chain, center_atom_str, partner_chain, partner_atom_str, da_distance)

    关键逻辑：
    - 先用 pat 正则匹配每一行，取出两端信息 + 距离
    - 只保留 residue number 为 312 的记录（任一端为 312 即保留）
    - 统一方向：让 312 永远在左边（center）
    """
    rows = []

    # 读取文本，errors="ignore" 用于容错（遇到奇怪编码字符直接跳过）
    for line in path.read_text(errors="ignore").splitlines():
        m = pat.search(line)
        if not m:
            # 行格式不匹配 pat 的，直接跳过
            continue

        # 将正则命名分组转为 dict，字段包括：
        # m1,c1,r1,n1,a1,m2,c2,r2,n2,a2,da
        d = m.groupdict()

        # 只要两端都不是 312，就跳过
        if not (is_312(d["n1"]) or is_312(d["n2"])):
            continue

        # ------------------------------------------------------------
        # 统一方向：让 312 永远在左边（center）
        # 若 312 在右边（n2），就把左右两端交换
        # ------------------------------------------------------------
        if is_312(d["n2"]):
            d = {
                "m1": d["m2"], "c1": d["c2"], "r1": d["r2"], "n1": d["n2"], "a1": d["a2"],
                "m2": d["m1"], "c2": d["c1"], "r2": d["r1"], "n2": d["n1"], "a2": d["a1"],
                "da": d["da"],
            }

        # ------------------------------------------------------------
        # 构造“中心原子”和“伙伴原子”的字符串标识
        # 例子：
        # center  = ASN312ND2 或 ASP312N（由残基名+编号+原子名拼出来）
        # partner = TYR308O
        # ------------------------------------------------------------
        center = f"{d['r1']}{d['n1']}{d['a1']}"
        partner = f"{d['r2']}{d['n2']}{d['a2']}"

        # D..A 距离（Å）转 float
        da = float(d["da"])

        # 保存一条记录
        rows.append((label, d["c1"], center, d["c2"], partner, da))

    return rows

def main():
    """
    主流程：
    1) 选择突变体与 WT 的 H-bonds txt
    2) 解析两份文件，收集所有涉及 312 的氢键记录
    3) 去重并保留最短距离
    4) 输出明细 CSV + 汇总 CSV
    """
    # 选择输入文件：提示语里写明了你默认的模型编号语义（#1/#2）
    mut = pick_file("选择突变体 H-bonds txt（ASN312 / model #1）")
    wt  = pick_file("选择 WT H-bonds txt（ASP312 / model #2）")

    # 输出目录：默认放到突变体 txt 所在文件夹（通常两份 txt 会在同一目录）
    out_dir = mut.parent

    # rows：收集原始解析结果（含重复）
    rows = []
    rows += parse_hbonds(mut, "Mut_ASN312(#1)")
    rows += parse_hbonds(wt,  "WT_ASP312(#2)")

    # ------------------------------------------------------------
    # 去重策略：
    # key = (state, center_chain, center_atom, partner_chain, partner_atom)
    # value = 该 key 出现过的最小 da 距离
    #
    # 这样做的意义：
    # - 同一对原子（center-partner）可能在文本里重复出现（不同几何、不同记录行等）
    # - 我们只保留最“强/近”的那一次（最短距离）
    # ------------------------------------------------------------
    best = {}
    for state, c1, center, c2, partner, da in rows:
        key = (state, c1, center, c2, partner)
        best[key] = min(best.get(key, 1e9), da)

    # 将 best 字典转回列表，并排序，便于输出稳定（同输入可复现一致顺序）
    detail = sorted(
        [(k[0], k[1], k[2], k[3], k[4], v) for k, v in best.items()],
        key=lambda x: (x[0], x[2], x[4], x[5])  # 按 state、center、partner、距离排序
    )

    # ------------------------------------------------------------
    # 输出明细表：hbonds_312_detail.csv
    # ------------------------------------------------------------
    out_detail = out_dir / "hbonds_312_detail.csv"
    with out_detail.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["State", "CenterChain", "CenterAtom(312)", "PartnerChain", "PartnerAtom", "Min_D..A(Å)"])
        for r in detail:
            # r[5] 为 da 距离，输出格式保留 3 位小数
            w.writerow([r[0], r[1], r[2], r[3], r[4], f"{r[5]:.3f}"])

    # ------------------------------------------------------------
    # 汇总统计：
    # - count：该状态下涉及 312 的“去重后配对数”
    # - min_da：该状态下所有配对里的最小 D..A 距离
    # ------------------------------------------------------------
    summary = defaultdict(lambda: {"count": 0, "min_da": 999.0})
    for state, *_ , da in detail:
        summary[state]["count"] += 1
        summary[state]["min_da"] = min(summary[state]["min_da"], da)

    # 输出汇总表：hbonds_312_summary.csv
    out_sum = out_dir / "hbonds_312_summary.csv"
    with out_sum.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["State", "HBonds_involving_312", "Min_D..A(Å)"])
        for state in summary:
            w.writerow([state, summary[state]["count"], f"{summary[state]['min_da']:.3f}"])

    # 控制台提示
    print("Done ✅")
    print("输出：", out_detail)
    print("输出：", out_sum)

# 脚本入口：直接运行此文件时执行 main()
if __name__ == "__main__":
    main()