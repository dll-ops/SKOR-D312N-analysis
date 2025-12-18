"""
plot_env_composition_gui.py

功能说明：
本脚本用于对比两个“env 文件”（txt/rtf，任意文本也行）中出现的残基集合，
并把这些残基按理化性质分类计数，最后画一张柱状对比图并保存为 PNG。

核心思路：
1) 让用户通过 GUI 选择两个输入文件（文件通常是某个工具/脚本导出的“环境残基列表”）。
2) 在文件全文中用正则 PAT 扫描形如 “链ID  三字母残基  编号” 的片段，例如：
   A ASP 312
3) 将匹配到的三字母残基名转为单字母（AA3_TO_1），并按 (chain, resnum) 去重：
   - uniq[(chain, resnum)] = aa1
   - 同一条链同一编号只保留一次（最后一次覆盖前一次，等价于去重）
4) 将去重后的残基集合按 GROUPS 分类计数（Positive/Negative/Polar/Hydrophobic）。
5) 根据 312 位点残基是 D 还是 N 自动标注状态：
   - D => WT(ASP312)
   - N => Mut(ASN312)
   - 其他/缺失 => Unknown312(...)
6) 绘制两组柱状图并保存为 plot_env_composition_312.png（保存到第一个文件所在目录）。

重要限制：
- 统计对象 = 输入文件里能被 PAT 抓到的所有残基条目（不是“严格的 8Å 邻域”计算结果）
  也就是说：邻域的真实性完全取决于你输入文件本身是否真的是“8Å 内残基列表”。
- PAT 只支持“单字符链ID”（[A-Za-z0-9]）+ “三字母残基名” + “编号”这种布局。
- GROUPS 分类完全由你定义的单字母集合决定，且残基只会落入第一个匹配到的组。
"""

import re                     # 正则表达式：用于从文本中提取“链 残基 编号”
from pathlib import Path       # 路径对象：跨平台文件路径拼接/读写
from collections import Counter # 计数器：按分类统计残基数量
import matplotlib.pyplot as plt # 画图：输出柱状图 PNG


def pick_file(title: str) -> Path:
    """
    让用户弹窗选择文件，并返回文件路径（Path）。

    实现策略：
    - 优先使用 tkinter 的文件对话框（跨平台）
    - 若 tkinter 不可用或报错，则在 macOS 上用 osascript 调用系统 choose file

    参数：
        title: 对话框标题文本
    返回：
        用户选择的文件路径（Path）
    退出：
        用户取消选择会直接 SystemExit 终止脚本（避免后续空路径导致的隐性错误）
    """
    try:
        import tkinter as tk                 # GUI 库（部分 Python 环境可能不带）
        from tkinter import filedialog       # 文件选择对话框
        root = tk.Tk()                       # Tk 根窗口
        root.withdraw()                      # 隐藏根窗口（只弹文件选择框）
        p = filedialog.askopenfilename(      # 弹窗选择文件
            title=title,
            filetypes=[("Text/RTF files", "*.txt *.rtf"), ("All files", "*.*")]
        )
        if not p:
            raise SystemExit("取消选择，退出。")
        return Path(p)                       # 转 Path，便于后面 .parent / 拼路径
    except Exception:
        import subprocess                    # 调用系统命令（用于 macOS osascript）
        script = f'POSIX path of (choose file with prompt "{title}")'  # AppleScript
        out = subprocess.check_output(["osascript", "-e", script], text=True).strip()
        if not out:
            raise SystemExit("取消选择，退出。")
        return Path(out)


# 三字母 -> 单字母（用于分类）
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 残基理化性质分类（以单字母集合表示）
# 注意：这里的分类是“你定义的规则”，不是任何“自动推断”。
GROUPS = {
    "Positive": set("KRH"),          # 正电：Lys/Arg/His
    "Negative": set("DE"),           # 负电：Asp/Glu
    "Polar": set("NQSTY"),           # 极性：Asn/Gln/Ser/Thr/Tyr（你没把 C 算进来）
    "Hydrophobic": set("AVLIFMWPCG") # 疏水/非极性：含 P/C/G（你的定义就是这样）
}

# 从任意文本（含 RTF）里抓 “链ID  三字母残基  编号”
# 例如：A ASP 312
# PAT 的三个捕获组分别对应：
#   1) 链ID：单字符（字母/数字）
#   2) 三字母残基名：长度为 3 的字母串
#   3) 残基编号：整数
PAT = re.compile(r"\b([A-Za-z0-9])\s+([A-Za-z]{3})\s+(\d+)\b")


def read_unique_residues(file_path: Path) -> dict:
    """
    扫描输入文件全文，提取并去重残基，返回字典：

        uniq[(chain, resnum)] = aa1

    去重语义：
    - 同一 (chain, resnum) 在文本里出现多次，只保留最后一次写入的 aa1
    - 因为 aa1 对同一编号理论上应一致，所以这等价于“去重”

    参数：
        file_path: 输入文件路径（txt/rtf/任意文本）
    返回：
        uniq: dict，key 为 (chain, int(resnum))，value 为单字母 aa1
    """
    text = file_path.read_text(errors="ignore")  # 读取全文；遇到乱码直接跳过（不抛异常）
    uniq = {}                                    # 去重字典： (chain,num) -> aa1
    for chain, aa3, num in PAT.findall(text):    # PAT.findall 返回所有匹配到的三元组
        aa3 = aa3.upper()                        # 三字母残基名统一大写，匹配 AA3_TO_1
        aa1 = AA3_TO_1.get(aa3)                  # 查表转单字母；未知残基会得到 None
        if aa1:
            uniq[(chain, int(num))] = aa1        # 保存去重结果（按链+编号唯一）
    return uniq


def count_groups(uniq_dict: dict) -> Counter:
    """
    将去重后的残基集合按 GROUPS 分类计数。

    参数：
        uniq_dict: read_unique_residues() 的返回值
    返回：
        Counter，其中键为分类名（Positive/Negative/Polar/Hydrophobic），值为计数
        初始化时保证每个分类键都存在（即使为 0）
    """
    c = Counter({k: 0 for k in GROUPS})           # 预置所有分类为 0，保证输出稳定
    for aa1 in uniq_dict.values():                # 遍历所有单字母残基
        for g, s in GROUPS.items():               # 逐类匹配
            if aa1 in s:
                c[g] += 1                         # 计入第一个匹配的分类
                break                             # 一旦匹配就停止，避免重复计数
    return c


def identify_state(uniq_dict: dict) -> str:
    """
    根据 312 位点的残基类型自动给样本打标签。

    逻辑：
    - 在 uniq_dict 中查找 num == 312 的残基（不关心链，只取第一个遇到的）
    - 如果 aa1 为 D -> WT(ASP312)
    - 如果 aa1 为 N -> Mut(ASN312)
    - 其他情况 -> Unknown312(...)

    参数：
        uniq_dict: read_unique_residues() 的返回值
    返回：
        状态标签字符串
    """
    aa312 = None                                  # 用于保存 312 位点的单字母残基（找不到则保持 None）
    for (chain, num), aa1 in uniq_dict.items():   # 遍历去重后的 (chain,num)->aa1
        if num == 312:
            aa312 = aa1
            break                                 # 找到一个就停（注意：不同链同时出现 312 时只取第一个）
    if aa312 == "D":
        return "WT(ASP312)"
    if aa312 == "N":
        return "Mut(ASN312)"
    return f"Unknown312({aa312})"                 # aa312 为 None 时会输出 Unknown312(None)


def main():
    """
    主流程：
    1) 选择两份 env 文件（文本/rtf）
    2) 提取去重残基集合
    3) 自动识别样本状态（WT / Mut / Unknown）
    4) 分类计数
    5) 绘图并保存 PNG
    """
    f1 = pick_file("选择第一个 env 文件（txt/rtf 都行）")  # 第一个输入文件路径
    f2 = pick_file("选择第二个 env 文件（txt/rtf 都行）")  # 第二个输入文件路径

    d1 = read_unique_residues(f1)                # 第一个文件的去重残基字典
    d2 = read_unique_residues(f2)                # 第二个文件的去重残基字典

    label1, c1 = identify_state(d1), count_groups(d1)  # 样本1：状态标签 + 各分类计数
    label2, c2 = identify_state(d2), count_groups(d2)  # 样本2：状态标签 + 各分类计数

    labels = list(GROUPS.keys())                 # 分类顺序（决定画图的 x 轴顺序）
    vals1 = [c1[l] for l in labels]              # 样本1每一类的计数（与 labels 对齐）
    vals2 = [c2[l] for l in labels]              # 样本2每一类的计数（与 labels 对齐）

    # 画图：并排柱状图（两组数据对比）
    x = list(range(len(labels)))                 # x 轴位置：0..(n-1)
    w = 0.35                                     # 单个柱子的宽度

    plt.figure()                                 # 新建画布（默认尺寸）
    plt.bar(x, vals1, w, label=label1)           # 第一组柱：位置 x
    plt.bar([i + w for i in x], vals2, w, label=label2)  # 第二组柱：位置右移 w
    plt.xticks([i + w / 2 for i in x], labels)   # x tick 放在两组柱的中间
    plt.ylabel("Residue count (within 8 Å of 312)")  # y 轴标签（是否真 8Å 取决于输入文件）
    plt.legend()                                 # 图例：显示 WT/Mut 标签
    plt.tight_layout()                           # 自动调整边距避免遮挡

    out = f1.parent / "plot_env_composition_312.png"   # 输出路径：放到第一个文件所在目录
    plt.savefig(out, dpi=300)                    # ✅ 只保存这一次（300dpi）
    plt.close()                                  # 关闭画布，避免交互环境里图像堆积

    print("Done ✅ 输出：", out)                   # 终端提示输出文件路径
    print(label1, c1)                            # 打印样本1标签与计数（便于快速核对）
    print(label2, c2)                            # 打印样本2标签与计数（便于快速核对）


if __name__ == "__main__":
    main()                                       # 直接运行脚本时执行主流程