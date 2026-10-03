"""Generate the Chinese PDF report for the a1 RBSM reproduction."""

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "pdf" / "a1_rbsm_reproduction_report.pdf"


def p(text: str, style):
    return Paragraph(text, style)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D9E2F3"))
    canvas.line(1.7 * cm, 1.25 * cm, A4[0] - 1.7 * cm, 1.25 * cm)
    canvas.setFont("STSong-Light", 8)
    canvas.setFillColor(colors.HexColor("#667085"))
    canvas.drawString(1.7 * cm, 0.82 * cm, "a1 / adaptec1  RBSM 全局布局复现报告")
    canvas.drawRightString(A4[0] - 1.7 * cm, 0.82 * cm, f"第 {doc.page} 页")
    canvas.restoreState()


def table(rows, widths):
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "STSong-Light"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.2),
        ("LEADING", (0, 0), (-1, -1), 12),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17365D")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B8C7D9")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F5F8FC")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return t


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
    styles = getSampleStyleSheet()
    title = ParagraphStyle("title", parent=styles["Title"], fontName="STSong-Light", fontSize=24,
                           leading=34, alignment=TA_CENTER, textColor=colors.HexColor("#17365D"), spaceAfter=12)
    subtitle = ParagraphStyle("subtitle", parent=styles["Normal"], fontName="STSong-Light", fontSize=11,
                              leading=18, alignment=TA_CENTER, textColor=colors.HexColor("#475467"))
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontName="STSong-Light", fontSize=16,
                        leading=24, textColor=colors.HexColor("#17365D"), spaceBefore=8, spaceAfter=8)
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontName="STSong-Light", fontSize=11.5,
                        leading=18, textColor=colors.HexColor("#1F4E79"), spaceBefore=8, spaceAfter=5)
    body = ParagraphStyle("body", parent=styles["BodyText"], fontName="STSong-Light", fontSize=9.5,
                          leading=16, spaceAfter=6, alignment=TA_LEFT)
    small = ParagraphStyle("small", parent=body, fontSize=8.2, leading=12, textColor=colors.HexColor("#475467"))
    callout = ParagraphStyle("callout", parent=body, backColor=colors.HexColor("#EAF2F8"),
                             borderColor=colors.HexColor("#9CC2E5"), borderWidth=0.6, borderPadding=8,
                             leading=16, spaceBefore=5, spaceAfter=10)

    doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=1.7 * cm, rightMargin=1.7 * cm,
                            topMargin=1.55 * cm, bottomMargin=1.65 * cm, title="a1 RBSM 全局布局复现报告")
    s = []
    s += [Spacer(1, 2.4 * cm), p("a1 / adaptec1 全局布局算法复现报告", title),
          p("论文：An efficient stochastic subgradient method for the global placement problem in very large-scale integration circuits", subtitle),
          Spacer(1, 0.5 * cm),
          p("范围：独立复现 RBSM（随机批量次梯度）全局布局算法；在 a1 场景运行，并以论文 Table 5 的指标进行核对。", subtitle),
          Spacer(1, 1.0 * cm),
          p("结论摘要", h1),
          p("代码已实现为可运行的 <b>rbsm_place</b> 包，覆盖 Bookshelf 解析、原始 HPWL、边界约束、局部重叠惩罚、随机批量、空间网格候选对、迭代日志与 .pl 输出。a1 上可以达到甚至优于论文的 HPWL 数值，但未能同时达到密度溢出目标；因此本次结论是：<b>算法流程已复现并可运行，论文三项联合指标尚未复现成功</b>。", callout),
          p("最关键的证据：在三次完整的密度优先运行中，最优密度溢出仍约为 53.17%，明显高于论文的 26.8%；而 HPWL 稳定在 5.59e7 左右。", body),
          Spacer(1, 0.55 * cm),
          p("报告日期：2026-07-17　　运行环境：PyTorch 2.7.0+cu118，RTX 3050 Laptop GPU（4 GB）", small),
          PageBreak()]

    s += [p("1. 目标、数据与评测口径", h1),
          p("目标不是只获得一个较小的线长，而是尝试同时复现论文 Table 5 在 adaptec1 上报告的：HPWL 5.05e7、density overflow 26.8%、pre-legalization overlap 10.36%、85 个 epoch。", body),
          p("a1 输入规模：211,447 个对象（210,904 个可移动、543 个固定）、221,142 条网络、944,053 个引脚；芯片核心区域为 (459,459)-(11151,11139)。输入 .pl 的初始 HPWL 为 104,924,229。", body),
          p("本实现使用 pin offset 计算原始 HPWL；密度溢出定义为各 bin 超出可用行容量的面积占可移动总面积的比例，固定单元作为容量阻塞；重叠以可移动单元两两交叠面积/可移动面积诊断。密度评测在指定分箱下精确计算；高密度 bucket 的两两重叠会采用有界确定性采样，故只作诊断而不作为本次是否达标的决定依据。", body),
          p("2. 复现方案", h1),
          p("按照论文的建模，将目标拆为线长项、边界罚项与二维线性帽函数（hat）表示的局部单元重叠罚项。每轮先做 HPWL/均场子步，再做边界/局部重叠子步；局部对通过空间网格检索，并对网络与单元对做随机批量抽样。", body),
          table([
              [p("模块", small), p("实现决定", small), p("原因", small)],
              [p("线长", small), p("原始 HPWL，含 pin offset；按网络抽样并做重要性加权", small), p("保持与论文全局布局目标一致", small)],
              [p("重叠", small), p("二维 hat 重叠罚项 + 空间网格局部候选对", small), p("避免 O(n²) 枚举", small)],
              [p("随机性", small), p("按近似投影度加权抽样；softmax 温度 50,000；噪声 0.2/k³", small), p("复现论文随机批量与退火思想", small)],
              [p("稳定性", small), p("余弦学习率、稀疏活跃对缓存、最大位移裁剪", small), p("4 GB 显存和大规模 a1 下避免数值爆炸", small)],
              [p("固定对象", small), p("保持固定，仅作为密度容量阻塞", small), p("符合 Bookshelf 全局布局语义", small)],
          ], [2.1*cm, 8.2*cm, 5.0*cm]),
          Spacer(1, 0.3*cm),
          p("实现中记录了论文文本的两处歧义：式 (8) 的文字/分段定义对应 min，但定理中写成 max；算法 3 对 f1、f2 的注释与前文定义互换；更新式给出加号而收敛式是减号。本实现主路径选择与分段定义及最小化方向一致的 min 与减号，并保留审计模式以便对照。", callout),
          PageBreak()]

    s += [p("3. 调参与实验过程", h1),
          p("调参按“先确认实现正确，再确认可稳定跑全量，最后压低密度”的顺序进行。下表数值均来自实际运行日志。warm start 使用了已有布局，因此仅作为问题定位与可行性对照，并非严格的论文随机初始化复现。", body),
          table([
              [p("阶段/配置", small), p("关键设置", small), p("结果", small), p("观察", small)],
              [p("随机初始化 smoke", small), p("1 epoch，1 inner step，随机位置", small), p("HPWL 2.214e9；density 16.71%（64 bins）；overlap 16.26%", small), p("CUDA 路径和梯度可执行，但离可用线长很远", small)],
              [p("输入 .pl 诊断", small), p("从 Bookshelf 原始位置开始", small), p("5 epoch 后 HPWL 1.0546e8；重叠约 4.85e6", small), p("大量单元精确重合时 hat 次梯度可取零，局部斥力难以启动", small)],
              [p("线长优先 warm start", small), p("alpha=5，85 epoch，25 inner，lr=5e-4，batch=2%", small), p("HPWL 4.7610e7；density 57.24%", small), p("线长优于论文，但均场设置进一步恶化密度", small)],
              [p("密度优先短跑", small), p("alpha=0，15 epoch，10 inner，local pairs=200k", small), p("HPWL 5.5934e7；density 53.17%", small), p("密度很快停在约 53%，不是迭代次数不足", small)],
              [p("密度优先完整三 seed", small), p("alpha=0，85 epoch，10 inner，lr=1e-3，gamma=1000，batch=1%，200k pairs", small), p("seed 2026/27/28：density 53.1744% / 53.1743% / 53.1777%", small), p("结果高度一致，说明是系统性建模/实现差异而非随机种子问题", small)],
          ], [2.75*cm, 4.6*cm, 4.1*cm, 3.85*cm]),
          Spacer(1, 0.4*cm),
          p("已尝试的主要旋钮：均场权重 alpha（5 → 0）、学习率（5e-4 → 1e-3）、inner steps（25 → 10）、批量占比（2% → 1%）、局部对上限（5 万 → 20 万）、pair refresh（5 → 2）、最大位移（4 → 10）、迭代次数（15 → 85）、随机 seed（2026/2027/2028）。其中增大局部对覆盖和去除均场有助于避免进一步恶化，但都未突破约 53.17% 的密度平台。", body),
          PageBreak()]

    s += [p("4. 结果与论文指标核对", h1),
          table([
              [p("方案", small), p("HPWL", small), p("density overflow", small), p("overlap", small), p("判定", small)],
              [p("论文 Table 5", small), p("5.05e7", small), p("26.8%", small), p("10.36%", small), p("目标", small)],
              [p("线长优先 warm start", small), p("4.7610e7", small), p("57.24%", small), p("5.34e6（采样估计）", small), p("仅 HPWL 达标", small)],
              [p("密度优先，seed 2026", small), p("5.5872e7", small), p("53.1744%", small), p("3.08e5（采样估计）", small), p("未达标", small)],
              [p("密度优先，seed 2027", small), p("5.5909e7", small), p("53.1743%", small), p("3.10e5（采样估计）", small), p("未达标", small)],
              [p("密度优先，seed 2028", small), p("5.6095e7", small), p("53.1777%", small), p("2.27e5（采样估计）", small), p("未达标", small)],
          ], [3.5*cm, 2.7*cm, 3.5*cm, 3.8*cm, 2.8*cm]),
          Spacer(1, 0.4*cm),
          p("完整密度优先三次运行的中位数为 HPWL 5.591e7、density overflow 53.1744%。与论文相比，HPWL 高约 10.7%，密度溢出高约 26.4 个百分点。这一差距足以否定“达到论文联合指标”的说法。", callout),
          p("另一方面，线长优先配置达到 4.761e7 的 HPWL，说明 Bookshelf 解析、引脚偏移线长计算、GPU 执行与输出布局链路是有效的；失败集中在全局密度控制，而不是程序完全无法优化。", body),
          p("正确性检查：对精确 HPWL、二维重叠梯度与受限步长共 3 个单元测试均通过；三份完整 solution.pl 已输出，可供后续布局可视化或接入外部评测器。", body),
          PageBreak()]

    s += [p("5. 我认为问题主要在哪里", h1),
          p("<b>第一，局部两两重叠罚项与全局 bin density 指标并不等价。</b>论文优化项是局部 hat 斥力，Table 5 却报告 bin 容量溢出；二者相关但没有直接保证。当前实验恰好呈现出这个缺口：局部重叠诊断下降，并不意味着每个 bin 的容量溢出同步下降。", body),
          p("<b>第二，候选对截断限制了高密区域的有效斥力。</b>a1 的潜在近邻对非常多，为适配显存与时间，单次最多处理 20 万对并周期性刷新。密集区域中未被选中的单元没有获得足够的排斥梯度，容易形成稳定的高密度平台。", body),
          p("<b>第三，初始化对非凸问题过于关键。</b>随机初始化的线长极大；输入 .pl 又会出现单元精确重合、hat 次梯度为零的退化情况；warm start 能获得好线长，却不属于论文严格的随机起点。论文未给出可直接复现的 adaptec1 初始位置、随机流、预扩散或多级初始化细节。", body),
          p("<b>第四，论文遗漏了决定可比性的工程细节，并存在公式歧义。</b>包括密度 bin 的大小/容量定义、pair 采样温度和预算、gamma 调度、重叠评测器、固定块处理、初始化与停止准则。4 GB 显存也迫使本次使用 candidate 截断、较小批量和位移裁剪，虽保持了算法结构，但无法等同于论文作者的运行配置。", body),
          p("因此，目前更合理的表述是：<b>尚不能判断论文结论错误；现有公开描述不足以让独立实现可靠复现其 a1 的三项联合数字。</b>", callout),
          p("6. 下一步建议", h1),
          p("(1) 首先向作者/原始代码核对初始布局、bin 容量、重叠定义、pair 预算与 gamma 日程；(2) 用与论文一致的官方评测器复算；(3) 若目标转为工程可用而非严格复刻，在 RBSM 外加入可微的全局电势/密度力或容量约束，并使用多级扩散初始化；(4) 增加高密 bin 的定向候选对覆盖和无截断/分块并行实验，检查 53% 平台是否来自采样预算。", body),
          p("交付物：实现位于 dai-code/rbsm_place；可复现实验报告位于 dai-code/results/a1_reproduction_report.md；本 PDF 为该过程的可阅读总结。", small)]

    doc.build(s, onFirstPage=footer, onLaterPages=footer)
    print(OUT)


if __name__ == "__main__":
    main()
