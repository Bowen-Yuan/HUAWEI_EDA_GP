import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_multiple_blocks_per_subplot(block_file, npy_files, titles, output_filename,
                                     Ws, Hs, wirelength, overlap):
    """
    每个子图使用各自的 (W, H) 轴范围；不共享 y 轴。
    """
    blocks = np.load(block_file)

    # 不共享 y，因为每个 H 不同
    fig, axes = plt.subplots(
        1, len(npy_files),
        sharey=False,
        figsize=(len(npy_files) * 5, 6)   # 每张图宽 3.2, 高 3.5，整体更紧凑
    )
    plt.subplots_adjust(wspace=1.1)

    # 兼容 len(npy_files)==1 的情况
    if len(npy_files) == 1:
        axes = [axes]

    for idx, (npy_file, title) in enumerate(zip(npy_files, titles)):
        pos = np.load(npy_file)
        ax = axes[idx]

        # 画矩形
        for (bw, bh), (xc, yc) in zip(blocks, pos):
            ax.add_patch(
                patches.Rectangle((xc - bw/2, yc - bh/2), bw, bh,
                                  linewidth=1, edgecolor='blue',
                                  facecolor='skyblue', alpha=0.5)
            )

        # 每个子图用各自的 W/H
        ax.set_xlim(0, Ws[idx])
        ax.set_ylim(0, Hs[idx])
        ax.set_aspect('equal', adjustable='box')

        # 标题 & 轴标签
        ax.set_title(f'HPWL = {wirelength[idx]:.1f}, Overlap = {overlap[idx]:.1f}',
                     fontsize=16, color='black', fontweight='bold')
        ax.set_xlabel('Width', fontsize=15)
        if idx == 0:
            ax.set_ylabel('Height', fontsize=15)
        ax.tick_params(axis='both', labelsize=16)

        # 在 x 轴下方放 (a)(b)(c)(d)
        ax.text(0.5, -0.2, titles[idx], fontsize=20,
                ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(output_filename, format='eps', bbox_inches='tight')  # 建议导出 png/pdf
    plt.show()


# ------------ 示例调用 ------------
block_file = './result/n200/block.npy'
npy_files = ['./result/wh/legal800.npy',
             './result/wh/legal700.npy',
             './result/wh/legal600.npy',
             './result/wh/legal500.npy']

titles = ['(a) W,H = 800', '(b) W,H = 700', '(c) W,H = 600', '(d) W,H = 500']

# 每个子图各自的 W/H
Ws = [800, 700, 600, 500]
Hs = [800, 700, 600, 500]

wirelength = [557147.4, 566231.1, 562323.9, 652456.0]  # 举例
overlap    = [2246.2,   2272.22,   2361.61,   4692.72]    # 举例
wirelength1 = [560782.6, 569241.1, 564395.4, 658194.1]  # 举例
overlap1    = [0.0,   0.0,   0.0,   0.0]    # 举例
draw_multiple_blocks_per_subplot(block_file, npy_files, titles,
                                 'wh.eps',  # 建议导出 png/pdf
                                 Ws, Hs, wirelength, overlap)