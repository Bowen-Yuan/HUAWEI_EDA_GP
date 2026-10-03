# 编译
make clean && make

# 运行全部 8 个 ISPD 2005 数据集（默认 ePlace QP 初始点）
./nsp_placer ispd2005/adaptec1/adaptec1
./nsp_placer ispd2005/adaptec2/adaptec2
./nsp_placer ispd2005/adaptec3/adaptec3
./nsp_placer ispd2005/adaptec4/adaptec4
./nsp_placer ispd2005/bigblue1/bigblue1
./nsp_placer ispd2005/bigblue2/bigblue2
./nsp_placer ispd2005/bigblue3/bigblue3
./nsp_placer ispd2005/bigblue4/bigblue4

## 可选优化过程可视化

默认运行不会写任何可视化文件。为默认 NSP 求解器按固定间隔导出快照：

    ./nsp_placer ispd2005/adaptec1/adaptec1 \
      --visualize visualizations/adaptec1 \
      --visualize-every 10

然后离线渲染 GIF（不计入求解器运行时间）：

    python visualize_nsp.py \
      --snapshots visualizations/adaptec1 \
      --benchmark ispd2005/adaptec1/adaptec1 \
      --output visualizations/adaptec1/optimization_animation.gif

`--visualize` 对默认 NSP 生成完整迭代动画；其他可选求解模式至少导出初始与最终帧。

# 可选参数
#   --init <strategy>   初始点策略: eplace|uniform|random|row|scaled|gaussian|quadrant
#   --cwtr              Phase A + CWTR 细化
#   --p0-iters N        Phase 0 迭代数 (默认 150, 0=跳过)
#   --p1-iters N        Phase 1 迭代数 (默认 100)
#   --p2-iters N        Phase 2 迭代数 (默认 150)
#   --p1-lambda X       Phase 1 初始 λ (默认 0.005)
#   --short-edge        仅保留密度梯度主导方向
#   --bo                贝叶斯优化 λ 选择

# 一键运行全部（bash）
for bench in adaptec1 adaptec2 adaptec3 adaptec4 bigblue1 bigblue2 bigblue3 bigblue4; do
  ./nsp_placer "ispd2005/$bench/$bench"
done
