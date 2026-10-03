# 编译
make clean && make

> 坐标与评估约定：Alg 内部使用 cell centre；Bookshelf `.pl` 输入/输出自动在
> centre 与 lower-left 之间转换。正式 HPWL/overflow 仍须用仓库的官方脚本复算。

# 运行全部 8 个 ISPD 2005 数据集（默认 ePlace QP 初始点）
./nsp_placer ispd2005/adaptec1/adaptec1
./nsp_placer ispd2005/adaptec2/adaptec2
./nsp_placer ispd2005/adaptec3/adaptec3
./nsp_placer ispd2005/adaptec4/adaptec4
./nsp_placer ispd2005/bigblue1/bigblue1
./nsp_placer ispd2005/bigblue2/bigblue2
./nsp_placer ispd2005/bigblue3/bigblue3
./nsp_placer ispd2005/bigblue4/bigblue4

# 可选参数
#   --init <strategy>   初始点策略: eplace|uniform|random|row|scaled|gaussian|quadrant
#   --cwtr              Phase A + CWTR 细化
#   --p0-iters N        Phase 0 迭代数 (默认 150, 0=跳过)
#   --p1-iters N        Phase 1 迭代数 (默认 100)
#   --p2-iters N        Phase 2 迭代数 (默认 150)
#   --p1-lambda X       Phase 1 初始 λ (默认 0.005)
#   --density-target R  ISPD density target，默认 0.60；不要使用旧版自动推导的 0.75×全局密度
#   --short-edge        仅保留密度梯度主导方向
#   --bo                贝叶斯优化 λ 选择

# 一键运行全部（bash）
for bench in adaptec1 adaptec2 adaptec3 adaptec4 bigblue1 bigblue2 bigblue3 bigblue4; do
  ./nsp_placer "ispd2005/$bench/$bench"
done
