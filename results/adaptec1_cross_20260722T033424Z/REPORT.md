# adaptec1 双向交叉求解实验

所有 HPWL 与 overflow 均由 `adaptive_pareto_soft` 的统一精确评估器复算。

| 实验 | 阶段 | HPWL | strict overflow (%) | legacy overflow (%) |
|---|---:|---:|---:|---:|
| Adaptive init → Alg | solver input | 2179986600.691 | 0.509763 | 0.294603 |
| Adaptive init → Alg | pre-legalize | 50946063.000 | 55.017078 | 51.828272 |
| Adaptive init → Alg | final | 998520928.000 | 41.481313 | 34.494631 |
| a1 default → Adaptive | raw file | 104924229.000 | 0.000000 | 0.000000 |
| a1 default → Adaptive | solver input | 96669138.000 | 57.328821 | 57.328821 |
| a1 default → Adaptive | final | 100654162.000 | 57.091429 | 57.091429 |

原始数据：`D:\codex_project\HUAWEI_EDA\combine\results\adaptec1_cross_20260722T033424Z\summary.json`
