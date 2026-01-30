# ASH xTB 显式分步工作流：经验笔记（Opt / NEB(CINEB) / TSOpt / NumFreq / IRC）

> 目的：把一次“单反应验证”的关键调用方式、输出文件、常见坑（尤其是 NumPy 2.x 与 KNARR/NEB 的兼容性）记录下来，后续不用再逐页翻文档。

## 1. 目标工作流（显式分步 A）

对应你原先 `do_orca_validation.py` 的概念流程：

1) 反应物优化（IM opt）  
2) 产物优化（IM opt）  
3) CI-NEB / CINEB（带 `TS_guess_file` 插入）  
4) TS 优化（OptTS / TSOpt）  
5) 频率（NumFreq / AnFreq）  
6) IRC（基于 TS + Hessian）
7) IRC 两端点提取 + 端点再优化（forward/backward）  
8) IRC 两端点再做频率（应无虚频；若有则终止）
9) IRC 两端点与原始 R/P 做交叉同构性判断（canonical SMILES）

ASH 里主要用到：

- `Optimizer(...)`：几何优化 / TSOpt / IRC（geomeTRIC）
- `NEB(...)`：CI-NEB / CINEB（KNARR）
- `NumFreq(...)`：数值频率（可产生 Hessian）
- `xTBTheory(...)`：xTB 量化理论（runmode=`inputfile` 或 `library`）

文档入口（你给的）：

- NEB/CINEB：`https://ash.readthedocs.io/en/latest/neb.html`
- Optimizer/TSOpt/IRC：`https://ash.readthedocs.io/en/latest/Geometry-optimization.html`
- NumFreq：`https://ash.readthedocs.io/en/latest/module_freq.html`

## 2. 可复现的测试脚本

仓库根目录新增了 smoke test：

- `ash_xtb_explicit_workflow_smoketest.py`

默认读取：

- `xyz/reactant.xyz`
- `xyz/product.xyz`
- `xyz/ts_ini.xyz`（用于 `NEB(..., TS_guess_file=...)`）

运行（在 conda env `ash` 下）：

```bash
conda run -n ash python -u ash_xtb_explicit_workflow_smoketest.py
```

运行结束会生成一个带时间戳的目录，例如：

- `ash_xtb_smoketest_runs/run_20260129_173431/`
- 关键汇总：`ash_xtb_smoketest_runs/run_20260129_173431/workflow_summary.json`

## 3. 本次示例的运行时间（量级）

在示例体系（12 atoms，GFN2-xTB，1 核，NEB images=7(总9张)，NEB maxiter=100，NumFreq npoint=1）上：

- **一次完整跑完（含 IRC）约 30 秒量级**（实际强依赖机器/负载/NEB 是否收敛/IRC 步数等）。

经验上耗时主要来自：

- NEB（图像数 × 迭代数 × 每步 E+Grad）
- IRC（步数 × 每步 E+Grad）
- NumFreq（位移点数 npoint 与原子数相关；npoint=1 更快）

## 4. 每步常见输出文件（方便你后续解析/对齐）

以 `run_20260129_173431` 为例：

### 4.1 Opt（reactant / product）

目录：`01_opt_reactant/`、`02_opt_product/`

常见文件：

- `Fragment-optimized.xyz`：优化后结构（脚本也会复制到 run 根目录 `R_opt.xyz` / `P_opt.xyz`）
- `ASH_Optimizer.result`：ASH 的结果对象（JSON）
- `geometric_OPTtraj.log` / `geometric_OPTtraj_optim.xyz`：geomeTRIC 的轨迹与日志
- `xtb_.out` / `xtb_.engrad`：xTB 输出与梯度

### 4.2 NEB（CI-NEB / CINEB）

目录：`03_neb/`

常见文件：

- `TSguess.xyz`：ASH 将 `TS_guess_file` 复制/写出的插入构型（KNARR 会用它做插值/IDPP）
- `initial_guess_path.xyz`：初始路径
- `idpp_optimization.xyz`、`idpp.interp`：IDPP 路径优化过程与插值文件
- `knarr_optimization.xyz`、`knarr_current.xyz`、`knarr_last_iter.xyz`：NEB 迭代轨迹与最后一次迭代路径
- `image_0/` ... `image_8/`：每个图像的单点/梯度（xTB 输入输出）
- `ASH_NEB.result`：NEB 结果（成功时会有 `saddlepoint_fragment`、`MEP_energies_dict` 等；失败时 label 类似 `"NEB-CI calc (fail)"`）

重要说明：

- 脚本现在采用“筛选/终止”模式：**NEB 未收敛时会直接终止并记录**（不再用 `knarr_last_iter.xyz` 的最高能量 image 作为 fallback 继续跑）。

### 4.3 TSOpt

目录：`04_tsopt/`

常见文件：

- `Fragment-optimized.xyz`：TSOpt 后结构（脚本复制到 run 根目录 `TS_opt.xyz`）
- `Hessian_from_xtb` / `hessian` / `vibspectrum` / `xtbhess.xyz`：因为 `TSOpt` 使用了 `hessian="xtb"`（xTB 直接算 Hessian）  
- `geometric_OPTtraj.log`：TS 优化日志

### 4.4 NumFreq（数值频率）

目录：`05_numfreq/`

常见文件：

- `Numfreq_dir/Hessian`：数值 Hessian
- `Numfreq_dir/orcahessfile.hess`：ORCA 风格 Hessian 文件（ASH 为可视化/兼容写出的）
- `ASH_NumFreq.result`：包含 `frequencies`、`hessian` 等

这次示例结果（见 `workflow_summary.json`）：

- `imaginary_freq_count = 1`
- 首个虚频约 `-472.9 cm^-1`

### 4.5 IRC

目录：`06_irc/`

常见文件：

- `geometric_OPTtraj_irc.xyz`：IRC 轨迹（关键）
- `Hessian_np`：脚本把 `NumFreq` 返回的 Hessian（numpy array）写成 geometric 可读的 Hessian 文件
- `Fragment-optimized.xyz`：IRC 结束时的最终结构（注意：IRC 会走路径，最终点不一定是你要的两端，需要你后处理/分段取两端）

脚本新增的端点处理（在 run 根目录）：

- `IRC_forward_end.xyz` / `IRC_backward_end.xyz`：从 `geometric_OPTtraj_irc.xyz` 直接取**第一帧/最后一帧**作为两端点（并用 log 的 forward/back step 做一致性校验，必要时交换 forward/back 标签）
- `IRC_forward_opt.xyz` / `IRC_backward_opt.xyz`：对两端点分别再优化后的结构（目录 `07_opt_irc_forward/`、`08_opt_irc_backward/`）
- `07a_freq_irc_forward/`、`08a_freq_irc_backward/`：对再优化后的两端点做频率；**要求无虚频**（若出现虚频脚本会终止并在 `workflow_summary.json` 写明原因）

重要坑点（这次踩到的）：

- `geometric_OPTtraj_irc.xyz` 的 comment 行里会有 `Iteration N`，但 **N 可能重复**（例如先 `48..0`，再 `0..72`）。因此不能用 `{iteration: block}` 这种 dict 去取端点，否则会被后面的重复 iteration 覆盖，导致 forward/back 端点取错甚至取成同一个。
- 当前脚本的端点抽取策略改为“**第一帧/最后一帧**”，这是对这种重复 iteration 情况最稳健的。

### 4.6 端点交叉匹配（SMILES / 同构性）

目录：`09_endpoint_match/`

常见文件：

- `endpoint_match.json`：使用 OpenBabel/pybel 将 XYZ 转为 canonical SMILES（含立体），并做交叉匹配，输出：
  - `endpoint_match`: `2-end match` / `1-end match` / `No match` / `Error`
  - `matches`: (trueR==ircR, trueP==ircP, trueR==ircP, trueP==ircR)
  - `rxn_status`: `Chemical reaction` / `Conformational change`

## 4.7 结果持久化（PKL 数据结构）

`workflow_summary.json` 适合"人读/排错"，脚本同时会生成 `ash_xtb_result.pkl`（只在全流程成功时写出），用于程序化取用。

### 数据来源说明

**重要**：所有物种（reactant/product/ts）的数据**统一来自对应的 NumFreq 计算目录**，确保数值一致性。

| 物种 | 数据来源目录 | 说明 |
|------|-------------|------|
| `reactant` | `08a_freq_irc_backward/Numfreq_dir/` | IRC 后向端点优化后做频率 |
| `product` | `07a_freq_irc_forward/Numfreq_dir/` | IRC 前向端点优化后做频率 |
| `ts` | `05_numfreq/Numfreq_dir/` | 过渡态优化后做频率 |
| `irc` | `06_irc/geometric_OPTtraj_irc.xyz` | IRC 轨迹（多帧） |

**注意**：
- `reactant` 实际对应 `irc_backward`（IRC 反方向端点）
- `product` 实际对应 `irc_forward`（IRC 正方向端点）
- 这样设计是因为 IRC 验证了 TS 正确连接了两个端点

### PKL 数据结构

```python
{
    # ========== 顶层元数据 ==========
    "units": {
        "coords": "Å",
        "energy": "Eh",
        "gradient": "Eh/Å",
        "hessian": "Eh/Å²",
        "frequencies": "cm⁻¹",
        "thermo": "Eh",
        "temperature": "298.15 K",
    },
    "system": {
        "elements": ["C", "C", "C", "O", "H", ...],  # 原子元素列表
        "charge": 0,                                  # 电荷
        "mult": 1,                                    # 自旋多重度
        "natoms": 12,                                 # 原子数
    },

    # ========== 反应物（来自 IRC backward 端点 NumFreq） ==========
    "reactant": {
        "coords_A": np.ndarray,           # (N, 3) 坐标，单位 Å
        "energy_Eh": float,               # 电子能量，单位 Eh
        "gradient_Eh_per_A": np.ndarray,  # (N, 3) 梯度，单位 Eh/Å
        "hessian_Eh_per_A2": np.ndarray,  # (3N, 3N) Hessian 矩阵，单位 Eh/Å²
        "frequencies_cm-1": np.ndarray,   # (3N,) 频率，单位 cm⁻¹（应无虚频）
        "thermochemistry": {
            "ZPVE_Eh": float,             # 零点振动能矫正，单位 Eh
            "Hcorr_Eh": float,            # 焓矫正，单位 Eh
            "Gcorr_Eh": float,            # 吉布斯自由能矫正，单位 Eh
        },
    },

    # ========== 产物（来自 IRC forward 端点 NumFreq） ==========
    "product": {
        # 结构同 reactant
    },

    # ========== 过渡态（来自 TS NumFreq） ==========
    "ts": {
        # 结构同 reactant
        # 注意：frequencies_cm-1 应包含 1 个虚频（< 0）
    },

    # ========== IRC 轨迹 ==========
    "irc": {
        "coords_A": np.ndarray,           # (n_frames, N, 3) IRC 轨迹坐标
        "energy_Eh": np.ndarray,          # (n_frames,) 每帧能量
    },

    # ========== IRC 端点（向后兼容，与 reactant/product 相同） ==========
    "irc_forward": {
        # 与 product 完全相同
    },
    "irc_backward": {
        # 与 reactant 完全相同
    },
}
```

### 读取示例

```python
import pickle
import numpy as np

with open('ash_xtb_smoketest_runs/run_xxx/ash_xtb_result.pkl', 'rb') as f:
    data = pickle.load(f)

# 获取系统信息
elements = data['system']['elements']
charge = data['system']['charge']

# 获取过渡态能量和热力学矫正
ts_energy = data['ts']['energy_Eh']
ts_gibbs = ts_energy + data['ts']['thermochemistry']['Gcorr_Eh']

# 获取反应能和能垒
r_energy = data['reactant']['energy_Eh']
p_energy = data['product']['energy_Eh']
barrier = data['ts']['energy_Eh'] - r_energy
reaction_energy = p_energy - r_energy

# 验证 TS 有且仅有 1 个虚频
imag_freqs = data['ts']['frequencies_cm-1'][data['ts']['frequencies_cm-1'] < 0]
assert len(imag_freqs) == 1

# 检查端点无虚频
assert np.all(data['reactant']['frequencies_cm-1'] >= 0)
assert np.all(data['product']['frequencies_cm-1'] >= 0)

# 获取 IRC 轨迹
irc_coords = data['irc']['coords_A']  # (n_frames, N, 3)
irc_energies = data['irc']['energy_Eh']  # (n_frames,)
```

### 文件位置

- 文件名：`ash_xtb_result.pkl`
- 位置：`ash_xtb_smoketest_runs/run_YYYYMMDD_HHMMSS/ash_xtb_result.pkl`
- 生成条件：**仅在全流程成功时生成**，任一步骤失败则不生成（避免不完整数据）

## 5. 为什么 “不改代码” 会在你当前环境出问题（NumPy 2.x）

你当前 conda env `ash` 里 NumPy 是：

- `numpy 2.4.1`

在这个版本上，NumPy 对“(1,) 数组 → 标量”的隐式转换已经变得非常严格，典型表现是：

- `float(np.array([1.0]))` 会报：`TypeError: only 0-dimensional arrays can be converted to Python scalars`
- `a[0] = np.array([1.0])` 会报：`ValueError: setting an array element with a sequence.`

而 KNARR/NEB 里很多地方历史上用的是列向量 `(N,1)`（或 `array([x])`）并把它当标量/1D 来打印/赋值/写文件，于是在 `numpy 2.4.1` 下会直接跑崩。

这也是为什么虽然 ASH 逻辑上支持 `TS_guess_file`，但在你这个环境里 **NEB 路径生成/IDPP/打印/写轨迹** 很容易触发异常。

补充：我们另外建的 env `ash-py310` 里 NumPy 是 `2.2.6`，同样的隐式转换**还没有硬报错**，只是会出现：

- `DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future.`

因此在 `ash-py310 (numpy 2.2.6)` 下，“不改代码”目前可以跑通，但未来升级 NumPy 仍有再次触发硬错误的风险。

## 6. 不修改代码的替代方案：降级 NumPy（推荐优先尝试）

最直接的“零代码改动”方案有两类：

1) **固定到一个目前可工作的 NumPy 版本**（例如我们验证过的 `numpy 2.2.6`）  
2) **把 NumPy 固定到 `<2`**（例如 1.26.x），通常更稳，也更不容易踩到 2.x 的行为变化

建议做法（两种思路）：

### 6.1 在现有 env 里降级（简单，但可能触发依赖重解）

```bash
conda install -n ash "numpy<2"
python -c "import numpy as np; print(np.__version__)"
```

如果求稳，可以直接 pin：

```bash
conda install -n ash "numpy=1.26.*"
```

如果你想走“固定到可工作版本”的路线，也可以 pin 到我们已验证的版本（示例）：

```bash
conda install -n ash "numpy=2.2.6"
```

### 6.2 新建一个“专用 env”（更稳，不影响你现有 env）

```bash
conda create -n ash_np1 python=3.12 "numpy<2"
# 然后在该环境里安装 ash / geometric / xtb 等依赖
```

> 说明：我在这个沙盒里无法联网替你实际执行 `conda install`/`conda create` 的求解与下载，但从报错特征与 NumPy 2.4.1 的行为来看，降级是最合理的无代码 workaround。

## 7. 如果你允许改代码：我做过的兼容性补丁（仅记录）

我们已经做过一次 A 验证：把补丁撤掉（恢复到 ASH/KNARR 原版代码）后，在 `ash-py310 (numpy 2.2.6)` 下 smoketest 依然可以完整跑通；说明你现在这套环境确实可以先走“零代码改动”路线。

但为了在 `numpy 2.4.1` 下也能直接跑通 NEB/IDPP/输出，我之前做过一组“把 (N,1) flatten、打印/写文件处强制 float/int”的兼容性补丁，涉及：

- `ash/ash/knarr/KNARRatom/utilities.py`
- `ash/ash/knarr/KNARRcalculator/utilities.py`
- `ash/ash/knarr/KNARRio/io.py`
- `ash/ash/knarr/KNARRio/output_print.py`
- `ash/ash/knarr/KNARRjobs/path.py`
- `ash/ash/knarr/KNARRjobs/utilities.py`
- `ash/ash/knarr/KNARRjobs/neb.py`

当前这些改动**没有保留在工作区**，而是放在 `ash/` 子仓库的 git stash 里，方便需要时再取用：

- 查看 stash：`git -C ash stash list`
- 查看补丁内容：`git -C ash stash show -p stash@{0}`
- 应用补丁（把改动恢复到工作区）：`git -C ash stash apply stash@{0}`
- 如果不需要补丁、想保持原版代码：`git -C ash restore .`（或直接不 apply）

使用建议：

- 如果你长期固定在 `numpy 2.2.6`（或 `<2`），可以先不管这个补丁。
- 如果你之后需要回到 `env ash (numpy 2.4.1)` 跑 NEB/IDPP，建议要么 apply 这个补丁，要么把 NumPy 降级/固定版本。

## 8. 下一步你可能需要的改进点

- **让 NEB 真正收敛**：提高 `--neb-maxiter`，或调 `images`、收敛阈值等（否则 TS_guess 只能算“最高能量 image”的近似）。
- **TS 有效性判断**：现在只记录虚频数；后续可以把 TS 的虚频对应的振动模式导出并对关键键长变化做判断。
- **端点判定更严格**：目前端点验证是用 canonical SMILES 交叉匹配；后续可以考虑在 SMILES 失败时增加更稳健的备选方案（如图同构/反应中心检查等）。
