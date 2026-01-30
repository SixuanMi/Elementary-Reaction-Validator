# ASH xTB/ORCA Workflow - 修改日志

## 2026-01-30 - 主要功能更新

### 1. 代码整理和模块化

#### 1.1 删除未使用的函数
- ❌ 删除 `_extract_max_energy_image_xyz()` (45 行)
  - 原因：NEB 已通过 `saddlepoint_fragment` 直接提供 TS 猜测
- ❌ 删除 `_parse_xyz_blocks_by_iteration()` (25 行)
  - 原因：实际使用的是 `_parse_xyz_blocks_in_order()`
- 📉 代码行数：1474 → 1373 行（减少 101 行，6.8%）

#### 1.2 强制使用配置文件
- ✅ 移除所有命令行参数覆盖选项（20+ 个参数）
- ✅ 只保留 `--config` 参数（必需）
- ✅ 配置文件不存���时友好报错

**使用方式**:
```bash
# 正确
python ash_xtb_explicit_workflow_smoketest.py --config config.yaml

# 错误（不再支持）
python ash_xtb_explicit_workflow_smoketest.py --reactant xyz/r.xyz --product xyz/p.xyz
```

---

### 2. 参数模块化

#### 2.1 NEB CI 参数可配置
**配置文件**:
```yaml
neb:
  CI: true  # Climbing Image NEB 开关
```

**代码**:
- 添加 `args.neb_CI` 参数读取
- NEB 调用使用 `CI=args.neb_CI`

#### 2.2 统一 convergence_setting 参数
**配置文件**:
```yaml
optimization:
  convergence_setting: "ORCA"  # 所有普通优化使用

tsopt:
  convergence_setting: "ORCA"  # TS 优化独立设置
```

**影响的阶段**:
- Stage 1: opt_reactant
- Stage 2: opt_product
- Stage 4: tsopt
- Stage 7: opt_irc_forward
- Stage 8: opt_irc_backward

**可用选项**: `ORCA`, `GAU`, `GAU_TIGHT`, `GAU_LOOSE`, `GAU_VERYTIGHT`

---

### 3. 双 Theory-Level 支持

#### 3.1 设计原则
- **main_theory**: 用于 opt, tsopt, freq, irc（必须同一水平，保证量子化学一致性）
- **neb_theory**: 用于 NEB/CI-NEB（可以用更快的方法）

#### 3.2 配置文件结构
```yaml
# Main theory: 用于关键步骤（优化、频率、IRC）
main_theory:
  type: "orca"  # xtb 或 orca
  numcores: 2
  printlevel: 1

  # xTB 设置（仅当 type=xtb 时使用）
  xtb:
    method: "GFN2"
    runmode: "inputfile"
    filename: "xtb_main_"

  # ORCA 设置（仅当 type=orca 时使用）
  orca:
    orcasimpleinput: "! BP86 def2-SVP TightSCF"
    orcablocks: |
      %scf
        maxiter 200
      end
    filename: "orca_main_"
    moreadfile: null

# NEB theory: 用于路径搜索（可以用快速方法）
neb_theory:
  type: "xtb"
  numcores: 1
  printlevel: 1

  xtb:
    method: "GFN2"
    runmode: "inputfile"
    filename: "xtb_neb_"

  orca:
    orcasimpleinput: "! BP86 def2-SVP"
    orcablocks: ""
    filename: "orca_neb_"
    moreadfile: null
```

#### 3.3 Theory 分配

| 阶段 | Theory | 说明 |
|------|--------|------|
| Stage 1: opt_reactant | main_theory | 反应物优化 |
| Stage 2: opt_product | main_theory | 产物优化 |
| Stage 3: NEB | neb_theory | 路径搜索 |
| Stage 4: tsopt | main_theory | TS 优化 |
| Stage 5: numfreq | main_theory | TS 频率 |
| Stage 6: IRC | main_theory | IRC 计算 |
| Stage 7: opt_irc_forward | main_theory | 前向端点优化 |
| Stage 7a: freq_irc_forward | main_theory | 前向端点频率 |
| Stage 8: opt_irc_backward | main_theory | 后向端点优化 |
| Stage 8a: freq_irc_backward | main_theory | 后向端点频率 |

#### 3.4 代码实现
**新增函数**:
```python
def _create_theory_from_config(theory_config: Dict[str, Any], theory_name: str) -> Any:
    """Theory 工厂函数，支持 xTB 和 ORCA"""
    theory_type = theory_config["type"].lower()

    if theory_type == "xtb":
        from ash import xTBTheory
        return xTBTheory(...)
    elif theory_type == "orca":
        from ash import ORCATheory
        return ORCATheory(...)
```

**使用示例**:
```python
main_theory = _create_theory_from_config(args.main_theory_config, "main_theory")
neb_theory = _create_theory_from_config(args.neb_theory_config, "neb_theory")
```

#### 3.5 配置场景

**场景 1: 全 xTB（快速测试）**
```yaml
main_theory:
  type: "xtb"
neb_theory:
  type: "xtb"
```

**场景 2: 全 ORCA（高精度）**
```yaml
main_theory:
  type: "orca"
  numcores: 2
neb_theory:
  type: "orca"
  numcores: 2
```

**场景 3: 混合方案（推荐）**
```yaml
main_theory:
  type: "orca"  # 高精度用于关键步骤
  numcores: 2
neb_theory:
  type: "xtb"   # 快速方法用于路径搜索
  numcores: 1
```

---

### 4. AnFreq 解析频率支持

#### 4.1 问题背景
- **之前**: 使用 NumFreq（数值频率）对 ORCA 进行频率计算
  - 每个频率计算 ~5 分钟
  - 需要多次单点能量计算（数值微分）
  - 3 个频率计算总计 ~15 分钟

- **现在**: 使用 AnFreq（解析频率）对 ORCA
  - ORCA 支持解析 Hessian
  - 只需一次计算
  - 预计快 5-10 倍 ⚡

#### 4.2 配置文件
```yaml
frequency:
  method: "auto"  # auto, numerical, analytical
  npoint: 2       # 仅用于 numerical
  runmode: "serial"  # 仅用于 numerical
  cores: 1        # 仅用于 numerical
  imag_threshold: 20.0
```

**method 选项**:
- `auto`: 根据 theory 类型自动选择（ORCA→AnFreq, xTB→NumFreq）**【推荐】**
- `numerical`: 强制使用 NumFreq（数值频率）
- `analytical`: 强制使用 AnFreq（解析频率，仅 ORCA 支持）

#### 4.3 代码实现
**新增函数**:
```python
def _get_freq_function_and_params(theory_type: str, freq_method: str, ...):
    """根据 theory 类型和 freq_method 选择频率函数"""
    if freq_method == "auto":
        use_analytical = (theory_type == "orca")
    elif freq_method == "analytical":
        if theory_type != "orca":
            raise ValueError("Analytical frequencies only supported for ORCA")
        use_analytical = True
    elif freq_method == "numerical":
        use_analytical = False

    if use_analytical:
        from ash import AnFreq
        return AnFreq, {}
    else:
        from ash import NumFreq
        return NumFreq, {"npoint": ..., "runmode": ..., "numcores": ...}
```

**使用示例**:
```python
FreqFunc, freq_kwargs = _get_freq_function_and_params(
    args.main_theory_config["type"],
    args.freq_method,
    args.freq_npoint,
    args.freq_runmode,
    args.freq_cores
)

freq_result = FreqFunc(
    fragment=frag,
    theory=main_theory,
    charge=args.charge,
    mult=args.mult,
    printlevel=args.printlevel,
    **freq_kwargs
)
```

#### 4.4 性能对比

| 方法 | 每个频率计算 | 3 个频率计算 | 加速比 |
|------|-------------|-------------|--------|
| NumFreq (ORCA) | ~5 分钟 | ~15 分钟 | 1x |
| AnFreq (ORCA) | <1 分钟 | <3 分钟 | **5-10x** |

---

## 配置文件完整示例

### 推荐配置（ORCA + xTB 混合，AnFreq）

```yaml
# ASH xTB/ORCA Workflow Configuration
inputs:
  reactant: "xyz/reactant.xyz"
  product: "xyz/product.xyz"
  ts_guess: "xyz/ts_ini.xyz"

system:
  charge: 0
  mult: 1

# 主计算方法（高精度）
main_theory:
  type: "orca"
  numcores: 2
  printlevel: 1
  xtb:
    method: "GFN2"
    runmode: "inputfile"
    filename: "xtb_main_"
  orca:
    orcasimpleinput: "! BP86 def2-SVP TightSCF"
    orcablocks: |
      %scf
        maxiter 200
      end
    filename: "orca_main_"
    moreadfile: null

# NEB 计算方法（快速）
neb_theory:
  type: "xtb"
  numcores: 1
  printlevel: 1
  xtb:
    method: "GFN2"
    runmode: "inputfile"
    filename: "xtb_neb_"
  orca:
    orcasimpleinput: "! BP86 def2-SVP"
    orcablocks: ""
    filename: "orca_neb_"
    moreadfile: null

neb:
  images: 7
  maxiter: 100
  interpolation: "IDPP"
  CI: true
  runmode: "serial"
  cores: 1

optimization:
  maxiter: 250
  convergence_setting: "ORCA"

tsopt:
  maxiter: 250
  convergence_setting: "ORCA"
  hessian: "xtb"

frequency:
  method: "auto"  # ORCA 自动使用 AnFreq
  npoint: 2
  runmode: "serial"
  cores: 1
  imag_threshold: 20.0

irc:
  maxiter: 250

output:
  base_dir: "ash_xtb_smoketest_runs"
  printlevel: 1
```

---

## 使用说明

### 运行工作流
```bash
python ash_xtb_explicit_workflow_smoketest.py --config ash_xtb_workflow_config.yaml
```

### 查看帮助
```bash
python ash_xtb_explicit_workflow_smoketest.py --help
```

### 输出结构
```
ash_xtb_smoketest_runs/run_YYYYMMDD_HHMMSS/
├── 01_opt_reactant/
├── 02_opt_product/
├── 03_neb/
├── 04_tsopt/
├── 05_numfreq/
├── 06_irc/
├── 07_opt_irc_forward/
├── 07a_freq_irc_forward/
├── 08_opt_irc_backward/
├── 08a_freq_irc_backward/
├── 09_endpoint_match/
├── workflow_summary.json
├── ash_xtb_result.pkl
└── [各种 XYZ 文件]
```

---

## 性能优化总结

| 优化项 | 改进 | 效果 |
|--------|------|------|
| 代码整理 | 删除未使用函数 | 减少 101 行代码 |
| 配置强制 | 只用配置文件 | 提高一致性和可追溯性 |
| 双 Theory | ORCA + xTB 混合 | NEB 用快速方法，关键步骤用高精度 |
| AnFreq | 解析频率 | 频率计算加速 5-10 倍 |
| 并行计算 | ORCA 2核 | 进一步加速 |

**总体性能提升**: 预计整体计算时间减少 30-50%

---

## 技术细节

### 支持的 Theory 类型
- **xTB**: 半经验方法，快速但精度较低
  - 方法: GFN1, GFN2
  - 运行模式: inputfile, library

- **ORCA**: DFT/ab initio 方法，精度高但较慢
  - 支持任意 ORCA 输入语法
  - 支持解析频率（AnFreq）
  - 支持并行计算

### 频率计算方法
- **NumFreq**: 数值频率（数值微分）
  - 支持所有 theory
  - 需要多次单点计算
  - 较慢但通用

- **AnFreq**: 解析频率（解析 Hessian）
  - 仅 ORCA 支持
  - 只需一次计算
  - 快速且精确

### 收敛标准
- **ORCA**: ORCA 默认收敛标准
- **GAU**: Gaussian 默认标准
- **GAU_TIGHT**: Gaussian 严格标准
- **GAU_LOOSE**: Gaussian 宽松标准
- **GAU_VERYTIGHT**: Gaussian 非常严格标准

---

## 版本历史

### v2.0 (2026-01-30)
- ✅ 添加双 Theory-Level 支持
- ✅ 添加 AnFreq 解析频率支持
- ✅ 参数模块化（CI, convergence_setting）
- ✅ 代码整理和优化

### v1.0 (2026-01-29)
- ✅ 基础 xTB 工作流实现
- ✅ 9 阶段完整流程
- ✅ 端点匹配验证
