# ASH xTB 工作流配置说明

## 快速开始

```bash
# 使用默认配置运行
python ash_xtb_explicit_workflow_smoketest.py

# 使用自定义配置文件
python ash_xtb_explicit_workflow_smoketest.py --config my_config.yaml

# 命令行参数覆盖配置文件
python ash_xtb_explicit_workflow_smoketest.py --neb-images 10 --neb-maxiter 200
```

## 配置文件说明

配置文件 `ash_xtb_workflow_config.yaml` 包含所有可调参数：

- **inputs**: 输入文件路径（reactant, product, ts_guess）
- **system**: 体系属性（charge, mult）
- **xtb**: xTB 理论设置（method, runmode, cores）
- **neb**: NEB 计算参数（images, maxiter, interpolation, runmode, cores）
- **optimization**: 优化参数（maxiter, convergence_setting）
- **tsopt**: TS 优化参数（maxiter, hessian）
- **frequency**: 频率计算参数（npoint, runmode, cores, imag_threshold）
- **irc**: IRC 计算参数（maxiter）
- **output**: 输出设置（base_dir, printlevel）

## 参数优先级

命令行参数 > 配置文件 > 默认值

## 输出

运行成功后会生成：
- `ash_xtb_smoketest_runs/run_YYYYMMDD_HHMMSS/` - 运行目录
- `workflow_summary.json` - 工作流摘要
- `ash_xtb_result.pkl` - 结果数据（可用于后续分析）

## 常用参数调整

```bash
# 快速测试（减少迭代次数）
--neb-maxiter 50 --opt-maxiter 100 --freq-npoint 1

# 高精度计算（增加迭代次数和图像数）
--neb-images 11 --neb-maxiter 300 --opt-maxiter 500

# 使用更多核心
--theory-cores 4
```
