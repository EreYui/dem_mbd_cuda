# EDEM 复现工具

本目录存放可复用于 `DATA*` 算例的 EDEM 2024 转换、运行和打开脚本。算例生成物放在 `Data/<算例名>/edem/`。

## DATA 算例

在项目根目录运行：

```bat
edem\generate_case.bat DATA
edem\run_case.bat DATA -OneStep
edem\open_case.bat DATA
```

- `generate_case.bat` 读取原程序的 `world.par`、`Particles.bt` 和 `CouplePolyFileName.csv`，生成 EDEM 工程。
- `run_case.bat -OneStep` 用原时间步执行 1 步，适合快速检查工程和许可证。
- 去掉 `-OneStep` 后按原算例的 2.5 s 总时长计算；默认使用 CUDA 设备 0，可用 `-Device 1` 选择其他设备。
- `open_case.bat` 在 EDEM 2024 GUI 中打开工程。

脚本使用 EDEM 安装目录自带的 Altair Python 3.8 和 EDEMpy，不会把 EDEMpy 复制进仓库。

## 建模对应关系

- 25,215 个床层颗粒按 `Particles.bt` 的位置、速度、角速度和四元数导入。
- 五个无限平面墙转换成一个开口尺寸为 0.5 m × 0.5 m 的有限 EDEM 容器。
- 0.15 m 立方体转换成一个 EDEM 多面体颗粒，而非“受力几何体”。这样能保留完整 6 自由度、3 kg 质量、输入惯量和初始姿态。
- 接触模型采用 EDEM 内置 `Linear Spring` 与 `Standard Rolling Friction`。
- `MECH_PP_*` 用于床层颗粒之间，`MECH_PT_*` 用于床层与立方体，`MECH_PW_*` 用于颗粒/立方体与容器。

EDEM 内置 Linear Spring 不直接接收固定的 `kN`，而是由材料参数和特征速度推导刚度。转换器在 0.5 m/s 特征速度下匹配床层颗粒的 `MECH_PP_kN`。`mu_T` 和 `beta` 在这组 EDEM 内置模型中没有一一对应参数，因此两套程序应比较下落、接触和沉降趋势，不应预期逐时间步数值完全相同。

EDEM 2024 的多面体颗粒只支持 CUDA 求解器，因此该 EDEM 对照算例需要 EDEM 支持的 NVIDIA GPU；本项目的正式求解同样要求可用的 NVIDIA CUDA 设备。
