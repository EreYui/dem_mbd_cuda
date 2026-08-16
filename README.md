# DEM–MBD CUDA

[English](README_EN.md) | **中文**

通用、双精度的离散元（DEM）—多体动力学（MBD）耦合程序。正式求解后端仅为 NVIDIA CUDA；无可用 CUDA 设备时程序会清晰报错并退出，不会回退到 CPU 求解。

## 功能

- `CONTROL_MULTIBODY_FLAG = 0`：纯 DEM。
- `CONTROL_MULTIBODY_FLAG = 1`：通用 DEM–MBD；刚体、三角形和颗粒数量均由输入动态确定。
- 每个刚体可通过 `CouplePolyFileName.csv` 的末列 `state` 独立选择受力运动、固定或规定运动。
- GPU 常驻颗粒 SoA、接触历史、排序空间网格、力和力矩；GPU 执行邻域搜索、颗粒/墙面/三角面接触、载荷归约和颗粒 leapfrog 推进。
- CPU 保留配置与文件 I/O、通用小规模 MBD 方程求解，以及每步少量刚体状态/六分量载荷交换。
- 保持既有参数文件、输入文件和输出文件的名称、字段顺序及步编号方式。
- 输入路径中的 `/` 和 `\` 均可使用；`DATA`、`DATA1` 等目录选择方式保持兼容。
- 输出使用 `std::thread` 双缓冲。未启用某类输出时不会为其分配巨大刚体×颗粒数组，也不会下载对应完整快照。

CUDA 架构、数值验证和性能结果见 [迁移与验证报告](docs/CUDA_MIGRATION_VALIDATION.md)。

## 要求

- NVIDIA GPU、兼容驱动和 CUDA Toolkit。
- CMake 3.24 或更高版本。
- C++17/CUDA 17 编译器。
- Windows：Visual Studio 2022 x64 C++ 工具链。
- Ubuntu：Ubuntu 24.04、受当前 CUDA Toolkit 支持的 GCC/G++。

项目默认以 `CMAKE_CUDA_ARCHITECTURES=native` 编译，可在 CMake 配置时显式覆盖。未启用 `--use_fast_math`，颗粒与动力学计算默认使用 `double`。

## Windows 构建与运行

在项目根目录执行：

```bat
build.bat
```

脚本检查 CMake、`nvcc`、驱动/GPU 和 Visual Studio 2022 x64 环境，生成 `build-cuda`，构建 Release，运行 CTest，并执行 CUDA 设备诊断。

```bat
.\build-cuda\Release\dem_mbd.exe --cuda-info
.\build-cuda\Release\dem_mbd.exe DATA
.\build-cuda\Release\dem_mbd.exe DATA1
.\build-cuda\Release\dem_mbd.exe F:\cases\my_case
.\build-cuda\Release\dem_mbd.exe F:\cases\my_case\InputFile\SettingData\world.par
```

## Ubuntu 24.04 构建与运行

```bash
bash build_ubuntu.sh
./build-cuda-linux/dem_mbd DATA
./build-cuda-linux/dem_mbd DATA1
```

也可构建、测试后直接启动算例：

```bash
bash build_ubuntu.sh --run DATA1
```

覆盖 GPU 架构或构建目录示例：

```bash
cmake -S . -B build-cuda-linux -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=120 -DBUILD_TESTING=ON
cmake --build build-cuda-linux --parallel
ctest --test-dir build-cuda-linux --output-on-failure
```

## 算例与输出安全

命令行只给 `DATA1` 时，程序解析为 `Data/DATA1`；也可传入算例目录或 `world.par`。参数文件内旧式 `Data\DATA\...` 路径会重定向到已选择的算例目录。

程序按原规则写入该算例的 `OutputFile`。验证时请复制输入到新的验证目录，避免覆盖已有结果。本仓库的自动测试使用 `validation/cases` 和专用短步目录，不写入正式 `Data/DATA*` 输出。

## 刚体运动状态

`CouplePolyFileName.csv` 在原 24 列之后增加可选的第 25 列 `state`：

```csv
filepath,mass,x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,lambda_x,lambda_y,lambda_z,lambda_w,I1,I2,I3,I4,I5,I6,I7,I8,I9,state
Data\DATA\InputFile\BodySet\plate1.bt,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0
```

- `state=0`：动力刚体，由接触力、重力和其他载荷决定运动。
- `state=1`：固定刚体，质心和姿态保持不变；CSV 中的速度会被置零。接触载荷仍正常计算和输出。
- `state=2`：规定运动刚体，按 CSV 中 `vx,vy,vz` 和 `omega_x,omega_y,omega_z` 作恒定线速度、恒定角速度运动；接触载荷仍计算，但不改变该运动规律。

旧的 24 列文件仍可直接读取，未提供 `state` 时默认按 `state=0` 处理。其他整数或非整数值会报告带绝对路径和刚体编号的输入错误。

## 测试

```bat
ctest --test-dir build-cuda -C Release --output-on-failure
```

测试覆盖几何基准、路径解析、CUDA 设备、自由落体、单墙碰撞、两颗粒接触、单三角面接触、三种刚体运动状态、纯 DEM 短步和动态两刚体 DEM–MBD 短步。数值比较工具：

```bat
.\build-cuda\Release\output_compare.exe reference_output_dir cuda_output_dir
.\build-cuda\Release\state_diagnostics.exe particle_state_file body_state_file
```

## 错误诊断

输入错误会报告文件用途、解析后的绝对路径和具体原因。CUDA kernel 启动均检查错误；调试非法访存可运行：

```bat
compute-sanitizer --tool memcheck --error-exitcode 99 ^
  .\build-cuda\Release\dem_mbd.exe validation\cases\triangle_contact
```

## 文档

- [CUDA 迁移与验证报告](docs/CUDA_MIGRATION_VALIDATION.md)
- [变更记录](CHANGELOG.md)
- [EDEM 对照说明](edem/README.md)
