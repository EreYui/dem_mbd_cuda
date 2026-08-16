# ParaView 一键后处理

本目录把原来的两次 MATLAB 转换和多组 VTK 导入简化为一个流程：

```text
ph.*.bt + At.*.bt + CouplePolyFileName.csv 中声明的刚体 .bt 网格
                              |
                              | convert_case.py
                              v
                   particles.pvd + bodies.pvd
                              |
                              v
                          ParaView 时间轴
```

转换程序只使用 Python 3 标准库，不依赖 MATLAB、NumPy 或 ParaView Python 模块。颗粒输出为二进制 `.vtp` 点集；刚体数量从所选算例的 `InputFile/BodySet/CouplePolyFileName.csv` 第一行读取，并在每个时刻合并为一个二进制 `.vtu`。两个 `.pvd` 文件负责时间轴。

## 最简单的使用方式

安装 Python 3 和 ParaView 后，可直接指定算例：

```bat
paraview\open_case.bat DATA3
```

这会转换并打开 `Data/DATA3`。如果经常处理同一个算例，只需编辑 [`case.dat`](case.dat) 中的一行：

```text
CASE_NAME=DATA3
```

之后双击 `open_case.bat` 即可，不需要再修改任何 `.bat` 或 `.py`。`case.dat` 还可设置 `PARAVIEW_EXE`；当前默认值已经是 `C:\Program Files\ParaView 6.0.1\bin\paraview.exe`。

脚本会依次：

1. 自动扫描所选算例 `OutputFile/state_particles/ph.*.bt`；
2. 自动扫描所选算例 `OutputFile/state_bodys/At.*.bt`；
3. 从 `world.par` 读取 `ODE_StepSize`，生成真实物理时间；
4. 把所有颗粒生成一个 `particles.pvd` 时间序列；
5. 读取 `CouplePolyFileName.csv` 中的刚体数量和 `.bt` 网格，把当前算例的全部刚体合并成一个 `bodies.pvd` 时间序列；
6. 启动 ParaView，优先使用 GPU 实例化的 `3D Glyphs` 球体显示颗粒，载入全部刚体并更新时间轴；
7. 保存 `dem_mbd_visualization.pvsm`，以后也可直接打开该状态文件。

如果 ParaView 没有安装在默认目录，先设置：

```bat
set PARAVIEW_EXE=C:\Program Files\ParaView 6.0.1\bin\paraview.exe
paraview\open_case.bat
```

当前机器的 `C:\Program Files\ParaView 6.0.1\bin\paraview.exe` 符合脚本的 `%ProgramFiles%\ParaView*` 自动搜索规则，通常不需要手工设置 `PARAVIEW_EXE`；上面的命令可作为自动搜索失败时的显式覆盖。

命令行算例名优先于 `case.dat`。例如可运行 `paraview\open_case.bat DATA2`，也可附加 `--max-frames 5` 等转换参数。批处理会把项目根目录和最终选中的算例目录传给 ParaView。ParaView 6.0 使用 `--script` 时可能不定义 Python 的 `__file__`；场景脚本不再依赖它，并保留配置文件和当前目录回退。

## 只转换，不启动 ParaView

指定算例运行：

```bat
paraview\prepare_case.bat DATA3
```

或者从项目根目录执行：

```powershell
python paraview/convert_case.py DATA3
```

如果省略 `DATA3`，转换器读取 `case.dat`。高级用法仍可通过 `--case F:\other\case` 指定任意算例路径，但 `open_case.bat` 推荐始终使用简短算例名。

完成后只需在 ParaView 中打开：

```text
Data/DATA3/OutputFile/paraview/particles.pvd
Data/DATA3/OutputFile/paraview/bodies.pvd
```

如果手动使用 `Glyph -> Sphere` 过滤器，设置：

- Scale Array：`radius`；
- Scale Factor：`2`；
- Glyph Mode：`All Points`。

`open_case.py` 已自动完成对应设置，并优先使用 ParaView 的 `3D Glyphs` 显示方式，避免为 30 万个颗粒显式复制球面网格。球源默认半径为 0.5，因此 Scale Factor 必须为 2，球体实际半径才等于颗粒的 `radius`。旧版 ParaView 不支持该显示属性时，脚本才回退到普通 Glyph 过滤器。

## 常用参数

快速预览前 5 帧：

```powershell
python paraview/convert_case.py DATA3 --max-frames 5 --force
```

每隔 4 个已有输出帧转换一次：

```powershell
python paraview/convert_case.py DATA3 --stride 4
```

只处理指定步数范围：

```powershell
python paraview/convert_case.py DATA3 --start-step 20000 --end-step 80000
```

只转换颗粒或刚体：

```powershell
python paraview/convert_case.py DATA3 --no-bodies
python paraview/convert_case.py DATA3 --no-particles
```

强制覆盖已有转换结果：

```powershell
python paraview/convert_case.py DATA3 --force
```

默认会复用比状态文件、刚体 `.bt` 网格和转换脚本更新的 VTK 文件，因此重复运行通常很快。

如果只需要恢复完整 `.pvd` 时间轴、完全不重写已有 VTK 文件：

```powershell
python paraview/convert_case.py DATA3 --index-only
```

## 刚体数量、网格和四元数

转换器以当前算例的 `InputFile/BodySet/CouplePolyFileName.csv` 为准：

- 第一行第一个数值是刚体数量；
- `filepath` 列按顺序给出每个刚体的 `.bt` 网格；
- 即使 CSV 中仍写着旧路径 `Data\DATA\...`，也会从当前所选算例的 `BodySet` 目录读取同名文件；
- `.bt` 网格坐标已经是米，不再依赖从机器人项目复制来的中心体、腿或静态 VTK 模板；
- 每个 `At.*.bt` 的非空行数必须与 CSV 声明的刚体数量一致，否则转换器会明确报错。

程序按本项目的格式直接读取 `At.*.bt` 中的 `[qw,qx,qy,qz]`。只有处理确实采用共轭四元数的旧数据时，才需要显式添加：

```powershell
python paraview/convert_case.py DATA3 --conjugate-body-quaternions
```

## 输出字段

颗粒点数据包含：

- `number`、`status`；
- `mass`、`radius`；
- `speed`、`angular_speed`；
- `velocity`；
- `angular_velocity_body`；
- `quaternion_wxyz`。

刚体网格的点和单元数据都包含 `body_id`，范围为 `0` 到“刚体数量减 1”，顺序与 `CouplePolyFileName.csv` 的数据行一致。可在 ParaView 中用它着色或选择单个刚体。

## 注意事项

- 转换结果位于所选算例的 `OutputFile/paraview/`，例如 `Data/DATA3/OutputFile/paraview/`；它属于可重新生成的数据，不应提交到 Git。
- `--max-frames` 会让 `.pvd` 只引用选中的帧；要恢复全部时间步，重新执行不带该参数的命令即可。
- 如果只启用了颗粒或刚体状态输出，转换器会生成存在的那一类时间序列。
- 颗粒很多时，Glyph 会占用较多显存；预览可降低球体分辨率或用 `--stride` 减少时间帧，但不要把 Glyph Mode 改为随机抽样后再用于正式渲染。
