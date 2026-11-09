# CK8011D 轮对车削动力学模型

用 MATLAB R2025b 打开本目录后，运行 `run_ck8011d_case` 可复现刀具11、采集3的单工况仿真。采集文件、截取时段、转速、刀台移速、切深和切削力系数均集中在 `ck8011d_case_config.m`。原始 WAV 通过配置中的本地路径读取。

`wheel_lathe_lumped_demo` 提供独立的 33 自由度结构求解入口。单工况默认输出到 `matlab_result_tool11_acq03_calibrated`，使用前一段实测 PSD 拟合等效动态切削力，用后一段比较带内 RMS 和主频。刀杆后端径向、轴向、切向依次使用 CH2、CH1、CH3。CH4 为声音，保留 FS。图件导出为中英文 SVG，基线仿真位于 `matlab_result_tool11_acq03`。

工况来源、假设及复现步骤见 [MATLAB 使用说明](MATLAB_使用说明.md)。
