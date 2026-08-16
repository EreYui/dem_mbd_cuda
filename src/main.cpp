// CoupDyn.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// Created by zkm on 2024.12.13
// Modified from szj and yyu's program

#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <string>
#include "multibody.h"
#include "particle.h"
#include "force.h"
#include "constant.h"
#include"simulation.h"
#include "file_utils.h"
#include "cuda_backend.h"

namespace {

struct LaunchPaths {
    std::filesystem::path caseDirectory;
    std::filesystem::path parameterFile;
};

bool looksLikeParameterFile(const std::filesystem::path& path) {
    if (path.extension() == ".par") {
        return true;
    }
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error;
}

std::filesystem::path inferCaseDirectory(
    const std::filesystem::path& parameterFile)
{
    const auto settingDirectory = parameterFile.parent_path();
    const auto inputDirectory = settingDirectory.parent_path();
    if (settingDirectory.filename() == "SettingData"
        && inputDirectory.filename() == "InputFile") {
        return inputDirectory.parent_path();
    }
    return std::filesystem::path("Data") / "DATA";
}

LaunchPaths resolveLaunchPaths(int argc, char* argv[]) {
    if (argc == 1) {
        const auto caseDirectory = std::filesystem::path("Data") / "DATA";
        return {
            caseDirectory,
            caseDirectory / "InputFile" / "SettingData" / "world.par"
        };
    }

    const auto argument = file_utils::normalizePath(argv[1]);
    if (looksLikeParameterFile(argument)) {
        return {inferCaseDirectory(argument), argument};
    }

    const auto caseDirectory = argument.is_absolute() || argument.has_parent_path()
        ? argument
        : std::filesystem::path("Data") / argument;
    return {
        caseDirectory,
        caseDirectory / "InputFile" / "SettingData" / "world.par"
    };
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc > 2) {
        std::cerr << "Usage / 用法: " << argv[0]
                  << " [CASE_NAME|CASE_DIRECTORY|path/to/world.par]"
                  << std::endl;
        return EXIT_FAILURE;
    }

    try {
        if (argc == 2 && std::string(argv[1]) == "--cuda-info") {
            if (!CudaBackendAvailable()) {
                throw std::runtime_error(
                    "No usable NVIDIA CUDA device was found / 未检测到可用的 NVIDIA CUDA 设备");
            }
            std::cout << "CUDA device / CUDA 设备: "
                      << CudaDeviceDescription() << std::endl;
            return EXIT_SUCCESS;
        }
        const LaunchPaths launchPaths = resolveLaunchPaths(argc, argv);
        file_utils::setCaseDirectory(launchPaths.caseDirectory);

        std::cout << "Case directory / 算例目录: "
                  << file_utils::absolutePathForMessage(launchPaths.caseDirectory)
                  << '\n'
                  << "Parameter file / 参数文件: "
                  << file_utils::absolutePathForMessage(launchPaths.parameterFile)
                  << std::endl;

        Simulation sim;
        sim.init(launchPaths.parameterFile.u8string());
        sim.run();
    }
    catch (const std::exception& error) {
        std::error_code pathError;
        const auto workingDirectory = std::filesystem::current_path(pathError);
        std::cerr << "\n[ERROR] Simulation stopped / 仿真已终止\n"
                  << "Reason / 原因: " << error.what() << '\n';
        if (!pathError) {
            std::cerr << "Working directory / 当前工作目录: "
                      << workingDirectory.u8string() << '\n';
        }
        std::cerr << "Please check the configured path and file permissions."
                  << " / 请检查配置路径及文件权限。" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
