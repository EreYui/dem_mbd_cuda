#pragma once

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

namespace file_utils {

namespace fs = std::filesystem;

// Configuration files in older cases use Windows backslashes. Convert them
// before constructing std::filesystem::path so the same files also work on
// Linux, where a backslash is an ordinary filename character.
inline fs::path normalizePath(const std::string& rawPath) {
    std::string portablePath = rawPath;
    std::replace(portablePath.begin(), portablePath.end(), '\\', '/');
    return fs::u8path(portablePath).lexically_normal();
}

inline fs::path& caseDirectoryStorage() {
    static fs::path caseDirectory = fs::path("Data") / "DATA";
    return caseDirectory;
}

inline void setCaseDirectory(const fs::path& directory) {
    if (directory.empty()) {
        throw std::runtime_error(
            "Missing case directory / 算例目录不能为空");
    }
    caseDirectoryStorage() = directory.lexically_normal();
}

inline const fs::path& caseDirectory() {
    return caseDirectoryStorage();
}

inline bool hasPathPrefix(const fs::path& path, const fs::path& prefix) {
    auto pathPart = path.begin();
    auto prefixPart = prefix.begin();
    for (; prefixPart != prefix.end(); ++prefixPart, ++pathPart) {
        if (pathPart == path.end() || *pathPart != *prefixPart) {
            return false;
        }
    }
    return true;
}

// Resolve paths against the selected case directory. Existing world.par and
// CSV files may still contain Data/DATA/...; redirect that legacy prefix to
// the case selected on the command line. InputFile/... and OutputFile/...
// paths are also interpreted relative to the selected case directory.
inline fs::path resolveCasePath(const std::string& rawPath) {
    const fs::path path = normalizePath(rawPath);
    if (path.is_absolute()) {
        return path;
    }

    const fs::path legacyCaseDirectory = fs::path("Data") / "DATA";
    if (hasPathPrefix(path, legacyCaseDirectory)) {
        return (caseDirectory() / path.lexically_relative(legacyCaseDirectory))
            .lexically_normal();
    }

    const auto firstPart = path.begin();
    if (firstPart != path.end()
        && (*firstPart == "InputFile" || *firstPart == "OutputFile")) {
        return (caseDirectory() / path).lexically_normal();
    }
    return path;
}

inline std::string absolutePathForMessage(const fs::path& path) {
    std::error_code error;
    const fs::path absolutePath = fs::absolute(path, error);
    return (error ? path : absolutePath).lexically_normal().u8string();
}

inline std::runtime_error pathError(
    const std::string& action,
    const std::string& description,
    const fs::path& path,
    const std::string& detail = {}) {
    std::string message = action + " " + description + ": \""
        + absolutePathForMessage(path) + "\"";
    if (!detail.empty()) {
        message += " (" + detail + ")";
    }
    return std::runtime_error(message);
}

inline fs::path requireInputFile(
    const std::string& rawPath,
    const std::string& description) {
    if (rawPath.empty()) {
        throw std::runtime_error(
            "Missing path / 文件路径为空: " + description);
    }

    const fs::path path = resolveCasePath(rawPath);
    std::error_code error;
    const bool exists = fs::exists(path, error);
    if (error) {
        throw pathError(
            "Cannot inspect / 无法检查", description, path, error.message());
    }
    if (!exists) {
        throw pathError("File not found / 文件不存在", description, path);
    }
    if (!fs::is_regular_file(path, error) || error) {
        const std::string detail = error ? error.message() : "not a regular file";
        throw pathError(
            "Invalid input / 输入路径不是普通文件", description, path, detail);
    }
    return path;
}

inline FILE* openInputFile(
    const std::string& rawPath,
    const std::string& description) {
    const fs::path path = requireInputFile(rawPath, description);
    errno = 0;
    FILE* file = std::fopen(path.string().c_str(), "r");
    if (file == nullptr) {
        const std::string detail = errno == 0
            ? "unknown error"
            : std::strerror(errno);
        throw pathError("Cannot open / 无法打开", description, path, detail);
    }
    return file;
}

inline fs::path prepareOutputFile(
    const std::string& rawPath,
    const std::string& description) {
    if (rawPath.empty()) {
        throw std::runtime_error(
            "Missing output path / 输出文件路径为空: " + description);
    }

    const fs::path path = resolveCasePath(rawPath);
    const fs::path parent = path.parent_path();
    if (!parent.empty()) {
        std::error_code error;
        fs::create_directories(parent, error);
        if (error) {
            throw pathError(
                "Cannot create output directory / 无法创建输出目录",
                description,
                parent,
                error.message());
        }
    }
    return path;
}

}  // namespace file_utils
