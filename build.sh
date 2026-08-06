#!/usr/bin/env bash
# macOS/Linux equivalent of build.bat
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

mkdir -p build
cd build
# -DCMAKE_BUILD_TYPE is needed here (in addition to --config below) because the
# default single-config Makefiles generator on macOS/Linux ignores `--config`
# at build time - only multi-config generators (VS, Xcode, Ninja Multi-Config)
# honor it. Setting both keeps this working regardless of generator.
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
cd ..
