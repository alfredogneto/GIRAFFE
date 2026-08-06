#!/usr/bin/env bash
# macOS/Linux equivalent of install_dependencies.bat
#
# Prerequisites: Homebrew (macOS) - https://brew.sh
#
# Installs:
#   - exprtk    (header-only, cloned from GitHub, same as Windows)
#   - eigen     (header-only, cloned from GitLab, same as Windows - currently unused by Giraffe itself)
#   - arpack-ng (via Homebrew, pulls in OpenBLAS for BLAS/LAPACK)
#   - libomp    (via Homebrew, OpenMP runtime for Clang)
#
# Also sets GIRAFFE_INSTALL in the user's shell profile, mirroring what
# install_dependencies.bat does for the Windows user environment variable.

set -euo pipefail

GIRAFFE_INSTALL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$GIRAFFE_INSTALL"

echo "Giraffe install directory: $GIRAFFE_INSTALL"

if [[ "$(uname)" == "Darwin" ]]; then
	if ! command -v brew >/dev/null 2>&1; then
		echo "Homebrew not found. Install it from https://brew.sh and re-run this script." >&2
		exit 1
	fi
	echo
	echo "Installing arpack, libomp, eigen, pkg-config via Homebrew..."
	brew install arpack libomp eigen pkg-config
else
	echo "This script targets macOS. On Linux, install libopenblas-dev, libarpack2-dev," >&2
	echo "libomp-dev and pkg-config with your distro's package manager." >&2
fi

echo
echo "Searching for exprtk..."
if [[ ! -d dependencies/exprtk ]]; then
	echo "exprtk not found, cloning exprtk..."
	git clone https://github.com/ArashPartow/exprtk.git dependencies/exprtk
else
	echo "exprtk found."
fi

echo
echo "Searching for eigen..."
if [[ ! -d dependencies/eigen ]]; then
	echo "eigen not found, cloning eigen..."
	git clone https://gitlab.com/libeigen/eigen.git dependencies/eigen
else
	echo "eigen found."
fi

echo
echo "Setting GIRAFFE_INSTALL in shell profile..."
PROFILE_FILE="$HOME/.zprofile"
[[ "$SHELL" == *bash* ]] && PROFILE_FILE="$HOME/.bash_profile"

if grep -q "^export GIRAFFE_INSTALL=" "$PROFILE_FILE" 2>/dev/null; then
	sed -i '' -E "s#^export GIRAFFE_INSTALL=.*#export GIRAFFE_INSTALL=\"$GIRAFFE_INSTALL\"#" "$PROFILE_FILE"
else
	{
		echo ""
		echo "# Added by Giraffe install_dependencies.sh"
		echo "export GIRAFFE_INSTALL=\"$GIRAFFE_INSTALL\""
	} >> "$PROFILE_FILE"
fi
export GIRAFFE_INSTALL

echo
echo "Step 1 of the setup is ready."
echo "Restart your terminal (or 'source $PROFILE_FILE') so GIRAFFE_INSTALL is available."
echo "After this, proceed with the CMake build by running ./build.sh"
