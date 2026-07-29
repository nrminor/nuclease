#!/usr/bin/env bash

if [[ "$(uname -s)" == "Darwin" ]]; then
    if ! command -v xcrun >/dev/null 2>&1; then
        printf '%s\n' \
            "Xcode Command Line Tools are required to build Nuclease on macOS." \
            "Install them with: xcode-select --install" >&2
        return 1
    fi

    export SDKROOT
    SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"

    export CC
    CC="$(xcrun --sdk macosx --find clang)"

    export CXX
    CXX="$(xcrun --sdk macosx --find clang++)"

    export AR
    AR="$(xcrun --sdk macosx --find ar)"

    export CARGO_TARGET_AARCH64_APPLE_DARWIN_LINKER="$CC"
    export CARGO_TARGET_X86_64_APPLE_DARWIN_LINKER="$CC"
fi
