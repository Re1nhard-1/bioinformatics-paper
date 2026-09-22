set -euo pipefail
project="$1"
downloads="$project/.tools/linux-runtime"
runtime=/opt/pbmc
mkdir -p "$runtime/bootstrap" "$runtime/mamba/pkgs"
if [ ! -x "$runtime/bootstrap/bin/micromamba" ]; then
    mkdir -p "$runtime/bootstrap/bin"
    cp "$downloads/micromamba" "$runtime/bootstrap/bin/micromamba"
    chmod +x "$runtime/bootstrap/bin/micromamba"
fi
export MAMBA_ROOT_PREFIX="$runtime/mamba"
if [ ! -x "$runtime/salmon/bin/salmon" ]; then
    cp -n "$downloads/salmon-1.10.3-h45fbf2d_5.tar.bz2" "$MAMBA_ROOT_PREFIX/pkgs/"
    "$runtime/bootstrap/bin/micromamba" create --yes --no-rc --prefix "$runtime/salmon" --override-channels --channel conda-forge --channel bioconda --strict-channel-priority 'salmon=1.10.3=h45fbf2d_5'
fi
"$runtime/salmon/bin/salmon" --version
"$runtime/bootstrap/bin/micromamba" list --prefix "$runtime/salmon" --explicit > "$project/environment/salmon/linux-64-explicit.txt"
