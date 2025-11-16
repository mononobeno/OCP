#!/bin/bash
set -euo pipefail

# ========================================
# Build and Load GROMACS Prep Image
# ========================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

IMAGE_NAME="gromacs-prep:latest"
CONTEXT_DIR="${ROOT}/images/gromacs-prep"

echo "========================================="
echo "Building GROMACS Prep Image"
echo "========================================="
echo "Image: ${IMAGE_NAME}"
echo "Context: ${CONTEXT_DIR}"
echo ""

cd "${CONTEXT_DIR}"

# Build the image
echo "[1/2] Building Docker image..."
docker build -t "${IMAGE_NAME}" . 2>&1 | tail -20

if [[ $? -eq 0 ]]; then
    echo "  ✓ Build successful"
else
    echo "  ✗ Build failed"
    exit 1
fi

# Verify the image
echo ""
echo "[2/2] Verifying image..."
docker images | grep gromacs-prep || echo "  ✗ Image not found"

# Check if running in kind cluster context
if command -v kind &> /dev/null; then
    echo ""
    echo "[Optional] Loading into kind cluster..."
    kind load docker-image "${IMAGE_NAME}" --name ocp-cluster 2>&1 || {
        echo "  (Skip: kind cluster not found or not running)"
    }
fi

echo ""
echo "========================================="
echo "Build Complete!"
echo "========================================="
echo "Image: ${IMAGE_NAME}"
echo ""
echo "Test with:"
echo "  docker run --rm -v \$(pwd):/workspace -e RUN_ID=7 -e LIGAND_ID=ZINC001241749345_1 ${IMAGE_NAME}"
echo "========================================="
