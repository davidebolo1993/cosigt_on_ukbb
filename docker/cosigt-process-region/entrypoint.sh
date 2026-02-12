#!/bin/bash
set -euo pipefail

REPO_URL="${COSIGT_REPO_URL:-https://github.com/davidebolo1993/cosigt_on_ukbb.git}"
BRANCH="${COSIGT_BRANCH:-main}"

CLONE_DIR="/tmp/cosigt-scripts-$$"

echo "Cloning COSIGT scripts from $REPO_URL (branch: $BRANCH)..."
git clone -b "$BRANCH" --depth 1 "$REPO_URL" "$CLONE_DIR"

cd /work 2>/dev/null || cd /

# Run process.sh with passed arguments
# But override R script paths to use the ones in the image
export CLUSTER_R_SCRIPT="/usr/local/bin/cluster.r"
export OUTLIERS_R_SCRIPT="/usr/local/bin/coverage_outliers.r"

exec bash "$CLONE_DIR/src/process.sh" "$@"

