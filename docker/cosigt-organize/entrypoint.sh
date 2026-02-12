#!/bin/bash
set -euo pipefail

# Clone to /tmp which is always writable in Singularity
REPO_URL="${COSIGT_REPO_URL:-https://github.com/davidebolo1993/cosigt_on_ukbb.git}"
BRANCH="${COSIGT_BRANCH:-main}"

CLONE_DIR="/tmp/cosigt-scripts-$$"  # Use PID for uniqueness

echo "Cloning COSIGT scripts from $REPO_URL (branch: $BRANCH)..."
git clone -b "$BRANCH" --depth 1 "$REPO_URL" "$CLONE_DIR"

# Change to work directory (for relative paths)
cd /work 2>/dev/null || cd /

# Run organize.py with passed arguments
exec python "$CLONE_DIR/src/organize.py" "$@"

