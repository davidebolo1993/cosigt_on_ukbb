#!/bin/bash
set -euo pipefail

cd /work 2>/dev/null || cd /

# Run process script directly from the image
exec Rscript /usr/local/bin/nebula_run.r "$@"
