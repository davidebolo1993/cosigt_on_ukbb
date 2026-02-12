#!/bin/bash
set -e

# Pass all arguments to organize.sh
exec /usr/local/bin/organize.sh "$@"
