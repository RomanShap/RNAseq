#!/usr/bin/env bash
# run_as_pipeline.sh
#
# Shell wrapper for the AS_pipeline.  Executes the main R script with a
# configuration file.  When run without arguments it defaults to
# config/config.yml under the project root.

set -euo pipefail

# Determine directory of this script and pipeline root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "$SCRIPT_DIR")"

# If a configuration file is provided use it, otherwise fall back to
# config/config.yml relative to the project root.
CONFIG_FILE="${1:-config/config.yml}"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Configuration file not found: $CONFIG_FILE" >&2
  exit 1
fi

echo "Running AS_pipeline_v2 with config: $CONFIG_FILE"

# Invoke the R script.  Use Rscript here to ensure batch execution.
Rscript "$PIPELINE_ROOT/AS_main.R" "$CONFIG_FILE"
