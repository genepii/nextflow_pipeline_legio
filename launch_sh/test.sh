#!/bin/bash

echo "pwd : "
echo $(pwd)

echo "test :"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "$SCRIPT_DIR"

ls "${SCRIPT_DIR}/../config/qiime2.config"