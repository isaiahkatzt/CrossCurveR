#!/usr/bin/env bash

# File to update
GITIGNORE_FILE=".gitignore"

# Entries to add
ENTRIES=(
  ".DS_Store"
  "reconstructed-data/"
  "*.Rproj"
)

# Create .gitignore if it doesn't exist
if [ ! -f "$GITIGNORE_FILE" ]; then
  touch "$GITIGNORE_FILE"
fi

# Append entries if not already present
for entry in "${ENTRIES[@]}"; do
  if ! grep -qxF "$entry" "$GITIGNORE_FILE"; then
    echo "$entry" >> "$GITIGNORE_FILE"
  fi
done

echo ".gitignore updated."
