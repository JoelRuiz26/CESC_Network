#!/bin/bash

echo "Enter the path to the file you want to force add (e.g., path/to/file.rds):"
read file

if [ ! -f "$file" ]; then
  echo "Error: '$file' does not exist."
  exit 1
fi

echo "Enter a commit message (or press Enter for default):"
read msg

if [ -z "$msg" ]; then
  msg="Force add important file: $file"
fi

git add -f "$file"
git commit -m "$msg"
git push

echo "Done! '$file' has been force-added and pushed to the repository."
