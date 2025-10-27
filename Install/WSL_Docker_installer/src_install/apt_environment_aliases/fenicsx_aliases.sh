#!/bin/bash
set -e

ALIAS_FILE=~/.bash_aliases
# Use environment variable FULL_IMAGE passed from the batch script
FULL_IMAGE=${FULL_IMAGE:-fenicsx:v0.9.0}

# Define the aliases. The volume mount (-v) ensures the current working directory in WSL is accessible.
ALIAS_FENICSX='alias fenicsx="docker run -it --rm -v \$(pwd):/home/fenicsx/shared -w /home/fenicsx/shared '$FULL_IMAGE' bash"'
ALIAS_FENICSX_RUN='alias fenicsx-run="docker run --rm -v \$(pwd):/home/fenicsx/shared -w /home/fenicsx/shared '$FULL_IMAGE' bash -c"'

# Ensure the alias file exists
touch "$ALIAS_FILE"

## Replacement Logic: Remove existing aliases first ##
# Use sed to remove lines that start with the alias names.
# -i (or -i.bak on some systems) edits the file in place.
# The 'f' command is used to suppress the printing of lines after a successful substitution/deletion.
sed -i.bak '/^alias fenicsx=/d' "$ALIAS_FILE"
sed -i.bak '/^alias fenicsx-run=/d' "$ALIAS_FILE"
# Remove the backup files created by sed
rm -f "$ALIAS_FILE.bak"

## Add the new aliases ##
echo "$ALIAS_FENICSX" >> "$ALIAS_FILE"
echo "$ALIAS_FENICSX_RUN" >> "$ALIAS_FILE"

# Ensure .bashrc sources .bash_aliases (only if it doesn't already)
grep -qxF 'source ~/.bash_aliases' ~/.bashrc || echo 'source ~/.bash_aliases' >> ~/.bashrc

echo "FEniCSx aliases updated (replaced) in ~/.bash_aliases using image: $FULL_IMAGE."
echo "You may need to run 'source ~/.bash_aliases' or open a new terminal for changes to take effect."