#!/bin/bash
set -e
cd ~/scripts/informational-asymmetry
GREEN='\033[0;32m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'
echo ""
echo "═══════════════════════════════════════════"
echo "  Informational Asymmetry — Production Push"
echo "═══════════════════════════════════════════"
if git diff --quiet && git diff --cached --quiet && [ -z "$(git ls-files --others --exclude-standard)" ]; then
  echo -e "\n  No changes to push.\n"
  exit 0
fi
echo ""
git status --short
if [ -n "$1" ]; then
  MSG="$1"
else
  read -p "  Commit message: " MSG
fi
echo -e "\n${CYAN}  Running verification suite...${NC}"
cd verification_suite
if python3 run_all.py 2>&1 | grep -q "ALL TESTS PASSED"; then
  echo -e "  ${GREEN}✓ All tests passing${NC}"
else
  echo -e "  ${RED}✗ TESTS FAILED — aborting${NC}"
  exit 1
fi
cd ~/scripts/informational-asymmetry
git add -A
git commit -m "$MSG"
git push
echo -e "\n${GREEN}  ✓ Pushed. Site updates in ~90 seconds.${NC}\n"
