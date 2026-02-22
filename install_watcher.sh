#!/bin/bash
set -e

GREEN='\033[0;32m'
CYAN='\033[0;36m'
RED='\033[0;31m'
NC='\033[0m'

REPO=~/scripts/informational-asymmetry
VENV=$REPO/.venv

echo ""
echo "═══════════════════════════════════════════"
echo "  Informational Asymmetry — Watcher Setup"
echo "═══════════════════════════════════════════"
echo ""

# 1. Create venv if it doesn't exist
if [ ! -d "$VENV" ]; then
  echo -e "${CYAN}Creating virtualenv...${NC}"
  python3 -m venv "$VENV"
fi

# 2. Install watchdog
echo -e "${CYAN}Installing watchdog...${NC}"
"$VENV/bin/pip" install --quiet watchdog

# 3. Copy the watcher script into the repo
echo -e "${CYAN}Installing ia_watcher.py...${NC}"
cp ia_watcher.py "$REPO/ia_watcher.py"
chmod +x "$REPO/ia_watcher.py"

# 4. Install systemd user service
echo -e "${CYAN}Installing systemd user service...${NC}"
mkdir -p ~/.config/systemd/user
cp ia-watcher.service ~/.config/systemd/user/ia-watcher.service

# Update the service ExecStart to use the correct venv python
sed -i "s|ExecStart=.*|ExecStart=$VENV/bin/python $REPO/ia_watcher.py|" \
  ~/.config/systemd/user/ia-watcher.service

systemctl --user daemon-reload
systemctl --user enable ia-watcher.service
systemctl --user start ia-watcher.service

echo ""
echo -e "${GREEN}✓ Watcher installed and running!${NC}"
echo ""
echo "  Useful commands:"
echo "    systemctl --user status ia-watcher    # check status"
echo "    journalctl --user -u ia-watcher -f    # live logs"
echo "    systemctl --user restart ia-watcher   # restart"
echo "    systemctl --user stop ia-watcher      # stop"
echo ""
echo "  The watcher is monitoring ~/Downloads"
echo "  Files are routed to ~/scripts/informational-asymmetry/"
echo "  and auto-committed + pushed to GitHub."
echo ""

# 5. Quick test — show status
systemctl --user status ia-watcher --no-pager || true