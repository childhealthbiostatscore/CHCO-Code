# Lambda System Event Log
IP: 10.158.40.243

> Maintains a history of boots, reboots, crashes, and suspends.

---

## Template Entry
- **Date/Time (local):** YYYY-MM-DD HH:MM TZ
- **Event:** Boot | Reboot | Suspend/Resume | Crash | Shutdown
- **Detection:** (command you used, e.g. `journalctl -b -1`)
- **Symptoms:** (what you observed, e.g. SSH dropped, fans loud, screen black)
- **Logs:** (paste relevant lines from `journalctl`, `dmesg`, etc.)
- **Notes/Resolution:** (what you did, what you plan to change)

---

### 2025-09-30 10:30 PDT
- **Event:** Resume from suspend  
- **Detection:** `journalctl -u NetworkManager -b`  
- **Symptoms:** SSH disconnected until physical wake. Looked like NetworkManager was asleep based on log.
- **Logs:**  YC added this line `sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target` to mask suspend/hibernate targets in the Network Manager. To re-enable later, run `sudo systemctl unmask sleep.target suspend.target hibernate.target hybrid-sleep.target`.
- **Notes/Resolution:** HH noticed the ssh and remote desktop won't connect and noted in data-science Slack channel. YC physically woke the lambda and everything was working and connected again. YC investigated ways to install caffeine, but the app was unavailable for linux systems. "Sleep" was masked in NetworkManager instead in efforts to resolve. YC also created this .md log and set up automated logs for system booting events, suspend events, and resume events. These logs can be accessed `cat /var/log/custom_boots.log`. Team plans on rebooting lambda every Friday for maintenance.