#!/bin/bash

echo "=== Current Running Jobs ==="
squeue -u $USER --format="%.18i %.30j %.8T %.10M %.6C %.7m"

echo ""
echo "=== Completed Jobs Summary ==="
sacct -S $(date -d '1 hour ago' +%Y-%m-%d-%H:%M) \
      --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS \
      | grep nebula

echo ""
echo "=== Failed Jobs ==="
sacct -S $(date -d '1 hour ago' +%Y-%m-%d-%H:%M) \
      --format=JobID,JobName%30,State,ExitCode \
      | grep -E "FAILED|CANCELLED|TIMEOUT"

echo ""
echo "=== Output Files Generated ==="
ls -lh logs/output/*.out 2>/dev/null | tail -5

echo ""
echo "=== Check for Errors ==="
grep -l "Error" logs/error/*.err 2>/dev/null