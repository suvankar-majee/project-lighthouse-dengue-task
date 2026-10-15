set -euo pipefail
mkdir -p /app/output
python3 /solution/seasonal_r0.py /app/data /app/output
