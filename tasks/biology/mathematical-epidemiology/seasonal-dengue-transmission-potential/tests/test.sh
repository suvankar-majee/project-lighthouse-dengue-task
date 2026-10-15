mkdir -p /logs/verifier
cd /app
python3 -m pytest /tests/test_outputs.py -rA -p no:cacheprovider
if [ $? -eq 0 ]; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
