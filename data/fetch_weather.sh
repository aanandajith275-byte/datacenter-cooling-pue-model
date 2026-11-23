#!/usr/bin/env bash
# Downloads the three EnergyPlus (EPW) typical-year weather files used by the study into data/weather/.
# Source: the Sinergym project (github.com/ugr-sail/sinergym), which redistributes IWEC / TMY3 files.
# IWEC data are (c) ASHRAE, so they are fetched rather than committed here.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p weather
tmp=$(mktemp -d)
git clone -q --depth 1 --filter=blob:none --no-checkout https://github.com/ugr-sail/sinergym "$tmp/sg"
cd "$tmp/sg"
for f in JPN_Tokyo.Hyakuri.477150_IWEC.epw USA_AZ_Davis-Monthan.AFB.722745_TMY3.epw SWE_Stockholm.Arlanda.024600_IWEC.epw; do
  git checkout -q HEAD -- "sinergym/data/weather/$f"
  cp "sinergym/data/weather/$f" "$OLDPWD/weather/"
done
rm -rf "$tmp"
echo "weather files in $(cd "$OLDPWD/weather" && pwd)"
