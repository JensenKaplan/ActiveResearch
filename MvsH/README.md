# Magnetization (M vs H) Data
The "M vs H" notebook is used to display all M vs H runs (at different temperatures) collected for a given compound. It displays data in both (emu vs Oe) and (uB vs T). This is used to find the saturation point of our compounds.

## Current Issue
PCF magnetization calculation for single crystal produces results closer to experimentation than that of PCF with powder averaging. I've implemented both, in both LS and J basis to display the effect of implementing powder averaging. Essentially the powder average magnetization is now much greater than the measured, while the single crystal magnetization is right below measurements.