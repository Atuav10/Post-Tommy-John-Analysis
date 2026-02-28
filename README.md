# Post-Tommy-John-Analysis
Using quasi-experimental methods to analyze pitcher performance post-tommy john surgery
- Dartmouth College Presidential Scholar project
- Dartmouth College Independent Study project
- SABR Analytics poster presentation
- Collaboration with Professor Michael Herron (Program in Quantitative Social Science)

# Overview

## Inputs
- Player statistics from baseballR
- List of qualified pitchers who underwent surgery - Pulled from public Tommy John datasets (@MLBPlayerAnalys on X)
- FIP constants for each year

## Causal methodologies
- Synthetic Control Method (SCM)
- Difference-in-Differences (DiD)

## Pipeline
1) Create a ccompiled dataset of every treeated (underwent surgery) player
2) Obtain pre-treatment and post-treatment periods for the aforementioned players
3) Create a customized donor pool for each treated player
4) Apply both causal methodologies for each player
5) Calculate the Average Treatment Effect
