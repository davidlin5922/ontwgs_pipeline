#!/usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt

readLength = pd.read_csv(sys.argv[1], sep="\t")
mapQuality = pd.read_csv(sys.argv[2], sep="\t")
depth = pd.read_csv(sys.argv[3], sep="\t")
sample_id = sys.argv[4]

plt.figure(figsize=(10,5))
plt.bar(mapQuality.iloc[:, 0], mapQuality.iloc[:, 1], width=1.0)

plt.xlabel("Mapping Quality")
plt.ylabel("Count")
plt.title("Mapping Quality Histogram")

plt.tight_layout()
plt.savefig(f"./{sample_id}/mapQuality.png", dpi=150)
plt.close()



plt.figure(figsize=(10,5))
plt.bar(depth.iloc[:, 0], depth.iloc[:, 2], width=1.0)

plt.xlabel("Sequencing Depth")
plt.xticks(rotation=45)
plt.ylabel("Frequency")
plt.title("Alignment Depth Distribution")

plt.tight_layout()
plt.savefig(f"./{sample_id}/depth.png", dpi=150)
plt.close()



plt.figure(figsize=(10,5))
plt.bar(mapQuality.iloc[:, 0], mapQuality.iloc[:, 1], width=1.0)

plt.xlabel("Read Length")
plt.ylabel("Count")
plt.title("Read Length Distribution")

plt.tight_layout()
plt.savefig(f"./{sample_id}/readLength.png", dpi=150)
plt.close()
