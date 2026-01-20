import pySRProxy
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")
import uproot
import pandas as pd
import numpy as np

NDLArXLo = -346.9
NDLArXHi = 346.9
NDLArYLo = -215.5
NDLArYHi = 81.7
NDLArZLo = 418.2
NDLArZHi = 913.3

def isFV(vtx):
  return (NDLArXLo < vtx.x < NDLArXHi and
          NDLArYLo < vtx.y < NDLArYHi and
          NDLArZLo < vtx.z < NDLArZHi)

def allPrimariesContained(inx):
  for part in inx.part.dlp:
    if part.primary and not part.contained:
      return False
  return True

def longestTrackIsMuon(inx):
  longest_length = 0
  longest_length_pid = -999
  for part in inx.part.dlp:
    length = ( (part.end.x - part.start.x)**2 + 
                (part.end.y - part.start.y)**2 + 
                (part.end.z - part.start.z)**2 )**0.5
    if length > longest_length:
      longest_length = length
      longest_length_pid = part.pdg

  return longest_length_pid == 13

files_dir = "files/"

files = os.listdir(files_dir)

events = {
  "Elep": [],
  "LongestTrackLength": [],
  "TrackLengths": []
}

for f in tqdm(files):
  rdr = pySRProxy.StandardRecordFileReader("cafTree", ["%s/%s" % (files_dir, f)])

  for ent in tqdm(rdr, leave=False):
    for i, inx in tqdm(enumerate(ent.common.ixn.dlp), leave=False, total=ent.common.ixn.ndlp):
      if not isFV(inx.vtx):
        continue
      if not allPrimariesContained(inx):
        continue
      if not longestTrackIsMuon(inx):
        continue
      
      track_lengths = [track.len_cm for track in ent.nd.lar.dlp[i].tracks]
      events["TrackLengths"].append(track_lengths)
      events["LongestTrackLength"].append(max(track_lengths))
      
      track_energies = [track.E for track in ent.nd.lar.dlp[i].tracks]
      events["Elep"].append(sum(track_energies))

df = pd.DataFrame(events)

# bunch of additional with dummy values columns for easy integration with MaCh3 Tutorial
df["Q2"] = 1.
df["Enu_true"] = df["Elep"]
df["Mode"] = 1
df["tgt"] = 1000060120
df["PDGLep"] = 14
df["CosLep"] = 1.0
df["flagCC0pi"] = 1
df["flagCC1pip"] = 0
df["flagCC1pim"] = 0

# file output for MaCh3
with uproot.recreate("selected_events.root") as f:
  f["FlatTree_VARS"] = df

# quick plots with Matplotlib
plt.hist(events["LongestTrackLength"], bins=50)
plt.xlabel("Longest Track Length (cm)")
plt.ylabel("Counts")
plt.savefig("longest_track_lengths.pdf")
plt.savefig("longest_track_lengths.png")
plt.clf()

plt.hist(events["Elep"], bins=50)
plt.xlabel("Energy (MeV)")
plt.ylabel("Counts")
plt.savefig("energy.pdf")
plt.savefig("energy.png")
plt.clf()