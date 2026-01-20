# Using `pySRProxy` to perform simples cuts and produce a summary file

This example was generated as part of pre-workshop task for a joint LBL + ND-Sim-Reco meeting. The intention was to perform some selection on events from a CAF file, and then feed that into your 'favourite' analysis software, which in this case was MaCh3. That was achieved by outputting a ROOT file with a TTree in the same format as used in the MaCh3 [tutorial](https://github.com/mach3-software/MaCh3Tutorial).

## The selections

- Reconstructed vertex is in NDLArFV
- All primary reconstructed particles are fully contained
- The longest reconstructed particle track was PID'd as a muon (PDG==13)

NDLArFV is 25cm 'in' from each face of the following detector bounding box:

- NDLArXLo = -346.9;
- NDLArXHi = 346.9;
- NDLArYLo = -215.5;
- NDLArYHi = 81.7;
- NDLArZLo = 418.2;
- NDLArZHi = 913.3;

These selections are also specified [here](../../duneanasel/nd/ndlar/Selections.h) as part of this `duneanasel` repository/package.

## Running it

First step is to make sure you can use the `pySRProxy` module. At the time of writing, this is done by cloning a particular branch of `duneanaobj`, building it with the correct `cmake` options, and then exporting the relevant `LD_LIBRARY_PATH` and `PYTHONPATH`. 

```
git clone --branch feature/Clingification_rebase2024 https://github.com/luketpickering/duneanaobj.git
cd duneanaobj; mkdir build; cd build;
cmake .. -Dduneanaobj_PYTHON_ENABLED=ON -DCMAKE_INSTALL_PREFIX=$(pwd)/Linux
make install
export LD_LIBRARY_PATH=$(pwd)/Linux/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=$(pwd)/Linux/python/3.9:${PYTHONPATH}
```

If you struggle to find an environment that works well for this, you can use a container like:
```
docker run -it --rm \         
-v "$(pwd)":/work \
picker24/buildbox:alma9 \
bash
```
and then `cd /work` and then run the instructions above to clone and build `duneanaobj`. 

To run the `cuts_with_summary.py` file, you'll also need to install some python dependencies, those being
- tqdm
- matplotlib
- mplhep
- uproot
- pandas
- numpy
which is most easily done with a python virtual env.

```
python -m venv env
source env/bin/activate
pip install tqdm matplotlib mplhep uproot pandas numpy
```

You'll also need some CAF files to run over that exist in `files/`. Then do `python cuts_with_summary.py`.