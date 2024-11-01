# Investigating the impact of antagonistic selection in vector-transmitted parasitic diseases: Python code
### *Author: Aiden Lewis*
*Written for Python 3.11.5 and the latest version of Windows 11 as of October 2024.*

### To run the code
The code can be run by typing `python main.py` into a command line terminal within this directory, presuming that Python is already installed. It must be run in an activated conda environment, and the required packages are as follows:
- `matplotlib`
- `numpy`
- `pandas`
- `os`
- `math`
- `random`
- `scipy`
- `PIL`
- `io`
- `colorist`
- `typing`
- `pickle`
- `time`
- `cProfile`
- `pstats`

If any of those last three packages -- `time`, `cProfile`, and `pstats` -- are unavailable for whatever reason, they can be removed relatively easily if necessary. Removing the others would be variously difficult to outright impossible. Altering the imports is, however, very much not recommended.

### Output and meanings
When the code finishes running, there are three main sorts of files it will produce: host data (denoted as `h1`), vector data (denoted as `vec`), and `net_output.opt`. For the first two, the main image files are their genotype frequency plots over time. The `.csv` files record all relevant data explicitly. The image files with names containing `freq` are based on the same data as the main image files; they simply depict the results as lines instead of shaded areas. As such graphs are harder to read, they do not generally see use. `net_output.opt` is a plain text file that can be opened in any text editor, e.g. Notepad -- the unusual extension `.opt` is only there to keep it from being caught by discarded file deletion methods (`.txt` files are used for some purposes mid-code, as are `.dat` files, and so are cleared out regularly). If multiple simulations are being run, the results will be stored in the `full_outputs` directory.

The `hists` folder is where histogram plots are stored, if the code is told to produce them. They display the infected individuals' internal allele frequency distributions at various points in time. The `full_outputs` folder is where the results of multi-simulation executions are stored (under a directory as specified in the main file). All the data for each simulation is stored in its own folder, and high-level data (i.e. data integrating the results of multiple simulations) is stored in the folder itself.

If `pyinstaller` is available, `compile.ps1` can be run in a PowerShell terminal to compile the Python code into an executable file. This may offer some (minor) increases to speed, but it would also be fairly inflexible with respect to adjusting the simulation's parameters and initial conditions.

### File outline
- `allele.py`: The `allele` class. Primarily used to relate a character representing the allele to its selection and transmission biases.
- `color.py`: A library of basic functions that have to do with the colors used during the graph generation process.
- `data_lib.py`: A library of functions oriented around data analysis and statistics.
- `func_lib.py`: A library of general-purpose functions that are mostly, though not entirely, semantically independent of this particular project.
- `gen_funcs.py`: A library of functions semantically oriented around genetics.
- `individual.py`: The `individual` class representing an infected host or vector and managing its internal parasite behavior.
- `main.py`: The main file of the code. This is the only file one needs to look at or modify in order to run different simulations themselves.
- `model.py`: The `SIR` class. Manages the behavior of a particular strain (i.e. genotype) in a particular population. With the standard setup, that means **D**, **d** in the host and **DD**, **Dd**, **dd** in the vector.
- `population.py`: The `population` class. Handles host and vector populations as a whole, not just their strains of infection.
- `sim_lib.py`: A library containing the shell of the simulation. It translates the given parameters into the format required of the simulation, then runs it over time.

Refer to each file for further details on its functions and purposes. The main area for user inputs is denoted accordingly in `main.py`.