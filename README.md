# Left Ventricle flattening
Author: Marta Nuñez-Garcia (marnugar@gmail.com)

## About
Quasi conformal flattening with regional constraints. This algorithm is an extension (simplification) of the left atrial flattening method described in:[*Fast quasi-conformal regional flattening of the left atrium*. Marta Nuñez-Garcia, Gabriel Bernardino, Francisco Alarcón, Gala Caixal, Lluís Mont, Oscar Camara, and Constantine Butakoff.  IEEE Transactions on Visualization and Computer Graphics (2020)](https://ieeexplore.ieee.org/abstract/document/8959311). Please cite this reference when using this code. Preprint available at: [arXiv:1811.06896.](https://arxiv.org/pdf/1811.06896.pdf) The code runs in Linux and Windows. 

Example:

![Example image](https://github.com/martanunez/LV_flattening/blob/master/example.png)

## Code
[Python](https://www.python.org/) scripts depending (basically) on [VTK](https://vtk.org/) and [VMTK](http://www.vmtk.org/). 


## Instructions
Clone the repository:
```
git clone https://github.com/martanunez/LV_flattening

cd LV_flattening
```

## Usage
```
flat_LV.py [-h] [--lv_meshfile PATH] [--rv_meshfile PATH] [--n N]

Arguments:
  --lv_meshfile PATH  path to LV input mesh
  --rv_meshfile PATH  path to RV input mesh

optional arguments:
  -h, --help          show this help message and exit
  --n N               n for radial displacement
```


## Usage example
```
python flat_LV.py --lv_meshfile data/lv_mesh.vtk --rv_meshfile data/rv_mesh.vtk 

```

## Dependencies
The scripts in this repository were successfully run with:
1. Ubuntu 16.04
    - [Python](https://www.python.org/) 2.7.12
    - [VMTK](http://www.vmtk.org/) 1.4
    - [VTK](https://vtk.org/) 8.1.0
2. Windows 8.1
    - [Python](https://www.python.org/) 3.6.4
    - [VMTK](http://www.vmtk.org/) 1.4
    - [VTK](https://vtk.org/) 8.1.0
  

### Python packages installation
To install VMTK follow the instructions [here](http://www.vmtk.org/download/). The easiest way is installing the VMTK [conda](https://docs.conda.io/en/latest/) package (it additionally includes VTK, NumPy, etc.). It is recommended to create an environment where VMTK is going to be installed and activate it:

```
conda create --name vmtk_env
conda activate vmtk_env
```
Then, install vmtk:
```
conda install -c vmtk vtk itk vmtk
```
<!--Activate the environment when needed using:
```
source activate vmtk_env
```-->
You can also build VMTK from source if you wish, for example, to use a specific VTK version. Instructions can be found [here.](http://www.vmtk.org/download/)

## Important note
You may need to slightly modify vmtkcenterlines.py from the VMTK package if you encounter the following error when running 1_mesh_standardisation.py:

```
     for i in range(len(self.SourcePoints)/3):
TypeError: 'float' object cannot be interpreted as an integer
```

Find vmtkcenterlines.py file and edit as follows:

Line 128: change ***for i in range(len(self.SourcePoints)/3):*** by ***for i in range(len(self.SourcePoints)//3):***

Line 133: change ***for i in range(len(self.TargetPoints)/3):*** by ***for i in range(len(self.TargetPoints)//3):*** 


## License
The code in this repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details: [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)
