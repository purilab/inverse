

import pandas as pd
!pip install pymzml
import pymzml as pz

#%%mass convert of mass files
import os
os.chdir("insert the path to the metatlas package here") #change file path as needed
from metatlas.io import mzml_loader as ml


input_folder = "XXX/convert"  #create a folder called 'convert' and import all .mzML files to this folder that should be converted. Change file path as necessary. 
output_folder = "XXX/convert"  #change path as needed. Should be the same as in line above. 

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for file_name in os.listdir(input_folder):
    if file_name.endswith(".mzML"):
        input_file_path = os.path.join(input_folder, file_name)
        output_file_name = os.path.splitext(file_name)[0] + ".h5"
        output_file_path = os.path.join(output_folder, output_file_name)

        ml.mzml_to_hdf(input_file_path, out_file=output_file_path, debug=False)

