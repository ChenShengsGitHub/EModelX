# EModelX
EModelX is a method for automatic cryo-EM protein complex structure modeling.
![EModelX](method.png)

## Environment
`conda env create -f EModelX.yml`  
For EModelX(+AF), you may need to run AlphaFold following <https://github.com/deepmind/alphafold> or get AlphaFold-predicted single-chain structures from AlphaFoldDB (<https://alphafold.ebi.ac.uk/>).  

## Download NN model weights
Download Nn model weights from <https://drive.google.com/file/d/13BKzEBfL0uubYgcJJTGSZ-PkO9hXAB9X/view?usp=drive_link> and build directory `./models`.  

## Minimal Example: Main Chain Modeling for new EM maps

For EModelX:   
`python run_com_modeling.py --protocol=temp_free --EM_map=./inputs/emd_32336.map.gz --fasta=./inputs/7w72.fasta --output_dir=./outputs`  
For EModelX(+AF):   
`python run_com_modeling.py --protocol=temp_flex --EM_map=./inputs/emd_32336.map.gz --fasta=./inputs/7w72.fasta --template_dir=./inputs/templates --output_dir=./outputs`   
, where you can replace `--EM_map` with your target EM map   
, and `--fasta` with your target fasta   
, and `--template_dir`: directory of the template folder, only needed when --protocol == temp_flex   
, and `--output_dir`: the output directory for modeling results  
Notice: If you want to run EModelX(+AF), please place your AlphaFold-predicted single-chain structures in your `--template_dir`， the format should follow our example in ./inputs/templates.  

## Postprocess for all-atom construction and refinement
### Environment
Install pulchra following <https://www.pirx.com/pulchra/> into a directory, e.g. `modules/pulchra304`  
Install phenix following <https://phenix-online.org/> into a directory, e.g. `modules/phenix-1.20.1-4487`  

### Example
For EModelX:   
`python run_com_modeling.py --protocol=temp_free --EM_map=./inputs/emd_32336.map.gz --fasta=./inputs/7w72.fasta --output_dir=./outputs --all_atom  --all_atom --pulchra_path modules/pulchra304/src/pulchra --phenix_act modules/phenix-1.20.1-4487/phenix_env.sh`  
For EModelX(+AF):   
`python run_com_modeling.py --protocol=temp_flex --EM_map=./inputs/emd_32336.map.gz --fasta=./inputs/7w72.fasta --template_dir=./inputs/templates --output_dir=./outputs --all_atom --pulchra_path modules/pulchra304/src/pulchra --phenix_act modules/phenix-1.20.1-4487/phenix_env.sh`   
The modeling results would be `outputs/Result_7w72.pdb`

## Web Server
EModelX's web server is accessible at <https://bio-web1.nscc-gz.cn/app/EModelX>   
