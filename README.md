# NanoNet


<p align="center"><img src="https://drive.google.com/uc?id=1DdACpv5loaOnKbrIIlJSygUUmt9dRUut" height="250"/></p>

NanoNet - a rapid nanobody modeling tool. 
for citation, please cite our paper: https://www.frontiersin.org/articles/10.3389/fimmu.2022.958584/full

How to run NanoNet from google Colaboratory:

    1. Open the Colab notebook (NanoNet.ipynb) - [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dina-lab3D/NanoNet/blob/main/NanoNet.ipynb).
    2. Select protein type (Nb/mAb heavy chain or TCR VB).
    3. Select input type (sequence or path to a fasta file)
    4. Provide a Nb sequence/fasta.
    5. Select whether or not you want to reconstruct the side chains.
    6. Press the 'Run all' option.

How to run NanoNet locally:

    1. Clone the git repository : git clone "https://github.com/dina-lab3D/NanoNet"
    2. Make sure you have the following libraries installed in your environment:

            - numpy
            - tensorflow (2.4.0 or higher)
            - Bio
            - modeller (requires license - https://salilab.org/modeller/)

    3. If you want to construct the side chains using Scwrl4, make sure you have it installed on your computer (requires license - http://dunbrack.fccc.edu/SCWRL3.php/).

    4. Run the following command (with python 3):

            python NanoNet.py <fasta file path>

            this will produce a backbone + cb pdb named '<record name>_nanonet_backbone_cb.pdb' for each record in the fasta file.

            options:

                    -s : write all the models into a single PDB file, separated with MODEL and ENDMDL (reduces running time when predicting many structures), default is False.
                    -o <output directory> : path to a directory to put the generated models in, default is './NanoNetResults'
                    -m : run side chains reconstruction using modeller, default is False. Output it to a pdb file named '<record name>_nanonet_full_relaxed.pdb'
                    -c <path to Scwrl4 executable>: run side chains reconstruction using scwrl, default is False. Output it to a pdb file named '<record name>_nanonet_full.pdb'
                    -t : use this parameter for TCR V-beta modeling, default is False

Running times for 1,000 structures on a single standard CPU: 

only backbone + Cb - less than 15 seconds (For better preformance use GPU and cuda).
backbone + SCWRL - about 20 minutes. 
backbone + Modeller - about 80 minutes.
