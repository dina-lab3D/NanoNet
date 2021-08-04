# NanoNet

NanoNet - a rapid nanobody modeling tool


How to run NanoNet from google Colaboratory:

    1. Just make a copy of the Colab notebook (NanoNet.ipynb) in your drive.
    2. select protein type (Nb/mAb heavy chain or TCR VB).
    3. select input type (sequence or path to a fasta file)
    4. Provide Nb sequence/fasta.
    5. Press the 'Run all' option.

How to run NanoNet locally:

    1. Clone the git repository : git clone "https://github.com/dina-lab3D/NanoNet"
    2. Make sure you have the following libraries installed in your environment:

            - numpy
            - tensorflow (2.4.0 or higher)
            - Bio

    3. If you want to construct the side chains, make sure you have 'pulchra' installed on your computer (https://www.pirx.com/pulchra/)

    4. Run the following command (with python 3):

            python NanoNet.py <fasta file path>

            this will produce a Ca pdb named '<record name>_nanonet_ca.pdb' for each record in the fasta file.

            options:

                    -n <path to trained NanoNet> : a path to the trained model, default is 'NanoNet'.
                    -s : write all the models into a single PDB file, separated with MODEL and ENDMDL (reduces running time when predicting many structures), default is False.
                    -o <output directory> : path to a directory to put the generated models in, default is './NanoNetResults'
                    -r : reconstruct side chains with pulchra, default is False
                    -p <pulchra> : path to pulchra executable, in order to reconstruct the side chains, default is 'pulchra'


Running time: under 20 seconds for 5,000 models written into a single PDB file, or under a minute written into separate PDB files on a standart CPU computer. for better preformance use GPU and cuda.
