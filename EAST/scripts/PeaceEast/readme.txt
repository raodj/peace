Please replace 'east' with the newest EAST executable file, replace 'peace' with the newest PEACE executable.
The source files for 'SepReverseCluster' are located in the folder 'SepReverseCluster_Scr'.
The source files for 'ProcessRC' are located in the folder 'ProcessRC_Scr'.

To run the pipeline, use the command:
python PeaceEastPipeline.py <est file name> -p|-s options
-p: run as parallel
-s: run as serial
options: all the switches for EAST.

If estFile is the name of the input est file, the pipeline will generate two files: estFile.fa.east.contigs and estFile.fa.east.singlets.
