There are three test pipelines in the folder: TestPipeline, TestPipeline-hybrid data, Memory footprint. The first one is used to test single type data, the second one is used to test hybrid data, and the last one is used to test memory usage.

If there is newer PEACE version, please replace the file "peace/peace" with the newest one.

1) TestPipeline
Create five folders under the folder "TestPipeline" with the following names: cap3, ESTAssemblyC++, tgicl, velvet, mira. Put the newest EAST executable file into the folder "ESTAssemblyC++", name it as "Main".

To run the pipeline, use the command:
python EastPipelineParallel.py <batch file name> <output file name>

To set running parameters, open the file "EastPipelineParalle.pl", go to line 21-29, we can change the following parameters:
num_genes: the number of genes used in a test;
min_length: the minimal length of genes in a test;
max_length: the maximal length of genes in a test;
error_rate: integer number, it's only used for sanger test; for short reads test, please set it to [0].
use_quality: whether use quality files or not. 0-not use; 1-use.
coverage: coverage depth.
data_type: there are three values - sanger, 454 and illumina.


num_trials: the number of trials for each test.


2) TestPipeline-hybrid data
Create five folders under the folder "TestPipeline-hybrid data" with the following names: cap3, ESTAssemblyC++, tgicl, velvet, mira. Put the newest EAST executable file into the folder "ESTAssemblyC++", name it as "Main".

To run the pipeline, use the command:
python EastPipelineParallel.py <batch file name> <output file name>

To set running parameters, open the file "EastPipelineParalle.pl", go to line 21-29, we can change the following parameters:
num_genes, min_length, max_length, use_quality, coverage, num_trials are same as those in "TestPipeline".
error_rate: integer number. The error rate is only applied to sanger data.
data_type: only one value "hybrid".

3) Memory footprint
This pipeline generates sanger data and test memory usage of each assembler. 

Create five folders under the folder "Memory footprint" with the following names: cap3, east, tgicl, velvet, mira. Put the newest EAST executable file into the folder "east", name it as "Main".

To run the pipeline, use the command:
python Memory.py

To set running parameters, open the file "Memory.py", go to line 7-9, we can change the following parameters.
minGeneLen: the minimal allowed length of the selected genes;
maxGeneLen: the maximal allowed length of the selected genes;
error rate: the error rate. Note it is a float number.
Other parameters can be changed at line 170, main(coverage, uniqueId, num_genes).
