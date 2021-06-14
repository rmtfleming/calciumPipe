Code for the paper:

Neuronal hyperactivity in a LRRK2-G2019S cellular model of Parkinson’s Disease

Edinson Lucumi Moreno1,2∗, Siham Hachi1∗, Sarah L. Nickels1,
Khalid I.W. Kane1, Masha Moein1, Jens C. Schwamborn1,
Alex Skupin1, Pieter Vanden Berghe3, Ronan M.T. Fleming1,4,5†

1Luxembourg Centre for Systems Biomedicine, University of Luxembourg, 6 avenue du Swing, L-4367 Belvaux, Luxembourg.
2Department of Neurology, Brigham and Women’s Hospital, Harvard Medical School, Harvard University, 75 Francis St, Boston, USA.
3Laboratory for Enteric Neuroscience, Translational Research in Gastrointestinal Disorders, Department of Clinical and Experimental Medicine, University of Leuven, Leuven Herestraat 49 3000, Leuven,Belgium.
4Division of Systems Biomedicine and Pharmacology, Leiden Academic Centre for Drug Research, Leiden University, Einsteinweg 55, Leiden, The Netherlands.
5School of Medicine, National University of Ireland, University Rd, Galway, Ireland.

Abstract
Monogenic Parkinson’s Disease can be caused by a mutation in the leucine-rich repeat kinase 2
(LRRK2) gene, causing a late-onset autosomal dominant inherited form of Parkinson’s Disease. The
function of the LRRK2 gene is incompletely understood, but several in vitro studies have reported that
LRRK2-G2019S mutations affect neurite branching, calcium homeostasis and mitochondrial function,
but thus far, there have been no reports of effects on electrophysiological activity. We assessed the
neuronal activity of induced pluripotent stem cell derived neurons from Parkinson’s Disease patients
with LRRK2-G2019S mutations and isogenic controls. Neuronal activity of spontaneously firing neuronal
populations was recorded with a fluorescent calcium-sensitive dye (Fluo-4) and analysed with a
novel image analysis pipeline that combined semi-automated neuronal segmentation and quantification
of calcium transient properties. Compared with controls, LRRK2-G2019S mutants have shortened
inter-spike intervals and an increased rate of spontaneous calcium transient induction.


Required library: SPAMS

Download the matlab interface from: http://spams-devel.gforge.inria.fr/downloads.html

- Unzip the file

- Compile the software

	- if you are using Windows open in Matlab the file compile.m, then inside this file change the variable compiler to 'mex' instead of gcc then run

	- if you are using Linux just run the script

===============================================================
			Running the pipeline
===============================================================

I - Segmentation and spike train generation


- the functions are ordered in folders by analysis step

- in the folder electrophysiologicalAnalysis open the file startupPipeline.m, modify the paths according to your machine and run

- to run the segmentation and calcium transient analysis 
	- got to the folder pipeline.m open the file intialise_all.m and tune the parameters
	- open the file runPipeline.m and run it
		- The algorithm will create 
			- a "mat" folder containing a MAT file of the data 
			- a "results" folder containing the segmented image, the traces and spike trains

- Organize resulting folders, wells, replicates, cell-lines etc. Put all wells from the same replicate in one folder and add cell-line name e.g. K7_wellE7

 	ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ

II - Statistical Analysis


- got to folder statisticalAnalysis
	- open get_signalProfiles.m
	- modify in the first two lines the key words of folder names
	- run script comparison_WT_MUT.m


