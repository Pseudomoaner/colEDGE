//Pre-processes all the data in the given BF .czi dataset using bioformats.

run("Bio-Formats Macro Extensions");

BFPath = "Name\\Of\\Input\\File.czi";
outBFPath = "Name\\Of\\Output\\Directory\\BF_Scene";
outBFFileName = "\\Frame_";
format = ".tif"

channelNo = 0; //The channel index in which the BF image can be found. Note - indexes from zero.

//Couple of options that will be exclusive to a given file
normalizeBlockSize = 40;

setBatchMode(true);

Ext.setId(BFPath);
Ext.getSeriesCount(NoFrames);
for (i=0;i<NoFrames;i++) {
	Ext.setSeries(i);
	Ext.getSizeT(sizeT);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);

	if (!(File.isDirectory(outBFPath + (i+1)))) {
		File.makeDirectory(outBFPath + (i+1));
	}

	for (t=0;t<sizeT;t++) {
		Ext.getIndex(0,channelNo,t,currInd);
		Ext.openImage("currBF",currInd);

		run("Normalize Local Contrast", "block_radius_x=" + normalizeBlockSize + " block_radius_y=" + normalizeBlockSize + " standard_deviations=3 center stack");
		run("Invert","stack");

		saveAs("Tiff", outBFPath + (i+1) + outBFFileName + (t+1) + format);

		close();
		
		run("Collect Garbage");
		run("Collect Garbage");
	}

	run("Collect Garbage");
	run("Collect Garbage");
}

setBatchMode(false);
