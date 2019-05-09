//Pre-processes all brightfield data in the given .czi dataset using bioformats. FINAL VERSION FOR ROSETTE PAPER - 08/05/19

run("Bio-Formats Macro Extensions");

BFPath = "D:\\SlabAssays\\CocultureExpansions\\10_11_17_PilH_WT_YFP_CFP_Expansion\\WT_CFP_PilH_YFP_10_11_17.czi";
outBFPath = "C:\\Users\\olijm\\Desktop\\cellOnEdgeTest\\Chan1_Scene";
outBFFileName = "\\Frame_";
format = ".tif"

channelNo = 0; //The channel index in which the BF image can be found. Note - indexes from zero.

//Couple of options that will be exclusive to a given file
rollingBallSize = 15;
normalizeBlockSize = 40;

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
		run("Subtract Background...", "rolling=" + rollingBallSize + " sliding");

		saveAs("Tiff", outBFPath + (i+1) + outBFFileName + (t+1) + format);

		close();
		
		run("Collect Garbage");
		run("Collect Garbage");
	}

	run("Collect Garbage");
	run("Collect Garbage");
}
