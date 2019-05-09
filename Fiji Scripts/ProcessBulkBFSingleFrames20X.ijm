//Pre-processes all the brightfield data in the given BF .czi dataset using bioformats. FINAL VERSION FOR ROSETTE PAPER - 08/05/19

run("Bio-Formats Macro Extensions");

BFPath = "D:\\SlabAssays\\ColonyExpansions\\12_04_17_WT_PilB_PilH_ColonyExpansion\\12_04_17_WT_PilB_PilH_ColonyExpansion.czi";
outBFPath = "C:\\Users\\olijm\\Desktop\\cellOnEdgeTest\\BF_Scene";
outBFFileName = "\\Frame_";
format = ".tif"

channelNo = 0; //The channel index in which the BF image can be found. Note - indexes from zero.

//Couple of options that will be exclusive to a given file
rollingBallSize = 10;

Ext.setId(BFPath);
Ext.getSeriesCount(NoFrames);
for (i=0;i<NoFrames;i++) {
	Ext.setSeries(i);
	Ext.getSizeT(sizeT);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);

	//sizeT = 2;

	if (!(File.isDirectory(outBFPath + (i+1)))) {
		File.makeDirectory(outBFPath + (i+1));
	}

	for (t=0;t<sizeT;t++) {
		Ext.getIndex(0,channelNo,t,currInd);
		Ext.openImage("currBF",currInd);
		
		run("Invert","stack");
		run("Subtract Background...", "rolling=" + rollingBallSize + " disable");

		saveAs("Tiff", outBFPath + (i+1) + outBFFileName + (t+1) + format);

		close();
		
		run("Collect Garbage");
		run("Collect Garbage");
	}

	run("Collect Garbage");
	run("Collect Garbage");
}
