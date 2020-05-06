function finalRename(openDir, saveDir, filename){
	// just renaming files
	nameParts=split(filename, ".");
	nameParts=split(nameParts[0], "-");
	refilename = nameParts[1] + "_"+nameParts[2]+"-"+nameParts[0] +".tif";
	//print(refilename);
	File.rename(openDir + filename, saveDir + refilename);
}

print("renaming started");
run("Close All");
setBatchMode(false);

baseDir=getDirectory("Choose a Directory");
dataSetList = getFileList(baseDir);

// loops through all data sets
for (a = 0; a < dataSetList.length; a++) {
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDirB=baseDir + dataSet + "\o5-binaryBF_videos\\";
	openDirF=baseDir + dataSet + "\o4-cropped_videos\\";
	saveDir=baseDir + dataSet + "\o6-final_videos\\";
	File.makeDirectory(saveDir);
	
	listB = getFileList(openDirB); 
	listF = getFileList(openDirF);
	
	// puts all binary in one folder and renames such that the channel number comes last
	// this way, stacks are grouped according to position when sorting by name in the file browser
	for (i = 0; i < listB.length; i++){
    	finalRename(openDirB, saveDir, listB[i]);
	}
	// same for cropped videos, i.e. BF + fluorescence
	for (i = 0; i < listF.length; i++){
		finalRename(openDirF, saveDir, listF[i]);
	}
}
print("renaming finished");