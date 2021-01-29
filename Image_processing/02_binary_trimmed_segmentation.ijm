function binaryBF(openDir, saveDir, filename){
	open(openDir + filename);

	// moderate filtering
	run("Median...", "radius=2 stack");
	run("Duplicate...", "duplicate");

	// blurred edge image
	run("Find Edges", "stack");
	run("Gaussian Blur...", "sigma=5 stack");

	// thresholding
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=100-Infinity pixel show=Masks stack");
	run("Invert LUT");
	run("Invert", "stack");

	// binary operations to improve robustness of segmentation
	run("Options...", "iterations=5 count=1 black do=Open stack");
	run("Adjustable Watershed", "tolerance=5 stack");
	run("Invert", "stack");
	run("Skeletonize", "stack");
	run("Invert", "stack");
	run("Options...", "iterations=10 count=1 black do=Erode stack");
	
	saveAs("Tiff", saveDir + replace(filename,"C1","C0"));
	run("Close All");
}

print("binary BF started");
run("Close All");
setBatchMode(true);

baseDir=getDirectory("Choose a Directory");
dataSetList = getFileList(baseDir);

// loop over all data sets
for (a = 0; a < dataSetList.length; a++) {
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDir=baseDir + dataSet +"\o4-cropped_videos\\";
	saveDir=baseDir + dataSet +"\o5-binaryBF_videos\\";
	
	File.makeDirectory(saveDir);
	list = getFileList(openDir); 
	
	// start timer to estimate finishing time
	startTime=getTime()/1000/60;
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	stackN = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], ".tif") && startsWith(list[i], "C1")){
			stackN +=1;
		}
	}
	current = 0;

	// loop through all position stacks
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){
			// creates trimmed binary masks using trimmed segmentation based on edge detection
	    	binaryBF(openDir, saveDir, list[i]);

	    	// the timer
	    	current +=1;
	    	currentTime=getTime()/1000/60;
			timeElapsed=currentTime-startTime;
			predictedTotalTime=timeElapsed*stackN/current;
			print(current + " of " + stackN + " stacks  was processed");
			print("Time elapsed: " + timeElapsed + " min; Time left: " + (predictedTotalTime-timeElapsed) + " min; Total time: " + predictedTotalTime + " min");
			print("estimated finishing time: "+ hour + ":" + (minute + predictedTotalTime));
		}
	}
}
setBatchMode(false);
print("binary BF finished");


