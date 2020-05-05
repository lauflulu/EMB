function binaryBF(openDir, saveDir, filename){
	open(openDir + filename);
	
	run("Median...", "radius=2 stack");
	run("Duplicate...", "duplicate");
	run("Find Edges", "stack");
	run("Gaussian Blur...", "sigma=5 stack");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=100-Infinity pixel show=Masks stack");
	run("Invert LUT");
	run("Invert", "stack");
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
setBatchMode(false);

baseDir=getDirectory("Choose a Directory");
dataSetList = getFileList(baseDir);

for (a = 0; a < dataSetList.length; a++) {
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDir=baseDir + dataSet +"\o4-cropped_videos\\";
	saveDir=baseDir + dataSet +"\o5-binaryBF_videos\\";
	
	File.makeDirectory(saveDir);
	list = getFileList(openDir); 
	
	// estimate finishing time
	startTime=getTime()/1000/60;
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	stackN = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], ".tif") && startsWith(list[i], "C1")){
			stackN +=1;
		}
	}
	current = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], ".tif") && startsWith(list[i], "C1")){ 
	    	current +=1;
	    	binaryBF(openDir, saveDir, list[i]);
	    	currentTime=getTime()/1000/60;
			timeElapsed=currentTime-startTime;
			predictedTotalTime=timeElapsed*stackN/current;
			print(current + " of " + stackN + " stacks  was processed");
			print("Time elapsed: " + timeElapsed + " min; Time left: " + (predictedTotalTime-timeElapsed) + " min; Total time: " + predictedTotalTime + " min");
			print("estimated finishing time: "+ hour + ":" + (minute + predictedTotalTime));
		}
	}
}
//setBatchMode(false);
print("binary BF finished");


