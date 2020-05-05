function excludeBG(openDir, saveDir, filename){
	open(openDir + filename);
	nameL=lengthOf(filename);
	Stack.getDimensions(w, h, channels, slices, frames);

	run("Bandpass Filter...", "filter_large=20 filter_small=2 suppress=None tolerance=5 process");
	run("Convert to Mask", "method=Otsu background=Light calculate black");
	
	//fill holes, also at edges
	run("Invert", "stack");
	run("Canvas Size...", "width="+(w+1)+" height="+(h+1)+" position=Top-Left");
	run("Invert", "stack");
	run("Fill Holes", "stack");	
	run("Canvas Size...", "width=&w height=&h position=Top-Left zero");
	
	run("Invert", "stack");
	run("Canvas Size...", "width="+(w+1)+" height="+(h+1)+" position=Bottom-Right");
	run("Invert", "stack");
	run("Fill Holes", "stack");	
	run("Canvas Size...", "width=&w height=&h position=Bottom-Right zero"); 

	run("Invert", "stack");
	run("Canvas Size...", "width="+(w+1)+" height="+(h+1)+" position=Top-Right");
	run("Invert", "stack");
	run("Fill Holes", "stack");	
	run("Canvas Size...", "width=&w height=&h position=Top-Right zero"); 

	run("Invert", "stack");
	run("Canvas Size...", "width="+(w+1)+" height="+(h+1)+" position=Bottom-Left");
	run("Invert", "stack");
	run("Fill Holes", "stack");	
	run("Canvas Size...", "width=&w height=&h position=Bottom-Left zero"); 
	
	// dilate 10 px
	run("Options...", "iterations=10 count=1 black do=Dilate stack");
	run("Invert", "stack");

	saveAs("Tiff", saveDir + substring(filename,0,nameL-6)+"B1.tif");
	run("Close All");
}

print("exclude started");
run("Close All");
setBatchMode(true);

baseDir=getDirectory("Choose a Directory");
dataSetList = getFileList(baseDir);

for (a = 0; a < dataSetList.length; a++) {
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDir=baseDir + dataSet +"\o6-final_videos\\";
	saveDir=baseDir + dataSet +"\o7-exclude_videos\\";
	
	File.makeDirectory(saveDir);
	list = getFileList(openDir); 
	
	// estimate finishing time
	startTime=getTime()/1000/60;
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	stackN = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){
			stackN +=1;
		}
	}
	
	current = 0;
	for (i = 0; i < list.length; i++){
		//if (endsWith(list[i], ".tif") && startsWith(list[i], "C1")){ 
		if (endsWith(list[i], "C1.tif")){
			current +=1;
	    	excludeBG(openDir, saveDir, list[i]);
	    	currentTime=getTime()/1000/60;
			timeElapsed=currentTime-startTime;
			predictedTotalTime=timeElapsed*stackN/current;
			print(current + " of " + stackN + " stacks  was processed");
			print("Time elapsed: " + timeElapsed + " min; Time left: " + (predictedTotalTime-timeElapsed) + " min; Total time: " + predictedTotalTime + " min");
			print("estimated finishing time: "+ hour + ":" + (minute + predictedTotalTime));
		}
	}
}
print("exclude finished");