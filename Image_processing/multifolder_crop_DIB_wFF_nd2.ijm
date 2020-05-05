function cropMulti(openDir, saveDir, filename,m){
	//open(openDir + filename);
	nameL=lengthOf(filename);
	titleName=substring(filename, 0, nameL-7);
	print(titleName);
	mergeString = "";
	mergeStringF = "";
	for (k = 0; k < list.length; k++){

		if (substring(filename, 4, nameL-6) == substring(list[k], 4, nameL-6)){
			open(openDir + list[k]);
			if (!endsWith(list[k], "C1.tif")) {
				//run("Subtract Background...", "rolling=300 stack");
			}
			//print(listF[k]);
			//mergeString = mergeString + " c" + substring(listF[k], 1, 2) + "=" + listF[k]; // not correct!
			mergeString = mergeString + " c" + substring(list[k], nameL-5, nameL-4) + "=" + list[k]; // not correct!
		}
		// choose correct flatfield image
		
	}
	Stack.getDimensions(width, height, channels, slices, frames);
	mergeString = mergeString + " create";
	for (i = 0; i < flatfieldList.length; i++){
			if (!startsWith(flatfieldList[i], "C1.tif")) {
			
				open(flatfieldDir + flatfieldList[i]);
				Stack.getDimensions(wF, hF, cF, sF, fF);
				run("Scale...", "x=- y=- width=" + width +" height=" + height + " interpolation=Bilinear average");
		    	mergeStringF = mergeStringF + " c" + (i+1) + "=" + flatfieldList[i];
		    	makeRectangle((wF-width)/2, (hF-height)/2, width,height);
				run("Crop");
			}	
		}
	mergeStringF = mergeStringF + " create";
	//print(mergeStringF);
	run("Merge Channels...", mergeStringF);
	rename('flatfieldstack');
	print(mergeString);
	run("Merge Channels...", mergeString);
	rename(titleName);
	
	//find bounding box regions for cropping
	run("Duplicate...", "duplicate  channels=1 frames=1");
	selectWindow(titleName+"-1");
	rename("segmentation");
	run("Bandpass Filter...", "filter_large=20 filter_small=2 suppress=None tolerance=5 process");
	run("Convert to Mask", "method=Mean background=Light calculate black");
	//run("Options...", "iterations=1 count=1 black do=Close stack");
	run("Analyze Particles...", "size=1000-Infinity pixel show=Masks exclude clear record add");
	//run("Invert LUT");
	//run("Fill Holes", "stack");
	nROI=roiManager("count");
	for (r = 0; r < nROI; r++){
		selectWindow("segmentation");
		roiManager("Select", r);
		getSelectionBounds(x, y, w, h);
		x=x-100;
		if (x<=0) {
			x=1;
		}
		y=y-100;
		if (y<=0) {
			y=1;
		}
		w=w+200;
		if (w>getWidth()) {
			w=getWidth();
		}
		h=h+200;
		if (w>getHeight()) {
			w=getHeight();
		}
		selectWindow(titleName);
		makeRectangle(x, y, w, h);
		run("Duplicate...", "duplicate");
		//split again and save
		Stack.getDimensions(width, height, channels, slices, frames);
		run("Split Channels");
		
		for (k = 1; k <= channels; k++){
			ch=k+1;
			titleSplit=split(titleName, "-");
			pos=titleSplit[titleSplit.length-1];
			selectWindow("C"+k+"-"+titleName+"-1");
			//waitForUser("click OK to proceed.");
			a=r+1;
			saveAs("Tiff", saveDir + "C"+k+"-"+pos+"-"+a+"_X"+x+"_Y"+y);
			}
		// also crop FF images
		selectWindow('flatfieldstack');
		makeRectangle(x, y, w, h);
		run("Duplicate...", "duplicate");
		//split again and save
		Stack.getDimensions(width, height, channels, slices, frames);
		run("Split Channels");
		
		for (k = 1; k <= channels; k++){
			selectWindow("C"+k+"-flatfieldstack-1");
			saveAs("Tiff", saveDir + "F"+(k+1)+"-"+pos+"-"+a+"_X"+x+"_Y"+y);
			}
		}
	run("Close All");
}

print("cropping started");
run("Close All");
setBatchMode(true);

baseDir=getDirectory("Choose a Directory");
dataSetList = getFileList(baseDir);

for (a = 0; a < dataSetList.length; a++) {
	
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDir=baseDir + dataSet + "\o2-binned_videos\\";
	flatfieldDir = baseDir + dataSet + "\o0-flatfield_images\\";
	saveDir=baseDir + dataSet+ "\o4-cropped_videos\\";
	
	File.makeDirectory(saveDir);
	flatfieldList = getFileList(flatfieldDir);
	list = getFileList(openDir); 
	
	print(openDir);
	print(saveDir);

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
	    	cropMulti(openDir, saveDir, list[i],current);
	    	currentTime=getTime()/1000/60;
			timeElapsed=currentTime-startTime;
			predictedTotalTime=timeElapsed*stackN/current;
			print(current + " of " + stackN + " stacks  was processed");
			print("Time elapsed: " + timeElapsed + " min; Time left: " + (predictedTotalTime-timeElapsed) + " min; Total time: " + predictedTotalTime + " min");
			print("estimated finishing time: "+ hour + ":" + (minute + predictedTotalTime));
		}
	}
}
print("cropping finished");
setBatchMode(false);