function cropMulti(openDir, saveDir, filename,m){
	nameL=lengthOf(filename);
	titleName=substring(filename, 0, nameL-7);
	print(titleName);
	mergeString = ""; // to combine BF and fluorescence channels
	mergeStringF = ""; // to combine flatfield channels

	// open all channels that match the current position name
	for (k = 0; k < list.length; k++){
		if (substring(filename, 4, nameL-6) == substring(list[k], 4, nameL-6)){
			open(openDir + list[k]);
			mergeString = mergeString + " c" + substring(list[k], nameL-5, nameL-4) + "=" + list[k];
		}
	}
	Stack.getDimensions(width, height, channels, slices, frames);
	mergeString = mergeString + " create";

	// open all relevant flatfield images and scale to the size of the position stacks
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
	
	// merge the flatfield and the position stacks
	run("Merge Channels...", mergeStringF);
	rename('flatfieldstack');
	print(mergeString);
	run("Merge Channels...", mergeString);
	rename(titleName);
	
	// segmentation based on BF contrast thresholding
	run("Duplicate...", "duplicate  channels=1 frames=1");
	selectWindow(titleName+"-1");
	rename("segmentation");
	run("Bandpass Filter...", "filter_large=20 filter_small=2 suppress=None tolerance=5 process");
	run("Convert to Mask", "method=Mean background=Light calculate black");
	run("Analyze Particles...", "size=1000-Infinity pixel show=Masks exclude clear record add");
	
	// loop over all segmented ROIs
	nROI=roiManager("count");
	for (r = 0; r < nROI; r++){
		
		// use the segmented regions to find bounding box regions for cropping (region bounds +100 pixels each side)
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
		if (x+w>getWidth()) {
			w=getWidth()-x;
		}
		h=h+200;
		if (y+h>getHeight()) {
			h=getHeight()-y;
		}
		selectWindow(titleName);
		makeRectangle(x, y, w, h);

		// make a cropped image by duplicating the selected rectangle
		run("Duplicate...", "duplicate");
		
		//split again and save
		Stack.getDimensions(width, height, channels, slices, frames);
		run("Split Channels");
		for (k = 1; k <= channels; k++){
			ch=k+1;
			titleSplit=split(titleName, "-");
			pos=titleSplit[titleSplit.length-1];
			selectWindow("C"+k+"-"+titleName+"-1");
			a=r+1;
			saveAs("Tiff", saveDir + "C"+k+"-"+pos+"-"+a+"_X"+x+"_Y"+y);
		}
		
		// also crop flatfield images
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

// loop through all data sets
for (a = 0; a < dataSetList.length; a++) {
	dataSet=substring(dataSetList[a],0,lengthOf(dataSetList[a])-1) + "\\";
	openDir=baseDir + dataSet + "\o2-binned_videos\\";
	flatfieldDir = baseDir + dataSet + "\o0-flatfield_images\\"; // folder containing the flatfield images
	saveDir=baseDir + dataSet+ "\o4-cropped_videos\\";
	print(openDir);
	print(saveDir);
	
	File.makeDirectory(saveDir);
	flatfieldList = getFileList(flatfieldDir);
	list = getFileList(openDir); 

	// start timer to estimate finishing time
	startTime=getTime()/1000/60;
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	stackN = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){
			stackN +=1;
		}
	}
	current = 0;
	// loop through all positions
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){
			current +=1;
			// this opens a binned BF stack, locates all droplet assemblies by BF contrast thresholding
			// and creates a rectangle around it for cropping
			// then the cropped BF+fluorescence+flatfield images are saved
	    	cropMulti(openDir, saveDir, list[i],current);

	    	// the timer
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