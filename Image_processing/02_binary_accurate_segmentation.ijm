// requires https://imagejdocu.tudor.lu/doku.php?id=plugin:segmentation:adjustable_watershed:start

function binaryBF(openDir, saveDir, filename){
	
	open(openDir + filename);
	rename("main image");
	mainImage=getTitle();

	// moderate filtering
	run("Median...", "radius=2 stack");

	// save a copy of the main image
	run("Duplicate...", "duplicate");
	rename("bilayer image");
	bilayerImage=getTitle();
	
	// create a binary mask by BF contrast thresholding the main image. this is the outline image
	selectWindow(mainImage);
	run("Bandpass Filter...", "filter_large=20 filter_small=2 suppress=None tolerance=5 process");
	run("Convert to Mask", "method=Mean background=Light calculate black");
	run("Analyze Particles...", "size=3000-Infinity show=Masks exclude pixel stack");
	run("Invert LUT");
	run("Fill Holes", "stack");
	run("Outline", "stack");
	mainImage=getTitle();
	// save a copy of the BF contrast binary image
	run("Duplicate...", "duplicate");
	rename("filter image");
	filterImage=getTitle();
	
		// process filter image. this is the main image but dilated. it is later used to trim off the outer droplet edges from the bilayer image
		selectWindow(filterImage);
		run("Options...", "iterations=15 count=1 black do=Dilate stack");
		run("Invert", "stack");
		filterImage=getTitle();
		
		// process bilayer image. segmentation based on thresholding the edge detection image
		selectWindow(bilayerImage);
		run("Find Edges", "stack");
		run("Gaussian Blur...", "sigma=5 stack");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		run("Analyze Particles...", "size=100-Infinity pixel show=Masks stack");
		run("Invert LUT");
		run("Invert", "stack");
		// binary operations to improve robustness
		run("Options...", "iterations=5 count=1 black do=Open stack");
		run("Adjustable Watershed", "tolerance=5 stack");
		run("Invert", "stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=100 count=7 black do=Erode stack");
		bilayerImage=getTitle();

		// trim off outer droplet edges from bilayer image
		imageCalculator("Multiply create 32-bit stack", bilayerImage,filterImage);
		run("8-bit");
		bilayerImage=getTitle();

	// combine the bilayer and the outline image	
	imageCalculator("OR create 32-bit stack", mainImage,bilayerImage);
	run("8-bit");
	run("Invert", "stack");
	// connect bilayers with the outline
	run("Adjustable Watershed", "tolerance=5 stack");
	run("Invert", "stack");
	// get rid of loose ends
	run("Options...", "iterations=100 count=7 black do=Erode stack");
	run("Invert", "stack");

	saveAs("Tiff", saveDir + replace(filename,"C1","B2"));
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
	openDir=baseDir + dataSet +"\o6-final_videos\\";
	saveDir=baseDir + dataSet +"\o9-size_videos\\";
	
	File.makeDirectory(saveDir);
	list = getFileList(openDir); 
	
	// starts timer to estimate finishing time
	startTime=getTime()/1000/60;
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	stackN = 0;
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){
			stackN +=1;
		}
	}
	current = 0;

	// loops through all position stacks
	for (i = 0; i < list.length; i++){
		if (endsWith(list[i], "C1.tif")){ 
	    	// creates binary videos based on accurate segmentation, using a combination of BF contrast and edge detection
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

