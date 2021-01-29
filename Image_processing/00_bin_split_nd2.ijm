// converts .nd2 files into individual .tif stacks for each channel with 2x2 binning
function binSplit(openDir, saveDir, filename, series){
	// FYI
	print(filename);
	print(series);
	
	// open current position from .nd2 file
	bioFormatsString = "open=" + openDir + filename + " color_mode=Default view=Hyperstack stack_order=XYCZT series_" + series;
	//print(bioFormatsString);
	run("Bio-Formats Importer", bioFormatsString);

	// bin 2x2 and split channels
	run("Bin...", "x=2 y=2 z=1 bin=Average");
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Split Channels");

	// save each channel individually
	for (k = 0; k < channels; k++){
		if (startsWith(getTitle(), "C")){
			saveAs("Tiff", saveDir + substring(filename, 0, lengthOf(filename)-4) +"-xy"+ series +"-" + substring(getTitle(),0,2) + ".tif");
			close(getTitle());
		}
	}
	
	// safety first
	run("Close All");
	run("Collect Garbage"); 
}

// before starting, the .nd2 file must be placed in baseDir\o1-original_videos\*.nd2
baseDir=getDirectory("Choose a Directory");
openDir=baseDir + "\o1-original_videos\\";
saveDir=baseDir + "\o2-binned_videos\\";
File.makeDirectory(saveDir);
 
print("binning started");
run("Close All");
print(openDir);
print(saveDir);

setBatchMode(true); 
list = getFileList(openDir); // should only have one entry

// get series count, i.e. number of acquired positions
id=openDir + list[0];
run("Bio-Formats Macro Extensions");
Ext.setId(id);
Ext.getSeriesCount(seriesCount);
print("Series count: " + seriesCount);
nSeries=seriesCount;

// start timer to estimate finishing time
startTime=getTime()/1000/60;
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);

// loop through all positions
for (i = 0; i < nSeries; i++){
	
	// binSplit will open one position at a time, bin 2x2 and save individual channels in saveDir
	binSplit(openDir, saveDir, list[0],i+1);
	
	// just the timer. minutes are not modulo(60)
	current = i+1;
	currentTime=getTime()/1000/60;
	timeElapsed=currentTime-startTime;
	predictedTotalTime=timeElapsed*nSeries/current;
	print(current + " of " + nSeries + " stacks  was processed");
	print("Time elapsed: " + timeElapsed + " min; Time left: " + (predictedTotalTime-timeElapsed) + " min; Total time: " + predictedTotalTime + " min");
	print("estimated finishing time: "+ hour + ":" + (minute + predictedTotalTime));
}

setBatchMode(false);
print("binning finished");