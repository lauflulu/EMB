baseDir=getDirectory("Choose a Directory");
openDir=baseDir + "\o1-original_videos\\";
saveDir=baseDir + "\o2-binned_videos\\";
File.makeDirectory(saveDir);

function binSplit(openDir, saveDir, filename, series){
	print(filename);
	print(series);
	
	bioFormatsString = "open=" + openDir + filename + " color_mode=Default view=Hyperstack stack_order=XYCZT series_" + series;
	print(bioFormatsString);
	run("Bio-Formats Importer", bioFormatsString);
	run("Bin...", "x=2 y=2 z=1 bin=Average");
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Split Channels");
	//print(channels);
	for (k = 0; k < channels; k++){
		if (startsWith(getTitle(), "C")){
			saveAs("Tiff", saveDir + substring(filename, 0, lengthOf(filename)-4) +"-xy"+ series +"-" + substring(getTitle(),0,2) + ".tif");
			close(getTitle());
		}
	}
	run("Close All");
	run("Collect Garbage");
}
 
print("binning started");
run("Close All");
print(openDir);
print(saveDir);

setBatchMode(true); 
list = getFileList(openDir);
// get series count
id=openDir + list[0];
run("Bio-Formats Macro Extensions");

Ext.setId(id);
Ext.getSeriesCount(seriesCount);
print("Series count: " + seriesCount);

nSeries=seriesCount;
// estimate finishing time
startTime=getTime()/1000/60;
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);

for (i = 0; i < nSeries; i++){
	binSplit(openDir, saveDir, list[0],i+1);
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
