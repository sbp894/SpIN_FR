function calib_fileName= get_lower_calibFile(PicFile)

[DataDir, fileName]=fileparts(PicFile);
picNum= getPicNum(fileName);

NH.calibFile= dir([DataDir filesep 'p*calib*']);

allCalibNums= cellfun(@(x) getPicNum(x), {NH.calibFile.name}');

calibPicNum= find(allCalibNums<picNum, 1, 'last');
curDir= pwd;
cd(DataDir);
calib_fileName= [DataDir filesep getFileName(allCalibNums(calibPicNum))];
cd(curDir);