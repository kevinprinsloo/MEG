function TrialClass=readClassFile(ClassFileName);

%   Reads a CTF ClassFile and stores information in a structure.
%   The class file allows a user to store a list of classsified trials in a data set.
%   The ClassFile format is defined in document CTF MEG File Formats, PN900-0088.
%   This format is rigid and readClassFile assumes that the ClassFile has the format
%   current in October 2006.

%   Inputs :
%      ClassFileName : marker file including the full path and extension .mrk.
%      trialList : List of trials to read.  Trial numbering : 1,2,...
%                  If omitted or empty, read all markers.

%  Output : Structure array marker.  Output trial numbering starts at 1.
%           See CTF MEG File Formats, (PN900-0088) for the meaning of the structure
%           fields.  A trial mat=y start before the t=0 point, so it is possible to have
%           markers with time<0 (see ds.res4.preTrigPts).

TrialClass=struct([]);

if exist(ClassFileName)~=2
  return     % File doesn't exist.
end

fid=fopen(ClassFileName,'r','ieee-be');

for k=1:5;fgetl(fid);end  % Skip 5 lines (including path info)
nClass=sscanf(fgetl(fid),'%d',1); %Read value and skip to the start of the next non-blank line.
if nClass<=0
  fprintf('readClassFile: File %s has %d classes.\n',nClass);
  return
end

TrialClass=struct('ClassGroupId',[],'Name',char([]),...
  'Comment',char([]),'Color',char([]),'Editable',char([]),'ClassId',[],'trial',[]);

for k=1:nClass
  %  Find the start of the next class identification
  %  There is no need to check for end of file because the loop ends before an attempt
  %  is made to read class nClass+1.
  while ~strcmp('CLASSGROUPID:',fgetl(fid));end
  ClassGroupId=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);
  Name=deblank(fgetl(fid));
  fgetl(fid);
  Comment=deblank(fgetl(fid));
  fgetl(fid);
  Color=deblank(fgetl(fid));
  fgetl(fid);
  Editable=deblank(fgetl(fid));
  fgetl(fid);
  ClassId=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);
  No_of_Trials=sscanf(fgetl(fid),'%d',1);
  fgetl(fid);fgetl(fid);
  if No_of_Trials>0
    trial=reshape(fscanf(fid,'%d',No_of_Trials),1,No_of_Trials);
  else
    trial=[];
  end
  %  Adjust trial numbering so it starts at 1.
  TrialClass(k)=struct('ClassGroupId',ClassGroupId,'Name',Name,...
    'Comment',Comment,'Color',Color,'Editable',Editable,'ClassId',ClassId,...
    'trial',trial+1);
end
fclose(fid);
return
%%%%%%%%%  End of readClassFile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
