function [kmlStruct,txt] = kml2struct_surface_18(kmlFile)
% kmlStruct = kml2struct(kmlFile)
%
% Import a .kml file as a vector array of shapefile structs, with Geometry, Name,
% Description, Lon, Lat, and BoundaryBox fields.  Structs may contain a mix
% of points, lines, and polygons.
%
% .kml files with folder structure will not be presented as such, but will
% appear as a single vector array of structs.
%
% 

[FID msg] = fopen(kmlFile,'rt');

if FID<0
    error(msg)
end

txt = fread(FID,'uint8=>char')';
fclose(FID);

expr = '<Placemark.+?>.+?</Placemark>';

objectStrings = regexp(txt,expr,'match');

Nos = length(objectStrings);
i=1;
for ii = 1:Nos
    % Find Object Name Field
    bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
    if isempty(bucket)
        name = 'undefined';
    else
        % Clip off flags
        name = regexprep(bucket{1},'<name.*?>\s*','');
        name = regexprep(name,'\s*</name>','');
    end
    
    % Find Object Description Field
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    if isempty(bucket)
        desc = '';
    else
        % Clip off flags
        desc = regexprep(bucket{1},'<description.*?>\s*','');
        desc = regexprep(desc,'\s*</description>','');
    end
    
    geom = 0;
    % Identify Object Type
    if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
        geom = 1;
    elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
        geom = 2;
    elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
        geom = 3;
    end
    
    switch geom
        case 1
            geometry = 'Point';
        case 2
            geometry = 'Line';
        case 3
            geometry = 'Polygon';
        otherwise
            geometry = '';
    end
    

    % Find code Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="code">.+?</SimpleData>','match');
    if isempty(bucket)
    code="0";
    else
    coordStr = regexp(bucket{1},'[0-9]+','match');
    code=string(coordStr);
    end
    
    if code=="2"
        continue
    end
    
    % Find Lanecode {Kind} Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="lanecode">.+?</SimpleData>','match');
    if isempty(bucket)
    lanecode="0";
    else
    coordStr = regexp(bucket{1},'[0-9]+','match');
    lanecode=string(coordStr);
    end
    
    % Find lno Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="lno">.+?</SimpleData>','match');
    if isempty(bucket)
    lno="0";
    else
    coordStr = regexp(bucket{1},'[0-9]+','match');
    lno=string(coordStr);
    end
    
    % Find id {ID} Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="id">.+?</SimpleData>','match');
    if isempty(bucket)
    id="0";
    else
    coordStr = regexp(bucket{1},'>[A-Z0-9]*<','match');
    if isempty(coordStr)
        id="0";
    else
    coordStr = regexp(coordStr{1},'[A-Z0-9]*','match');   
    id=string(coordStr);
    end
    end
    
    % Find r_linkid {R_LinkID} Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="r_linkid">.+?</SimpleData>','match');
    if isempty(bucket)
    r_linkid="0";
    else
    coordStr = regexp(bucket{1},'>[A-Z0-9]*<','match');
    coordStr = regexp(coordStr{1},'[A-Z0-9]*','match');
    r_linkid=string(coordStr);
    end
    
    % Find l_linkid {L_LinkID} Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="l_linkid">.+?</SimpleData>','match');
    if isempty(bucket)
    l_linkid="0";
    else
    coordStr = regexp(bucket{1},'>[A-Z0-9]*<','match');
    coordStr = regexp(coordStr{1},'[A-Z0-9]*','match');
    l_linkid=string(coordStr);
    end
    
    % Find lanetype Property
    bucket = regexp(objectStrings{ii},'<SimpleData name="lanetype">.+?</SimpleData>','match');
    if isempty(bucket)
    lanetype="0";
    else
    coordStr = regexp(bucket{1},'[0-9]+','match');
    lanetype=string(coordStr);
    end
    
    % Find Coordinate Field
    bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
    % Clip off flags
    coordStr = regexprep(bucket{1},'<coordinates.*?>(\s+)*','');
    coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
    % Split coordinate string by commas or white spaces, and convert string
    % to doubles
    coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
    % Rearrange coordinates to form an x-by-3 matrix
    [m,n] = size(coordMat);
    coordMat = reshape(coordMat,3,m*n/3)';
    
    % define polygon in clockwise direction, and terminate
    [Lat, Lon] = poly2ccw(coordMat(:,2),coordMat(:,1));
    if geom==3
        Lon = [Lon;NaN];
        Lat = [Lat;NaN];
    end
    
    % Create structure
    kmlStruct(i).Geometry = geometry;
    kmlStruct(i).Name = name;
    kmlStruct(i).Description = desc;
    kmlStruct(i).Lon = Lon;
    kmlStruct(i).Lat = Lat;
    kmlStruct(i).BoundingBox = [[min(Lon) min(Lat);max(Lon) max(Lat)]];
    kmlStruct(i).lanecode = lanecode;
    kmlStruct(i).id = id;
    kmlStruct(i).r_linkid = r_linkid;
    kmlStruct(i).l_linkid = l_linkid;
    kmlStruct(i).lanecode = lanecode;
    kmlStruct(i).code = code;
    kmlStruct(i).lanetype = lanetype;
    kmlStruct(i).lno = lno;
    i=i+1;
end