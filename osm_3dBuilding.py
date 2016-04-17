import bpy, bmesh

class Osm3DBuilding: 
    def __init__(self,verts,tags):
        self.source_vertices=verts
        self.tags=tags
        self.vertsCount=len(verts)
        self.__calcWalls()
        self.direction={}
        self.__fillDirectionDict()
        
        self.height=None
        self.min_height=None
        self.roof_height=None
        self.roof_angle=None
        self.roof_shape="flat" # Default roof Shape
        self.roof_orientagion="along" # Default roof Orientation
        self.length_along=None
        self.length_across=None
        self.roof_direction=None
        self.heightPerLevel=2.80
        self.min_level=None
        
        self.__readTags(tags)
        self.__fillMissingInformation()
        
    def __fillMissingInformation(self):
        pass
        
        
    def __readTags(self,tags):
        if "min_height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            self.min_height = parse_scalar_and_unit(tags["min_height"])

        if "height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            self.height,unit = parse_scalar_and_unit(tags["height"])
                  
        if "building:levels" in tags:
            # If no height is given, calculate height of Building from Levels
            self.levels = parse_scalar_and_unit(tags["building:levels"])
            
        if "roof:levels" in tags:
            # If no height is given, calculate height of Building from Levels
            self.roof_levels = parse_scalar_and_unit(tags["roof:levels"])
            
        if "building:min_level" in tags:
            # If no height is given, calculate height of Building from Levels
            self.min_level = parse_scalar_and_unit(tags["building:min_level"])
                                
        if "roof:height" in tags:
            print("roof:height="+str(tags["roof:height"]))
            self.roof_height,unit = parse_scalar_and_unit(tags["roof:height"])
                    
        if "roof:shape" in tags:
            self.roof_shape = tags["roof:shape"]
          

    def getMesh(self):
        vertices=[]
        edges=[]
        faces=[]
        
        for vert1 in self.source_vertices:
            vertices.append(vert1)
            
        # Base Vertices
        countvert=len(vertices)
            
        # Roof Vertices
        for vert2 in self.source_vertices:
            # Fixed Height 5 until Height Calculation is finished
            height=5 if not self.height else self.height
            roofVert=vert2[0],vert2[1],height
            vertices.append(roofVert)
        
        return vertices,edges,faces
                
    def dir(self):
        # Get the main Direction tuple of the Building
        # (maxlength, x vector, y vector)
        maxlen=0
        maxlenkey=None
        for key in self.direction:
            if (self.direction[key][0]>maxlen):
                maxlen=self.direction[key][0]
                maxlenkey=key
                print("newmax="+ str(key)+" "+str(self.direction[key]))
        return self.direction[maxlenkey]
    
    def __fillDirectionDict(self):
        # Dictionary Example:
        #  direction[100]=125,1,0
        #       100  Degrees orientation
        #       125  total length
        #       1,0  Vector of defining direction to calculate exact orientation
        for wall in self.Walls:
            if ( wall[0] in self.direction ):
                newlength=self.direction[wall[0]][0]+wall[1]
                x=self.direction[wall[0]][1]
                y=self.direction[wall[0]][2]
                self.direction[wall[0]]=newlength,x,y
            else:
                self.direction[wall[0]]=wall[1],wall[2],wall[3]
        
    def __calcWalls(self):
        cx=0
        cy=0
        
        wallList=[]
        length=len(self.source_vertices)
        walls=length-1.0
        print( "Source: "+str(self.source_vertices) )
        
        print("Verts: " + str(len(self.source_vertices)))
        for i in range(1,length):
            
            print( "__calcWalls: "+str(i)+"-"+str(length) )
            
            v1=self.source_vertices[i][0] - self.source_vertices[i-1][0]
            v2=self.source_vertices[i][1] - self.source_vertices[i-1][1]
            v3=self.source_vertices[i][2] - self.source_vertices[i-1][2]
            wallList.append(self.__getWallInfoTuple((v1,v2,v3)))
            cx=cx+(self.source_vertices[i][0]/walls)
            cy=cy+(self.source_vertices[i][1]/walls)
        self.Walls=wallList
        self.center=cx,cy
        
    @staticmethod
    def __getWallInfoTuple(wallDirection):
        if (wallDirection[1]<0):
            x=0-wallDirection[0];
            y=0-wallDirection[1];
        else:
            x=0+wallDirection[0];
            y=0+wallDirection[1];
        _length=math.sqrt(math.pow(x,2)+math.pow(y,2))
        _dir=math.trunc(math.atan2(y,x)/math.pi*180/5)*5
        return _dir,_length,wallDirection[0],wallDirection[1]
