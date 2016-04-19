import bpy, bmesh

class Osm3DBuilding: 
    def __init__(self,verts,tags):
    
        self.height=None
        self.min_height=None
        self.roof_height=None
        self.roof_angle=None
        self.roof_shape="flat" # default roof shape
        self.roof_orientagion="along" # default roof orientation along|across
        self.roof_levels=None # number of Levels in the Roof
        self.length_along=None # building length on the longer side
        self.length_across=None # building width on the shorter side
        self.levels=None
        self.length_across=None
        self.roof_direction=None
        self.heightPerLevel=2.80
        self.min_level=None
        
        self.Walls=None
        self.center=0,0
       
        
        self.source_vertices=verts
        self.tags=tags
        self.vertsCount=len(verts)
        self.__calcWalls()
        
        self.direction={}
        self.__fillDirectionDict()
                
        self.__readTags(tags)
        self.__fillMissingInformation()
        
        
    def __fillMissingInformation(self):
        # If height is not defined, calculate height from levels
        if not self.height:
            if self.levels:  
                self.height = float(self.levels) * self.heightPerLevel
            if self.roof_levels:  
                # Roof-Levels are on top of building levels, so add height here
                self.height += float(self.roof_levels) * self.heightPerLevel
        # if roof height is not defined, calculate from levels
        if not self.roof_height:
            if self.roof_levels:
                self.roof_height=float(self.roof_levels) * self.heightPerLevel
               
        # if min height is not defined, calculate from min level
        if not self.min_height:
            if self.min_level:
                self.min_height=float(self.min_level) * self.heightPerLevel
                
        # if height is still not defined set default
        if not self.height:
            self.height=5
        if not self.min_height:
            self.min_height=0
        if not self.roof_height:
            self.roof_height=0
            self.roof_shape="flat"
         
        
        
    def __readTags(self,tags):
        # read all tags and store relevant information in the class
        if "min_height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            self.min_height,unit = parse_scalar_and_unit(tags["min_height"])

        if "height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            self.height,unit = parse_scalar_and_unit(tags["height"])
                  
        if "building:levels" in tags:
            # If no height is given, calculate height of Building from Levels
            self.levels,unit = parse_scalar_and_unit(tags["building:levels"])
            
        if "roof:levels" in tags:
            # If no height is given, calculate height of Building from Levels
            self.roof_levels,unit = parse_scalar_and_unit(tags["roof:levels"])
            
        if "building:min_level" in tags:
            # If no height is given, calculate height of Building from Levels            
            self.min_level,unit = parse_scalar_and_unit(tags["building:min_level"])
                                
        if "roof:height" in tags:
            self.roof_height,unit = parse_scalar_and_unit(tags["roof:height"])
                    
        if "roof:shape" in tags:
            self.roof_shape = tags["roof:shape"]
            
            
            
    def getVertexTopHeight(self,vertex):
        """This function returns the height of each wall"""
        if self.roof_shape == "gabled":
            pass
        elif self.roof_shape == "skillon":
            pass
        elif self.roof_shape == "gambrel":
            pass
        elif self.roof_shape == "round":
            pass
        elif self.roof_shape == "saltbox":
            pass
            
        # per default the height of the vertices is the same
        # as the building height
        return self.height-self.roof_height

           
           
    def getRoofVertices(self):
        firstBaseVertex=0
        lastBaseVertex=self.vertsCount-1
        firstRoofVertex=self.vertsCount
        lastRoofVertex=(self.vertsCount*2)-1
        
        print("Base: "+str(firstBaseVertex)+"-"+str(lastBaseVertex))
        print("Roof: "+str(firstRoofVertex)+"-"+str(lastRoofVertex))
        
        vertices=[]
        edges=[]
        faces=[]
        
        elif self.roof_shape == "flat":
            roofFace=[]
            for i in range(firstRoofVertex,lastRoofVertex+1):
                roofFace.append(i)
            faces.append(roofFace)
        elif self.roof_shape == "skillon":
            pass
        elif self.roof_shape == "gabled":
            pass
        elif self.roof_shape == "half-hipped":
            pass
        elif self.roof_shape == "hipped":
            pass
        elif self.roof_shape == "pyramidal":
            pass
        elif self.roof_shape == "gambrel":
            pass
        elif self.roof_shape == "mansard":
            pass
        elif self.roof_shape == "dome":
            pass
        elif self.roof_shape == "onion":
            pass
        elif self.roof_shape == "round":
            pass
        elif self.roof_shape == "saltbox":
            pass
        
        print("Vertices: "+str(vertices))
        print("Edges:    "+str(edges))
        print("Faces:    "+str(faces))
        return vertices,edges,faces
          

    def getMesh(self):
        # create the building-mesh and return the vertices used
        vertices=[]
        edges=[]
        faces=[]
        
        # add base vertices
        for vert1 in self.source_vertices:
            vertices.append((vert1[0],vert1[1],self.min_height))
            
        # count Vertices of the base - for each base there will be one top vertex
        BaseVertCount=len(vertices)
            
        # Roof Vertices
        for vert2 in self.source_vertices:
            # Fixed Height 5 until Height Calculation is finished
            roofVert=vert2[0],vert2[1],self.getVertexTopHeight(vert2)
            vertices.append(roofVert)
            
        # ground face
        groundFace=[]
        for i in range(0,BaseVertCount):
            groundFace.append(i)
            
        faces.append(groundFace)
        
        # wall faces except last one
        for wverts in range(1,BaseVertCount):
            wall_face=[]
            wall_face.append(BaseVertCount+wverts-1)
            wall_face.append(BaseVertCount+wverts)
            wall_face.append(wverts)
            wall_face.append(wverts-1)
            faces.append(wall_face)
        # last wall consists of first and last vertex in row
        wall_face=[]
        wall_face.append(BaseVertCount+BaseVertCount-1)
        wall_face.append(BaseVertCount)
        wall_face.append(0)
        wall_face.append(BaseVertCount-1)
        faces.append(wall_face)
        
        # get roof vertices, edges, and faces
        rv,re,rf=self.getRoofVertices()
        
        # append vertex information
        vertices+=rv
        edges+=re
        faces+=rf        
        
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
        
        for i in range(1,length):
            
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
