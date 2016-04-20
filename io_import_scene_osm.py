# This is the release version of the plugin file io_import_scene_osm_dev.py
# If you would like to make edits, make them in the file io_import_scene_osm_dev.py and the other related modules
# To create the release version of io_import_scene_osm_dev.py, execute:
# python plugin_builder.py io_import_scene_osm_dev.py
bl_info = {
    "name": "Import OpenStreetMap (.osm)",
    "author": "Vladimir Elistratov <vladimir.elistratov@gmail.com> and gtoonstra",
    "version": (1, 1, 0),
    "blender": (2, 7, 4),
    "location": "File > Import > OpenStreetMap (.osm)",
    "description": "Import a file in the OpenStreetMap format (.osm)",
    "warning": "",
    "wiki_url": "https://github.com/vvoovv/blender-geo/wiki/Import-OpenStreetMap-(.osm)",
    "tracker_url": "https://github.com/vvoovv/blender-geo/issues",
    "support": "COMMUNITY",
    "category": "Import-Export",
}

import bpy, bmesh
# ImportHelper is a helper class, defines filename and invoke() function which calls the file selector
from bpy_extras.io_utils import ImportHelper

import os
import math

# see conversion formulas at
# http://en.wikipedia.org/wiki/Transverse_Mercator_projection
# and
# http://mathworld.wolfram.com/MercatorProjection.html
class TransverseMercator:
    radius = 6378137

    def __init__(self, **kwargs):
        # setting default values
        self.lat = 0 # in degrees
        self.lon = 0 # in degrees
        self.k = 1 # scale factor
        
        for attr in kwargs:
            setattr(self, attr, kwargs[attr])
        self.latInRadians = math.radians(self.lat)

    def fromGeographic(self, lat, lon):
        lat = math.radians(lat)
        lon = math.radians(lon-self.lon)
        B = math.sin(lon) * math.cos(lat)
        x = 0.5 * self.k * self.radius * math.log((1+B)/(1-B))
        y = self.k * self.radius * ( math.atan(math.tan(lat)/math.cos(lon)) - self.latInRadians )
        return (x,y)

    def toGeographic(self, x, y):
        x = x/(self.k * self.radius)
        y = y/(self.k * self.radius)
        D = y + self.latInRadians
        lon = math.atan(math.sinh(x)/math.cos(D))
        lat = math.asin(math.sin(D)/math.cosh(x))

        lon = self.lon + math.degrees(lon)
        lat = math.degrees(lat)
        return (lat, lon)
import xml.etree.cElementTree as etree
import inspect, importlib

def prepareHandlers(kwArgs):
    nodeHandlers = []
    wayHandlers = []
    # getting a dictionary with local variables
    _locals = locals()
    for handlers in ("nodeHandlers", "wayHandlers"):
        if handlers in kwArgs:
            for handler in kwArgs[handlers]:
                if isinstance(handler, str):
                    # we've got a module name
                    handler = importlib.import_module(handler)
                if inspect.ismodule(handler):
                    # iterate through all module functions
                    for f in inspect.getmembers(handler, inspect.isclass):
                        _locals[handlers].append(f[1])
                elif inspect.isclass(handler):
                    _locals[handlers].append(handler)
        if len(_locals[handlers])==0: _locals[handlers] = None
    return (nodeHandlers if len(nodeHandlers) else None, wayHandlers if len(wayHandlers) else None)

class OsmParser:
    
    def __init__(self, filename, **kwargs):
        self.nodes = {}
        self.ways = {}
        self.relations = {}
        self.minLat = 90
        self.maxLat = -90
        self.minLon = 180
        self.maxLon = -180
        # self.bounds contains the attributes of the bounds tag of the .osm file if available
        self.bounds = None
        
        (self.nodeHandlers, self.wayHandlers) = prepareHandlers(kwargs)
        
        self.doc = etree.parse(filename)
        self.osm = self.doc.getroot()
        self.prepare()

    # A 'node' in osm:  <node id="2599524395" visible="true" version="1" changeset="19695235" timestamp="2013-12-29T12:40:05Z" user="It's so funny_BAG" uid="1204291" lat="52.0096203" lon="4.3612318"/>
    # A 'way' in osm: 
    #  <way id="254138613" visible="true" version="1" changeset="19695235" timestamp="2013-12-29T12:58:34Z" user="It's so funny_BAG" uid="1204291">
    #  <nd ref="2599536906"/>
    #  <nd ref="2599537009"/>
    #  <nd ref="2599537013"/>
    #  <nd ref="2599537714"/>
    #  <nd ref="2599537988"/>
    #  <nd ref="2599537765"/>
    #  <nd ref="2599536906"/>
    #  <tag k="building" v="yes"/>
    #  <tag k="ref:bag" v="503100000022259"/>
    #  <tag k="source" v="BAG"/>
    #  <tag k="source:date" v="2013-11-26"/>
    #  <tag k="start_date" v="1850"/>
    #  </way>
    #
    # The parser creates two dictionaries: {}
    #   parser.nodes[ node_id ] = { "lat":lat, "lon":lon, "e":e, "_id":_id }
    #   parser.ways[ way_id ] = { "nodes":["12312312","1312313123","345345453",etc..], "tags":{"building":"yes", etc...}, "e":e, "_id":_id }
    # So the way data is stored can seem a little complex
    #
    # The parser is then passed, along with the 'way' or 'node' object of the parser
    # to the handler functions of buildings and highways, where they are used to convert
    # them into blender objects.
    #
    def prepare(self):
        allowedTags = set(("node", "way", "bounds"))
        for e in self.osm: # e stands for element
            attrs = e.attrib
            if e.tag not in allowedTags : continue
            if "action" in attrs and attrs["action"] == "delete": continue
            if e.tag == "node":
                _id = attrs["id"]
                tags = None
                for c in e:
                    if c.tag == "tag":
                        if not tags: tags = {}
                        tags[c.get("k")] = c.get("v")
                lat = float(attrs["lat"])
                lon = float(attrs["lon"])
                # calculating minLat, maxLat, minLon, maxLon
                # commented out: only imported objects take part in the extent calculation
                #if lat<self.minLat: self.minLat = lat
                #elif lat>self.maxLat: self.maxLat = lat
                #if lon<self.minLon: self.minLon = lon
                #elif lon>self.maxLon: self.maxLon = lon
                # creating entry
                entry = dict(
                    id=_id,
                    e=e,
                    lat=lat,
                    lon=lon
                )
                if tags: entry["tags"] = tags
                self.nodes[_id] = entry
            elif e.tag == "way":
                _id = attrs["id"]
                nodes = []
                tags = None
                for c in e:
                    if c.tag == "nd":
                        nodes.append(c.get("ref"))
                    elif c.tag == "tag":
                        if not tags: tags = {}
                        tags[c.get("k")] = c.get("v")
                # ignore ways without tags
                if tags:
                    self.ways[_id] = dict(
                        id=_id,
                        e=e,
                        nodes=nodes,
                        tags=tags
                    )
            elif e.tag == "bounds":
                self.bounds = {
                    "minLat": float(attrs["minlat"]),
                    "minLon": float(attrs["minlon"]),
                    "maxLat": float(attrs["maxlat"]),
                    "maxLon": float(attrs["maxlon"])
                }
        
        self.calculateExtent()

    def iterate(self, wayFunction, nodeFunction):
        nodeHandlers = self.nodeHandlers
        wayHandlers = self.wayHandlers
        
        if wayHandlers:
            for _id in self.ways:
                way = self.ways[_id]
                if "tags" in way:
                    for handler in wayHandlers:
                        if handler.condition(way["tags"], way):
                            wayFunction(way, handler)
                            continue
        
        if nodeHandlers:
            for _id in self.nodes:
                node = self.nodes[_id]
                if "tags" in node:
                    for handler in nodeHandlers:
                        if handler.condition(node["tags"], node):
                            nodeFunction(node, handler)
                            continue

    def parse(self, **kwargs):
        def wayFunction(way, handler):
            handler.handler(way, self, kwargs)
        def nodeFunction(node, handler):
            handler.handler(node, self, kwargs)
        self.iterate(wayFunction, nodeFunction)

    def calculateExtent(self):
        def wayFunction(way, handler):
            wayNodes = way["nodes"]
            for node in range(len(wayNodes)-1): # skip the last node which is the same as the first ones
                nodeFunction(self.nodes[wayNodes[node]])
        def nodeFunction(node, handler=None):
            lon = node["lon"]
            lat = node["lat"]
            if lat<self.minLat: self.minLat = lat
            elif lat>self.maxLat: self.maxLat = lat
            if lon<self.minLon: self.minLon = lon
            elif lon>self.maxLon: self.maxLon = lon
        self.iterate(wayFunction, nodeFunction)

import bpy, bmesh
import bpy, bmesh

def extrudeMesh(bm, thickness):
    """
    Extrude bmesh
    """
    geom = bmesh.ops.extrude_face_region(bm, geom=bm.faces)
    verts_extruded = [v for v in geom["geom"] if isinstance(v, bmesh.types.BMVert)]
    bmesh.ops.translate(bm, verts=verts_extruded, vec=(0, 0, thickness))


def assignMaterials(obj, materialname, color, faces):
    # Get material
    if bpy.data.materials.get(materialname) is not None:
        mat = bpy.data.materials[materialname]
    else:
        # create material
        mat = bpy.data.materials.new(name=materialname)
        mat.diffuse_color = color

    # Assign it to object
    matidx = len(obj.data.materials)
    obj.data.materials.append(mat) 

    for face in faces:
        face.material_index = matidx
import re

def assignTags(obj, tags):
    for key in tags:
        obj[key] = tags[key]


def parse_scalar_and_unit( htag ):
    for i,c in enumerate(htag):
        if not c.isdigit():
            return int(htag[:i]), htag[i:].strip()
    return int(htag), ""

import bpy, bmesh, math, mathutils

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
        self.roof_direction=None
        self.heightPerLevel=2.80
        self.min_level=None
        
        self.Walls=None
        self.center=0,0
       
        
        self.source_vertices=verts
        self.tags=tags
        self.vertsCount=len(verts)
        self.__calcCenter()
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
            
        # Roof Angle
        # translate roof directions into degrees
        if self.roof_direction == "N":
            self.roof_direction = 0/180.0*math.pi
        if self.roof_direction == "NE":
            self.roof_direction = 45/180.0*math.pi
        if self.roof_direction == "E":
            self.roof_direction = 90/180.0*math.pi
        if self.roof_direction == "SE":
            self.roof_direction = 135/180.0*math.pi
        if self.roof_direction == "S":
            self.roof_direction = 180/180.0*math.pi
        if self.roof_direction == "SW":
            self.roof_direction = 225/180.0*math.pi
        if self.roof_direction == "W":
            self.roof_direction = 270/180.0*math.pi
        if self.roof_direction == "NW":
            self.roof_direction = 315/180.0*math.pi
        
        if not self.roof_direction:
            self.roof_direction=self.getRoofDirection()
            
        print("Roof Direction ="+str(self.roof_direction /math.pi*180.0))
            
        self.getWidthLength()
                    
        
    def getWidthLength(self):
        """Calculates the Width and Length of the Building"""    
        # DirectionVec
        self.vec_along  = mathutils.Vector((math.sin(self.roof_direction),math.cos(self.roof_direction),0.0))
        self.vec_across = mathutils.Vector((self.vec_along.y,0.0-self.vec_along.x,0.0))
        
        minLengthAcross=0
        maxLengthAcross=0
        
        minLengthAlong=0
        maxLengthAlong=0
        
        for vert in self.source_vertices:
            along,across = self.getAlongAcrossOfPoint(vert)
            if along < minLengthAlong:
                minLengthAlong=along
            if along > maxLengthAlong:
                maxLengthAlong=along
                
            if across < minLengthAcross:
                minLengthAcross=across
            if across > maxLengthAcross:
                maxLengthAcross=across
              
        self.LengthAlong = maxLengthAlong-minLengthAlong
        self.LengthAcross = maxLengthAcross-minLengthAcross   
        print("along, across",str((self.LengthAlong,self.LengthAcross)))     
        
            
            
    def getAlongAcrossOfPoint(self,point):
        print("center: "+str(self.center))
        diff = mathutils.Vector((self.center[0]-point[0],self.center[1]-point[1],0.0))
        across = diff.dot(self.vec_along)                
        along = diff.dot(self.vec_across)
        return along,across
        
        
            
    def getRoofDirection(self):
        # Get Roof Angle from Longest Side
        len,x,y=self.dir()
        angle=math.atan2(x,y)
        
        print("angle =" + str(angle*180.0/math.pi))
        
        # If roof direction is along, rotate by 90 degrees
        if self.roof_orientagion=="along":
            angle+=(math.pi/2)
            
        if angle>(2*math.pi):
            angle-=(2*math.pi)
        
        return angle
            
        
        
        
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
            
        if "roof:direction" in tags:
            self.roof_direction = tags["roof:direction"]
            
          
    
    
    def getVertexTopHeight(self,vertex):
        """This function returns the height of each wall"""
        if self.roof_shape == "gabled":
            along,across = self.getAlongAcrossOfPoint(vertex)
            factor= abs(across/(self.LengthAcross/2))
            
            return self.height-(self.roof_height*factor)
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
                
        vertices=[]
        edges=[]
        faces=[]
        
        if self.roof_shape == "flat":
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
        wallList=[]
        length=len(self.source_vertices)
        walls=length-1.0
        for i in range(1,length):
            v1=self.source_vertices[i][0] - self.source_vertices[i-1][0]
            v2=self.source_vertices[i][1] - self.source_vertices[i-1][1]
            v3=self.source_vertices[i][2] - self.source_vertices[i-1][2]
            wallList.append(self.__getWallInfoTuple((v1,v2,v3)))
        self.Walls=wallList
        
        
    def __calcCenter(self):
        cx=0
        cy=0
        length=len(self.source_vertices)
        maxX,minX=self.source_vertices[0][0],self.source_vertices[0][0]
        maxY,minY=self.source_vertices[0][1],self.source_vertices[0][1]
        for i in range(1,length):
            if self.source_vertices[i][0] > maxX:
                maxX=self.source_vertices[i][0]
            if self.source_vertices[i][0] < minX:
                minX=self.source_vertices[i][0]
            if self.source_vertices[i][1] > maxY:
                maxY=self.source_vertices[i][1]
            if self.source_vertices[i][1] < minY:
                minY=self.source_vertices[i][1]
        self.center=(maxX+minX)/2.0,(maxY+minY)/2.0
        
        
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


class Buildings:
    @staticmethod
    def condition(tags, way):
        return "building" in tags
    
    @staticmethod
    def handler(way, parser, kwargs):
        wayNodes = way["nodes"]
        numNodes = len(wayNodes)-1 # we need to skip the last node which is the same as the first ones
        # a polygon must have at least 3 vertices
        if numNodes<3: return
        
        if not kwargs["bm"]: # not a single mesh
            tags = way["tags"]
            osmId = way["id"]
            # compose object name
            name = osmId
            if "addr:housenumber" in tags and "addr:street" in tags:
                name = tags["addr:street"] + ", " + tags["addr:housenumber"]
            elif "name" in tags:
                name = tags["name"]
        
        bm = kwargs["bm"] if kwargs["bm"] else bmesh.new()
        verts = []
        for node in range(numNodes):
            node = parser.nodes[wayNodes[node]]
            v = kwargs["projection"].fromGeographic(node["lat"], node["lon"])
            verts.append( bm.verts.new((v[0], v[1], 0)) )
        
        bm.faces.new(verts)
        
        if not kwargs["bm"]:
            tags = way["tags"]
            thickness = 0
            if "height" in tags:
                # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
                thickness,unit = parse_scalar_and_unit(tags["height"])
            else:
                thickness = kwargs["thickness"] if ("thickness" in kwargs) else 0

            # extrude
            if thickness>0:
                extrudeMesh(bm, thickness)
            
            bm.normal_update()
            
            mesh = bpy.data.meshes.new(osmId)
            bm.to_mesh(mesh)
            
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.objects.link(obj)
            bpy.context.scene.update()
            
            # final adjustments
            obj.select = True
            # assign OSM tags to the blender object
            assignTags(obj, tags)

            assignMaterials( obj, "roof", (1.0,0.0,0.0), [mesh.polygons[0]] )
            assignMaterials( obj, "wall", (1,0.7,0.0), mesh.polygons[1:] )


class BuildingParts:
    @staticmethod
    def condition(tags, way):
        return "building:part" in tags
    
    @staticmethod
    def handler(way, parser, kwargs):
        wayNodes = way["nodes"]
        numNodes = len(wayNodes)-1 # we need to skip the last node which is the same as the first ones
        # a polygon must have at least 3 vertices
        if numNodes<3: return
        
        tags = way["tags"]
        if not kwargs["bm"]: # not a single mesh
            osmId = way["id"]
            # compose object name
            name = osmId
            if "addr:housenumber" in tags and "addr:street" in tags:
                name = tags["addr:street"] + ", " + tags["addr:housenumber"]
            elif "name" in tags:
                name = tags["name"]

        min_height = 0
        height = 0
        if "min_height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            min_height,unit = parse_scalar_and_unit(tags["min_height"])

        if "height" in tags:
            # There's a height tag. It's parsed as text and could look like: 25, 25m, 25 ft, etc.
            height,unit = parse_scalar_and_unit(tags["height"])

        bm = kwargs["bm"] if kwargs["bm"] else bmesh.new()
        verts = []
        vertPositions = []
        
        for node in range(numNodes):
            node = parser.nodes[wayNodes[node]]
            v = kwargs["projection"].fromGeographic(node["lat"], node["lon"])
            #verts.append( bm.verts.new((v[0], v[1], min_height)) )
            
            vertPositions.append( (v[0],v[1],min_height) )
            # bm.faces.new(verts)
                
        if not kwargs["bm"]:
            tags = way["tags"]
            
            #bm.normal_update()
            
            mesh = bpy.data.meshes.new(osmId)           
            #bm.to_mesh(mesh)
            mesh.from_pydata(*Osm3DBuilding(vertPositions,tags).getMesh())
            
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.objects.link(obj)
            bpy.context.scene.update()
            
            # final adjustments
            obj.select = True
            # assign OSM tags to the blender object
            assignTags(obj, tags)

class Highways:
    @staticmethod
    def condition(tags, way):
        return "highway" in tags
    
    @staticmethod
    def handler(way, parser, kwargs):
        wayNodes = way["nodes"]
        numNodes = len(wayNodes) # we need to skip the last node which is the same as the first ones
        # a way must have at least 2 vertices
        if numNodes<2: return
        
        if not kwargs["bm"]: # not a single mesh
            tags = way["tags"]
            osmId = way["id"]
            # compose object name
            name = tags["name"] if "name" in tags else osmId
        
        bm = kwargs["bm"] if kwargs["bm"] else bmesh.new()
        prevVertex = None
        for node in range(numNodes):
            node = parser.nodes[wayNodes[node]]
            v = kwargs["projection"].fromGeographic(node["lat"], node["lon"])
            v = bm.verts.new((v[0], v[1], 0))
            if prevVertex:
                bm.edges.new([prevVertex, v])
            prevVertex = v
        
        if not kwargs["bm"]:
            mesh = bpy.data.meshes.new(osmId)
            bm.to_mesh(mesh)
            
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.objects.link(obj)
            bpy.context.scene.update()
            
            # final adjustments
            obj.select = True
            # assign OSM tags to the blender object
            assignTags(obj, tags)
class Naturals:
    @staticmethod
    def condition(tags, way):
        return "natural" in tags
    
    @staticmethod
    def handler(way, parser, kwargs):
        wayNodes = way["nodes"]
        numNodes = len(wayNodes) # we need to skip the last node which is the same as the first ones
    
        if numNodes == 1:
            # This is some point "natural".
            # which we ignore for now (trees, etc.)
            pass

        numNodes = numNodes - 1

        # a polygon must have at least 3 vertices
        if numNodes<3: return
        
        tags = way["tags"]
        if not kwargs["bm"]: # not a single mesh
            osmId = way["id"]
            # compose object name
            name = osmId
            if "name" in tags:
                name = tags["name"]

        bm = kwargs["bm"] if kwargs["bm"] else bmesh.new()
        verts = []
        for node in range(numNodes):
            node = parser.nodes[wayNodes[node]]
            v = kwargs["projection"].fromGeographic(node["lat"], node["lon"])
            verts.append( bm.verts.new((v[0], v[1], 0)) )
        
        bm.faces.new(verts)
        
        if not kwargs["bm"]:
            tags = way["tags"]
            bm.normal_update()
            
            mesh = bpy.data.meshes.new(osmId)
            bm.to_mesh(mesh)
            
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.objects.link(obj)
            bpy.context.scene.update()
            
            # final adjustments
            obj.select = True
            # assign OSM tags to the blender object
            assignTags(obj, tags)

            naturaltype = tags["natural"]
            color = (0.5,0.5,0.5)

            if naturaltype == "water":
                color = (0,0,1)

            assignMaterials( obj, naturaltype, color, [mesh.polygons[0]] )

import bpy, bmesh

def extrudeMesh(bm, thickness):
    """
    Extrude bmesh
    """
    geom = bmesh.ops.extrude_face_region(bm, geom=bm.faces)
    verts_extruded = [v for v in geom["geom"] if isinstance(v, bmesh.types.BMVert)]
    bmesh.ops.translate(bm, verts=verts_extruded, vec=(0, 0, thickness))


def assignMaterials(obj, materialname, color, faces):
    # Get material
    if bpy.data.materials.get(materialname) is not None:
        mat = bpy.data.materials[materialname]
    else:
        # create material
        mat = bpy.data.materials.new(name=materialname)
        mat.diffuse_color = color

    # Assign it to object
    matidx = len(obj.data.materials)
    obj.data.materials.append(mat) 

    for face in faces:
        face.material_index = matidx

class ImportOsm(bpy.types.Operator, ImportHelper):
    """Import a file in the OpenStreetMap format (.osm)"""
    bl_idname = "import_scene.osm"  # important since its how bpy.ops.import_scene.osm is constructed
    bl_label = "Import OpenStreetMap"
    bl_options = {"UNDO"}

    # ImportHelper mixin class uses this
    filename_ext = ".osm"

    filter_glob = bpy.props.StringProperty(
        default="*.osm",
        options={"HIDDEN"},
    )

    ignoreGeoreferencing = bpy.props.BoolProperty(
        name="Ignore existing georeferencing",
        description="Ignore existing georeferencing and make a new one",
        default=False,
    )
    
    singleMesh = bpy.props.BoolProperty(
        name="Import as a single mesh",
        description="Import OSM objects as a single mesh instead of separate Blender objects",
        default=False,
    )

    importBuildings = bpy.props.BoolProperty(
        name="Import buildings",
        description="Import building outlines",
        default=True,
    )

    importNaturals = bpy.props.BoolProperty(
        name="Import naturals",
        description="Import natural outlines",
        default=False,
    )

    importHighways = bpy.props.BoolProperty(
        name="Import roads and paths",
        description="Import roads and paths",
        default=False,
    )

    thickness = bpy.props.FloatProperty(
        name="Thickness",
        description="Set thickness to make OSM building outlines extruded",
        default=0,
    )

    def execute(self, context):
        # setting active object if there is no active object
        if context.mode != "OBJECT":
            # if there is no object in the scene, only "OBJECT" mode is provided
            if not context.scene.objects.active:
                context.scene.objects.active = context.scene.objects[0]
            bpy.ops.object.mode_set(mode="OBJECT")
        
        bpy.ops.object.select_all(action="DESELECT")
        
        name = os.path.basename(self.filepath)
        
        if self.singleMesh:
            self.bm = bmesh.new()
        else:
            self.bm = None
            # create an empty object to parent all imported OSM objects
            bpy.ops.object.empty_add(type="PLAIN_AXES", location=(0, 0, 0))
            parentObject = context.active_object
            self.parentObject = parentObject
            parentObject.name = name
        
        self.read_osm_file(context)
        
        if self.singleMesh:
            bm = self.bm
            # extrude
            if self.thickness>0:
                extrudeMesh(bm, self.thickness)
            
            bm.normal_update()
            
            mesh = bpy.data.meshes.new(name)
            bm.to_mesh(mesh)
            
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.objects.link(obj)
            
            # remove double vertices
            context.scene.objects.active = obj
            bpy.ops.object.mode_set(mode="EDIT")
            bpy.ops.mesh.select_all(action="SELECT")
            bpy.ops.mesh.remove_doubles()
            bpy.ops.mesh.select_all(action="DESELECT")
            bpy.ops.object.mode_set(mode="OBJECT")
            
            bpy.context.scene.update()
        else:
            # perform parenting
            context.scene.objects.active = parentObject
            bpy.ops.object.parent_set()
        
        bpy.ops.object.select_all(action="DESELECT")
        return {"FINISHED"}

    def read_osm_file(self, context):
        scene = context.scene
        
        wayHandlers = []
        if self.importBuildings:
            wayHandlers.append(Buildings)
            wayHandlers.append(BuildingParts)

        if self.importNaturals:
            wayHandlers.append(Naturals)

        if self.importHighways: wayHandlers.append(Highways)
        
        osm = OsmParser(self.filepath,
            # possible values for wayHandlers and nodeHandlers list elements:
            #    1) a string name for the module containing classes (all classes from the modules will be used as handlers)
            #    2) a python variable representing the module containing classes (all classes from the modules will be used as handlers)
            #    3) a python variable representing the class
            # Examples:
            # wayHandlers = [buildings, highways]
            # wayHandlers = [handlers.buildings]
            # wayHandlers = [handlers]
            # wayHandlers = ["handlers"]
            wayHandlers = wayHandlers
        )
        
        if "latitude" in scene and "longitude" in scene and not self.ignoreGeoreferencing:
            lat = scene["latitude"]
            lon = scene["longitude"]
        else:
            if osm.bounds and self.importHighways:
                # If the .osm file contains the bounds tag,
                # use its values as the extent of the imported area.
                # Highways may go far beyond the values of the bounds tag.
                # A user might get confused if higways are used in the calculation of the extent of the imported area.
                bounds = osm.bounds
                lat = (bounds["minLat"] + bounds["maxLat"])/2
                lon = (bounds["minLon"] + bounds["maxLon"])/2
            else:
                lat = (osm.minLat + osm.maxLat)/2
                lon = (osm.minLon + osm.maxLon)/2
            scene["latitude"] = lat
            scene["longitude"] = lon
        
        osm.parse(
            projection = TransverseMercator(lat=lat, lon=lon),
            thickness = self.thickness,
            bm = self.bm # if present, indicates the we need to create as single mesh
        )


# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(ImportOsm.bl_idname, text="OpenStreetMap (.osm)")

def register():
    bpy.utils.register_class(ImportOsm)
    bpy.types.INFO_MT_file_import.append(menu_func_import)

def unregister():
    bpy.utils.unregister_class(ImportOsm)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)

# This allows you to run the script directly from blenders text editor
# to test the addon without having to install it.
if __name__ == "__main__":
    register()


