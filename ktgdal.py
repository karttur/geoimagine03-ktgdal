'''
Created on 10 Mar 2021

@author: thomasgumbricht
'''

# Standard library imports

import os

from sys import exit

# Third party imports

from osgeo import gdal

from gdalconst import *

# Package application imports

from geoimagine.gis import GetRasterMetaData

from geoimagine.ktgdal import ogr2ogr

import geoimagine.support as supp

#from geoimagine.gis import  MjProj

import geoimagine.support.karttur_dt as mj_dt 
from plotnine.tests.test_geom_ribbon_area import width


def MakeMosaic(tileL, dstFPN, xRes=False, yRes=False, resampleAlg=False, extent=False):
    '''
    '''
    
    # Create text file listing the tiles to mosaic
    
    dstFP, mosaicFN = os.path.split(dstFPN)
    
    tileListFN = '%s-tilelist.txt' %(os.path.splitext(mosaicFN)[0])
    
    tileListFPN = os.path.join(dstFP, tileListFN)
    
    f = open(tileListFPN,'w')
     
    for tile in tileL:
        
        writeln = '%(t)s\n' %{'t':tile}
        
        f.write(writeln)
        
    f.close()
    
    #replace the target extension with vrt
    vrtFN = '%s-full.vrt' %(os.path.splitext(mosaicFN)[0])
    
    vrtFPN = os.path.join(dstFP, vrtFN)
    
    # See https://gdal.org/python/osgeo.gdal-module.html#BuildVRTOptions

    if xRes or resampleAlg:
        
        kwargs = {}
        
        if xRes:
        
            kwargs['xRes'] = xRes
            
            kwargs['yRes'] = yRes
            
        if resampleAlg:
        
            kwargs['resampleAlg'] = resampleAlg
            
    if extent:
        
        kwargs['outputBounds'] = extent
            
        ds = gdal.BuildVRT(vrtFPN, 
                        tileL,
                        **kwargs)
        
    else:
        
        ds = gdal.BuildVRT(vrtFPN, 
                        tileL)

    ds = None
    
    return vrtFPN

class GDALinternal():
    '''
    '''
    
    def __init__(self):
        ''' Initiation only
        '''
        
        self._SetConstantsAndDicts()
    
    def _SetConstantsAndDicts(self):
        '''
        '''
        

        self.dstSRS = 'EPSG:%(epsg)d' %{'epsg':self.pp.procsys.dstepsg}
            
        self.resampleD = {'near':gdal.GRA_NearestNeighbour, 'bili':gdal.GRA_Bilinear,
                     'cubi':gdal.GRA_Cubic,'aver':gdal.GRA_Average}
        
        self.cellTypeD = {'byte':'gdalconst.GDT_Byte','etc':'gdalconst.GDT_Byte'}
        
        self.gdalFormatD = {'.tif':'GTiff', '.vrt':'VRT'}
    
    def _GdalDemTPI(self, dstFPN, srcFPN):
        ''' TPI - see https://gdal.org/python/osgeo.gdal-module.html#DEMProcessing
        '''
        
        kwargs = {'format': 'GTiff'}
            
        if self.pp.process.parameters.compute_edges:
            
            kwargs['computeEdges'] = True
                
        ds = gdal.DEMProcessing(dstFPN, srcFPN, 'TPI', **kwargs)
        
        ds = None
            
    def _GdalDemTRI(self, dstFPN, srcFPN):
        ''' TPI - see https://gdal.org/python/osgeo.gdal-module.html#DEMProcessing
        '''
        
        kwargs = {'format': 'GTiff'}
            
        if self.pp.process.parameters.compute_edges:
            
            kwargs['computeEdges'] = True
            
        if self.pp.process.parameters.algorithm != 'default':
            
            kwargs['alg'] = self.pp.process.parameters.algorithm
            
        if hasattr(self.pp.process.parameters, 'alg'):
            
            kwargs['alg'] = self.pp.process.parameters.alg
            
        if self.verbose:
            
            infostr = '        Analysing TRI with params: ', kwargs
                
            print (infostr)
                
        ds = gdal.DEMProcessing(dstFPN, srcFPN, 'TRI', **kwargs)
        
        ds = None
        
    def _GdalDemHillshade(self, dstFPN, srcFPN):
        ''' TPI - see https://gdal.org/python/osgeo.gdal-module.html#DEMProcessing
        '''
        
        kwargs = {'format': 'GTiff', 'azimuth':self.pp.process.parameters.azimuth,
                  'altitude':self.pp.process.parameters.altitude,'zFactor':self.pp.process.parameters.zfac}
            
        if self.pp.process.parameters.compute_edges:
            
            kwargs['computeEdges'] = True
            
        if self.pp.process.parameters.algorithm != 'default':
            
            kwargs['alg'] = self.pp.process.parameter.algorithm
            
        if hasattr(self.pp.process.parameters, 'alg'):
            
            kwargs['alg'] = self.pp.process.parameters.alg
            
        if self.verbose:
            
            infostr = '        Analysing Hillshade with params: %s' %(kwargs)
                
            print (infostr)
            
            print ('        src:',srcFPN)
            
            print ('        dst:',dstFPN)
            
        ds = gdal.DEMProcessing(dstFPN, srcFPN, 'hillshade', **kwargs)
        
        ds = None
        
    def _GdalDemSlope(self, dstFPN, srcFPN):
        ''' Slope - see https://gdal.org/python/osgeo.gdal-module.html#DEMProcessing
        '''
        
        kwargs = {'format': 'GTiff'}
            
        if self.pp.process.parameters.compute_edges:
            
            kwargs['computeEdges'] = True
            
        if self.pp.process.parameters.algorithm != 'default':
            
            kwargs['alg'] = self.pp.process.parameter.algorithm
            
        if hasattr(self.pp.process.parameters, 'alg'):
            
            kwargs['alg'] = self.pp.process.parameters.alg
                
        if self.pp.process.parameters.trigonometric:
            
            kwargs['slopeFormat'] = 'degree'
            
        else:
            
            kwargs['slopeFormat'] = 'percent'
            
        if self.verbose:
            
            infostr = '        Analysing slope with params: ', kwargs
                
            print (infostr)
                
        ds = gdal.DEMProcessing(dstFPN, srcFPN, 'slope', **kwargs)
        
        ds = None
            
    def _GdalDemAspect(self, dstFPN, srcFPN):
        ''' Slope - see https://gdal.org/python/osgeo.gdal-module.html#DEMProcessing
        '''
        
        kwargs = {'format': 'GTiff'}
        
            
        if self.pp.process.parameters.compute_edges:
            
            kwargs['computeEdges'] = True
            
        if self.pp.process.parameters.algorithm != 'default':
            
            kwargs['alg'] = self.pp.process.parameter.algorithm
            
        if hasattr(self.pp.process.parameters, 'alg'):
            
            kwargs['alg'] = self.pp.process.parameters.alg
                
        if self.pp.process.parameters.trigonometric:
            
            kwargs['slopeFormat'] = 'degree'
            
        else:
            
            kwargs['slopeFormat'] = 'percent'
            
        if self.pp.process.parameters.zero_for_flat:
        
            kwargs['zeroForFlat'] = True
            
        if self.pp.process.parameters.trigonometric:
        
            kwargs['trigonometric'] = True
            
        if self.verbose:
            
            infostr = '        Analysing aspect with params: ', kwargs
                
            print (infostr)
                
        ds = gdal.DEMProcessing(dstFPN, srcFPN, 'aspect', **kwargs)
        
        ds = None
        
    def _GdalTranslate(self, dstFPN, srcFPN):
        '''
        '''
        gdalFormatD = {'.tif':'GTiff', '.vrt':'VRT'}
        
        resampleD = {'near':gdal.GRA_NearestNeighbour, 'bili':gdal.GRA_Bilinear,
                     'cubi':gdal.GRA_Cubic,'aver':gdal.GRA_Average}
        
        gdalFormat = gdalFormatD[os.path.splitext(dstFPN)[1]]
         
        kwargs = {'format': gdalFormat}
        
        if hasattr(self.pp.process.parameters, 'resample'):
            
            print (self.pp.process.parameters.resample[0:4])
                
            kwargs['resampleAlg'] = resampleD[ self.pp.process.parameters.resample[0:4] ]
          
        if hasattr(self.pp.process.parameters, 'dst_ulx'):
             
            if self.pp.process.parameters.dst_ulx or self.pp.process.parameters.dst_lry or self.pp.process.parameters.dst_lrx or self.pp.process.parameters.dst_uly:
                
                kwargs['projWin'] = (self.pp.process.parameters.dst_ulx, self.pp.process.parameters.dst_uly, 
                                          self.pp.process.parameters.dst_lrx, self.pp.process.parameters.dst_lry)
        
        if hasattr(self.pp.process.parameters, 'xsize'): 
               
            if self.pp.process.parameters.xsize or self.pp.process.parameters.ysize:
                
                spatialRef, srcLayer = GetRasterMetaData(srcFPN)
                            
                if self.pp.process.parameters.xsize == 0:
                    
                    xsize = srcLayer.cols
                    
                else:
                    
                    xsize = self.pp.process.parameters.xsize
                    
                if self.pp.process.parameters.ysize == 0:
                    
                    ysize = srcLayer.lins
                    
                else:
                    
                    ysize = self.pp.process.parameters.ysize
                                
                kwargs['srcWin'] = (self.pp.process.parameters.xoff,self.pp.process.parameters.yoff,
                                    xsize,ysize)
                
        if hasattr(self.pp.process.parameters, 'width'):
            
            kwargs['width'] = self.pp.process.parameters.width
            
            kwargs['height'] = self.pp.process.parameters.height
            
        elif hasattr(self.pp.process.parameters, 'tr_xres') and self.pp.process.parameters.tr_xres:
                
            kwargs['xRes'] = self.pp.process.parameters.tr_xres
                
            kwargs['yRes'] = self.pp.process.parameters.tr_yres
            
        if hasattr(self.pp.process.parameters, 'src_min'):
            
            if self.pp.process.parameters.src_min or self.pp.process.parameters.src_max:
                
                kwargs['scaleParams'] = (self.pp.process.parameters.src_min, self.pp.process.parameters.src_max, 
                                        self.pp.process.parameters.dst_min, self.pp.process.parameters.dst_max)
                                     
        if hasattr(self.pp.process.parameters, 'exponent') and self.pp.process.parameters.exponent:
            
            kwargs['exponents'] = [self.pp.process.parameters.exponent]
        
        if hasattr(self.pp.process.parameters, 'dstEPSG') and  self.pp.process.parameters.dstEPSG:
            
            kwargs['outputSRS'] = self.pp.process.parameters.dstEPSG
            
        print ('dstFPN', dstFPN)
        print ('srcFPN', srcFPN)
        print ('translatekwargs', kwargs)
        if not os.path.exists(srcFPN):
            SNULLE
        
        #ds = gdal.Open(srcFPN)
            
        ds = gdal.Translate(dstFPN, srcFPN, 
                            **kwargs)
        
        ds = None
       
    def _GdalWarp(self,dstFPN, srcFPN): 
        '''
        '''
        
        dstSRS = 'EPSG:%(epsg)d' %{'epsg':self.pp.procsys.dstepsg}
        
        gdalFormatD = {'.tif':'GTiff', '.vrt':'VRT'}
        
        resampleD = {'near':gdal.GRA_NearestNeighbour, 'bili':gdal.GRA_Bilinear,
                     'cubi':gdal.GRA_Cubic,'aver':gdal.GRA_Average}
        
        gdalFormat = gdalFormatD[os.path.splitext(dstFPN)[1]]
         
        
        kwargs = {'format': gdalFormat,
                  'dstSRS': dstSRS,
                  'resampleAlg': gdal.GRA_NearestNeighbour}
        
        kwargs['outputBounds'] = (self.pp.process.parameters.dst_ulx, self.pp.process.parameters.dst_lry, 
                    self.pp.process.parameters.dst_lrx, self.pp.process.parameters.dst_uly)
        
        if hasattr(self.pp.process.parameters, 'tr_xres') and self.pp.process.parameters.tr_xres:
                
            kwargs['xRes'] = self.pp.process.parameters.tr_xres
                
            kwargs['yRes'] = self.pp.process.parameters.tr_yres
        
        if hasattr(self.pp.process.parameters, 'resample'):
            
            kwargs['resampleAlg'] = resampleD[self.pp.process.parameters.resample.lower()[0:4]]
               
        if hasattr(self.pp.process.parameters, 'cellnull'):
            
            kwargs['dstNodata'] = self.pp.process.parameters.cellnull
  
        if hasattr(self.pp.process.parameters, 'errorthreshold') and  self.pp.process.parameters.errorthreshold > 0:
                                
            kwargs['errorTHreshold'] = self.pp.process.parameters.errorthreshold
                
        print ('warpkwargs', kwargs)        
          
        ds = gdal.Warp(dstFPN, 
                    srcFPN,
                    **kwargs)

        ds = None
            
class GDALexternal:
    '''
    '''
    
    def __init__(self):
        ''' Initiation only
        '''
        
        pass
    
    def _GdalDemTPI(self):
        ''' TPI
        '''
        
        pass
        
class GDALmosaicAdjacentTiles:
    ''' Mosaic adjacent tiles for neighbor processing
    ''' 
    
    def __init__(self):
        ''' Initiation only
        '''
        
        pass
    
    def _GetAdjacenTiles(self,locus, datum, comp):
        '''
        '''
        self.neighborTiles = [[False,False,False],[False,1,False],[False,False,False]]
        
        if self.pp.procsys.srcsystem == 'modis':
            pass
        
        elif self.pp.procsys.srcsystem[0:4] == 'ease':
            
            x = int(locus[1:3]); y = int(locus[4:6])
            
            centertile = (x,y)
            
            paramL = ['xtile','ytile']
          

        adjacentTiles = self.session._SelectAdjacentTiles(self.pp.defregion, centertile, paramL, self.pp.procsys.srcsystem, 'regions' )
            
        self.srcTileL = []
        
        for tile in adjacentTiles:
                        
            if self.pp.procsys.srcsystem == 'modis':
            
                pass
            
            elif self.pp.procsys.srcsystem[0:4] == 'ease':
                
                tcol = tile[0]-centertile[0]+1
                
                trow = centertile[1]-tile[1]+1
                
                print (adjacentTiles)
                
                print ('centretile',centertile)
                print (tile)
                
                print ('tcol',tcol)
                print ('trow',trow)
                
                
                self.neighborTiles[trow][tcol] = tile
                
                tileLocus = supp.ConvertXYinteger(tile[0],tile[1])
                
                
                
                
                
                tileFPN = self.pp.srcLayerD[tileLocus['prstr']][datum][comp].FPN
                
                if not os.path.exists(tileFPN):
                    
                    exitstr = 'Exiting - missing tile for mosaicAdjacentTiles;\n    %s' %(tileFPN)
                
                    exit(exitstr)
                    
                self.srcTileL.append(tileFPN)
                
        print (self.neighborTiles)
            

                
    def _GetMosaicCorners(self,srcFPN):
        ''' Get the corners (edges) of the mosaic to create
        '''
        
        spatialRef, srcLayer = GetRasterMetaData(srcFPN)
                        
        self.xsize = srcLayer.cols
                        
        self.ysize = srcLayer.lins
        
        cellsize = srcLayer.cellsize
        
        minx, miny, maxx, maxy = srcLayer.bounds
        
        overlap = self.pp.process.parameters.overlap
        
        print (self.neighborTiles)
        
        if any([item[0] for item in self.neighborTiles]):
        
            minx -= overlap*cellsize
            
            self.xsize += overlap
            
        if any (self.neighborTiles[2]):
        
            miny -= overlap*cellsize
            
            self.ysize += overlap
            
        if any([item[2] for item in self.neighborTiles]):
        
            maxx += overlap*cellsize
            
            self.xsize += overlap
            
        if any (self.neighborTiles[0]):
        
            maxy += overlap*cellsize
            
            self.ysize += overlap
        
        self.extent = (minx, miny, maxx, maxy)
        
    def _TranslateMosaic(self,dstLayer):
    
        xRes = self.pp.process.parameters.tr_xres
        
        yRes =  self.pp.process.parameters.tr_yres
        
        resampleAlg = self.resampleD[self.pp.process.parameters.resample.lower()[0:4]]
            
        vrtFPN = MakeMosaic(self.srcTileL, 
                             dstLayer.FPN,
                             xRes,yRes,
                             resampleAlg,
                             self.extent
                             ) 
                
    def _MosaicAdjacentTiles(self):
        ''' Get adjacent tiles
        '''
        
        for datum in self.pp.srcPeriod.datumL:
            
            for comp in self.pp.srcCompL:
                
                for locus in self.pp.srcLayerD:
                    
                    srcLayer = self.pp.srcLayerD[locus][datum][comp]
                    
                    dstLayer = self.pp.dstLayerD[locus][datum][comp]
                    
                    if dstLayer._Exists() and not self.pp.process.overwrite:
                        
                        continue
            
                    self._GetAdjacenTiles(locus, datum ,comp)
                    
                    self._GetMosaicCorners(srcLayer.FPN)
                                        
                    self._TranslateMosaic(dstLayer)
                    

class ProcessGDAL(GDALinternal,GDALmosaicAdjacentTiles):
    '''class for processing using GDAl
    '''   
    
    def __init__(self, pp, session):
        ''' Expects the pp processing instructions and an open postgres database session
        '''
                
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose)

        print ('        ProcessGDAL',self.pp.process.processid) 
        
        # Initiare the processes
        
        GDALinternal.__init__(self)
        
        GDALmosaicAdjacentTiles.__init__(self)

        # Direct to sub-processes
        
        if self.pp.process.processid.lower() == 'tileancillaryregion':
            
            self._SetConstantsAndDicts()
            
            self._TileRegion()
            
        elif self.pp.process.processid.lower() == 'mosaictiles':
            
            self._SetConstantsAndDicts()
            
            self._MosaicTiles()
            
        elif self.pp.process.processid.lower() == 'mosaicadjacenttiles':
             
            self._MosaicAdjacentTiles()
            
        elif self.pp.process.processid.lower() == 'reprojectancillaryregion':
            
            self._SetConstantsAndDicts()
            
            self._ReprojectRegion()
            
        elif self.pp.process.processid.lower() == 'translateregion':
            
            self._SetConstantsAndDicts()
            
            self._TranslateTilesRegions()
            
        elif self.pp.process.processid.lower() == 'translatetiles':
            
            self._SetConstantsAndDicts()
            
            self._TranslateTilesRegions()
            
        else:
            
            exitstr = 'Exiting, processid %(p)s missing in ProcessGDAL' %{'p':self.pp.process.processid}
            
            exit(exitstr)
       
    def _ReprojectRegion(self):
        '''
        '''
        self.scriptFD = {}; self.scriptFPND = {}
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            for comp in self.pp.dstCompD:
                
                scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[comp].source,'region','reprojectscript')
                
                if not os.path.exists(scriptFP):
                    
                    os.makedirs(scriptFP)
                    
                scriptFN = 'reproject_%(comp)s-%(today)s.sh' %{'comp':comp, 'today':today}
                
                scriptFPN = os.path.join(scriptFP,scriptFN)
                
                self.scriptFPND[comp] = scriptFPN
                
                self.scriptFD[comp] = open(scriptFPN,'w')
                
                writeln = '# Script created by Kartturs GeoImagine Framework for reprojection %s, created %s\n' %(comp, today)
        
                self.scriptFD[comp].write(writeln)
                      
        for srcLocus in self.pp.srcLayerD:
                        
            for datum in self.pp.srcLayerD[srcLocus]:
                          
                for comp in self.pp.srcLayerD[srcLocus][datum]:
                                        
                    if not os.path.exists(self.pp.srcLayerD[srcLocus][datum][comp].FPN):
                    
                        exitstr = 'EXITING - region layer for tiling missing\n    %s' %(self.pp.srcLayerD[srcLocus][datum][comp].FPN)
            
                        exit(exitstr)
                        
                    # Loop over all destination tiles with the same comp and datum
                    
                    for locus in self.pp.dstLayerD:
                        
                        dstLayer = self.pp.dstLayerD[locus][datum][comp]
                        
                        if not dstLayer._Exists(): 
                        
                            if self.verbose > 1:
                            
                                infostr = '            reproject region: %s' %( dstLayer.FPN )
                                
                                print (infostr)
                                
                            if self.pp.srcCompD[comp].ext.lower() in ['.shp','.json','.geojson']:
                                
                                self._ReProjectVector(locus, datum, comp, self.pp.srcLayerD[srcLocus][datum][comp].FPN, dstLayer)
                            
                            else:
                                
                                self._ReProjectRaster(locus, datum, comp, self.pp.srcLayerD[srcLocus][datum][comp].FPN, dstLayer)
 
 
    def _TranslateTilesRegions(self):
        '''
        '''
        self.scriptFD = {}; self.scriptFPND = {}
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            for comp in self.pp.dstCompD:
                
                scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[comp].source,'region','reprojectscript')
                
                if not os.path.exists(scriptFP):
                    
                    os.makedirs(scriptFP)
                    
                scriptFN = 'translate_%(comp)s-%(today)s.sh' %{'comp':comp, 'today':today}
                
                scriptFPN = os.path.join(scriptFP,scriptFN)
                
                self.scriptFPND[comp] = scriptFPN
                
                self.scriptFD[comp] = open(scriptFPN,'w')
                
                writeln = '# Script created by Kartturs GeoImagine Framework for translating %s, created %s\n' %(comp, today)
        
                self.scriptFD[comp].write(writeln)
                      
        for locus in self.pp.srcLayerD:
                        
            for datum in self.pp.srcLayerD[locus]:
                          
                for comp in self.pp.srcLayerD[locus][datum]:
                                        
                    if not os.path.exists(self.pp.srcLayerD[locus][datum][comp].FPN):
                    
                        exitstr = 'EXITING - tile missing\n    %s' %(self.pp.srcLayerD[locus][datum][comp].FPN)
            
                        exit(exitstr)
                        
                    # Loop over all destination tiles with the same comp and datum
                    
                    for locus in self.pp.dstLayerD:
                        
                        dstLayer = self.pp.dstLayerD[locus][datum][comp]
                        
                        if not dstLayer._Exists(): 
                        
                            if self.verbose > 1:
                            
                                infostr = '            translate tile: %s' %( dstLayer.FPN )
                                
                                print (infostr)

                                self._GdalTranslate(dstLayer.FPN, self.pp.srcLayerD[locus][datum][comp].FPN)

 
    def _TranslateRegionOld(self):
        '''
        '''
        self.scriptFD = {}; self.scriptFPND = {}
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            for comp in self.pp.dstCompD:
                
                scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[comp].source,'region','reprojectscript')
                
                if not os.path.exists(scriptFP):
                    
                    os.makedirs(scriptFP)
                    
                scriptFN = 'reproject_%(comp)s-%(today)s.sh' %{'comp':comp, 'today':today}
                
                scriptFPN = os.path.join(scriptFP,scriptFN)
                
                self.scriptFPND[comp] = scriptFPN
                
                self.scriptFD[comp] = open(scriptFPN,'w')
                
                writeln = '# Script created by Kartturs GeoImagine Framework for reprojection %s, created %s\n' %(comp, today)
        
                self.scriptFD[comp].write(writeln)
                      
        for locus in self.pp.srcLayerD:
                        
            for datum in self.pp.srcLayerD[locus]:
                          
                for comp in self.pp.srcLayerD[locus][datum]:
                                        
                    if not os.path.exists(self.pp.srcLayerD[locus][datum][comp].FPN):
                    
                        exitstr = 'EXITING - region layer for tiling missing\n    %s' %(self.pp.srcLayerD[locus][datum][comp].FPN)
            
                        exit(exitstr)
                        
                    # Loop over all destination tiles with the same comp and datum
                    
                    for locus in self.pp.dstLayerD:
                        
                        dstLayer = self.pp.dstLayerD[locus][datum][comp]
                        
                        if not dstLayer._Exists(): 
                        
                            if self.verbose > 1:
                            
                                infostr = '            translate region: %s' %( dstLayer.FPN )
                                
                                print (infostr)
                                
                                
                                
                                self._GdalTranslate(dstLayer.FPN, self.pp.srcLayerD[locus][datum][comp].FPN)


    def _TileRegion(self):
        '''
        '''
        
        self.scriptFD = {}; self.scriptFPND = {}
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            for comp in self.pp.dstCompD:
                
                scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[comp].source,'tiles','tilingscript')
                
                if not os.path.exists(scriptFP):
                    
                    os.makedirs(scriptFP)
                    
                scriptFN = 'tiling_%(comp)s-%(today)s.sh' %{'comp':comp, 'today':today}
                
                scriptFPN = os.path.join(scriptFP,scriptFN)
                
                self.scriptFPND[comp] = scriptFPN
                
                self.scriptFD[comp] = open(scriptFPN,'w')
                
                writeln = '# Script created by Kartturs GeoImagine Framework for tiling %s, created %s\n' %(comp, today)
        
                self.scriptFD[comp].write(writeln)
                      
        for srcLocus in self.pp.srcLayerD:
                        
            for datum in self.pp.srcLayerD[srcLocus]:
                          
                for comp in self.pp.srcLayerD[srcLocus][datum]:
                                        
                    if not os.path.exists(self.pp.srcLayerD[srcLocus][datum][comp].FPN):
                    
                        exitstr = 'EXITING - region layer for tiling missing\n    %s' %(self.pp.srcLayerD[srcLocus][datum][comp].FPN)
            
                        exit(exitstr)
                        
                    # Loop over all destination tiles with the same comp and datum
                    
                    for locus in self.pp.dstLayerD:
                        
                        dstLayer = self.pp.dstLayerD[locus][datum][comp]
                        
                        if not dstLayer._Exists(): 
                        
                            if self.verbose > 1:
                            
                                infostr = '            creating tile: %s' %( dstLayer.FPN )
                                
                                print (infostr)
                            
                            self._MakeTile(locus, datum, comp, self.pp.srcLayerD[srcLocus][datum][comp].FPN, dstLayer)
                             
                        if os.path.exists(dstLayer.FPN):
                            
                            pass
                        
                            #self.session._InsertLayer(self.process.dstLayerD[locus][datum][comp],self.process.overwrite,self.process.delete)
        
        if self.pp.process.parameters.asscript:
            
            for comp in self.pp.dstCompD:
                                
                infostr = 'To tile the data please run the script file',self.scriptFPND[comp]
                
                print (infostr)
                
                self.scriptFD[comp].close()
                      
    def _ReProjectVector(self, locus,  datum, comp, srcLayerFPN, dstLayer): 
        
        
        cmd = 'ogr2ogr -overwrite '
                        
        cmd += '-t_srs EPSG:%(epsg)s ' %{'epsg':self.pp.procsys.dstepsg}
            
        cmd += '%(dst)s %(src)s' %{'dst':dstLayer.FPN, 'src':srcLayerFPN}       
             
        if self.pp.process.parameters.asscript:
            
            self.scriptFD[comp].write(cmd)
            
        else:

            print (cmd)
            
            os.system(cmd)
                
    def _ReProjectRaster(self, locus,  datum, comp, srcLayerFPN, dstLayer):
        '''
        '''
        
        pass
                                                 
    def _MakeTile(self, locus,  datum, comp, srcLayerFPN, dstLayer):
        '''
        '''
        
        # Get the destination projection and extent for the new tile
        if self.pp.procsys.dstsystem == 'modis':
            queryD = {'hvtile':locus}
            paramL =['minxsin','minysin','maxxsin','maxysin','ullat','ullon','lrlat','lrlon','urlat','urlon','lllat','lllon']
            dstCoords = self.session._SelectSingleTileCoords(queryD, paramL)
            coordsD = dict(zip(paramL,dstCoords))

        elif self.pp.procsys.dstsystem[0:4] == 'ease':
                                               
            queryD = {'xytile':locus} 
            
            paramL =['minxease','minyease','maxxease','maxyease','ullat','ullon','lrlat','lrlon','urlat','urlon','lllat','lllon']

            dstCoords = self.session._SelectSingleTileCoords(queryD, paramL, self.pp.procsys.dstsystem)
                         
        paramL = ['minx','miny','maxx','maxy','ullat','ullon','lrlat','lrlon','urlat','urlon','lllat','lllon']

        coordsD = dict(zip(paramL,dstCoords))
                        
        if self.pp.process.parameters.asscript:
            
            
            
            '''
            GDALwarp = GDALstuff(self.process.srcLayerD[self.srcLocus][datum][comp].FPN, self.process.dstLayerD[locus][datum][comp].FPN, self.process.params)
            GDALwarp.SetClipBox(coordsD['minxsin'], coordsD['maxysin'], coordsD['maxxsin'], coordsD['minysin'])
            GDALwarp.SetTargetProj(self.process.params.epsg)
            gdalcmd = GDALwarp.WarpRaster(self.process.params.asscript)
            '''
            
            cmd = 'GDALwarp -t_srs %(t_srs)s '  %{'t_srs': self.dstSRS} 
            cmd += '-te %(xmin)f %(ymin)f %(xmax)f %(ymax)f ' %{'xmin':coordsD['minx'],'ymin':coordsD['miny'], 'xmax':coordsD['maxx'], 'ymax':coordsD['maxy']}
            cmd += '-tr %(xres)f %(yres)f ' %{'xres':self.pp.process.parameters.tr_xres, 'yres':self.pp.process.parameters.tr_yres}
            cmd += '-r %(resample)s ' %{'resample': self.pp.process.parameters.resample }
                        
            if self.pp.process.parameters.errorthreshold > 0:
                
                cmd += '-et %(et)f ' %{'et':self.pp.process.parameters.errorthreshold}
                
            if not self.pp.process.parameters.celltype == 'auto':
                
                cmd += '-of %(of)f ' %{'of':self.pp.process.parameters.celltype}
                
            if hasattr(self.pp.srcCompD[comp], 'cellnull'):
                
                cmd += '-srcnodata %(null)f -dstnodata %(null)f ' %{'null':self.pp.srcCompD[comp].cellnull}
                
            cmd += '%(src)s %(dst)s\n' %{'src':srcLayerFPN ,'dst':dstLayer.FPN}    
                
            self.scriptFD[comp].write(cmd)
            
        else:
            # see https://gdal.org/python/osgeo.gdal-module.html#WarpOptions
            # https://gis.stackexchange.com/questions/257257/how-to-use-gdal-warp-cutline-option
            '''
            ds = gdal.Open(srcLayerFPN)
            
            gdal.WarpOptions(options, 
            format, 
            outputBounds, 
            outputBoundsSRS, 
            xRes, 
            yRes, 
            targetAlignedPixels, 
            width, 
            height, 
            srcSRS, 
            dstSRS, 
            coordinateOperation, srcAlpha, dstAlpha, warpOptions, errorThreshold, warpMemoryLimit, creationOptions, outputType, workingType, resampleAlg, srcNodata, dstNodata, multithread, tps, rpc, geoloc, polynomialOrder, transformerOptions, cutlineDSName, cutlineLayer, cutlineWhere, cutlineSQL, cutlineBlend, cropToCutline, copyMetadata, metadataConflictValue, setColorInterpretation, overviewLevel, callback, callback_data)
            


            'outputBounds'= [ coordsD['minx'], coordsD['miny'], coordsD['maxx'], coordsD['maxy'] ]
            'xRes': self.pp.process.parameters.tr_xres
            'yRes': self.pp.process.parameters.tr_yres
            'dstSRS': self.pp.procsys.dstepsg
            'resampleAlg': self.pp.process.parameters.tr_xres
            '''
            '''
            options to test:
            - geoloc = raises ERROR 1: Unable to compute a GEOLOC_ARRAY based transformation between pixel/line and georeferenced coordinates 

            '''
            '''
            options = gdal.WarpOptions(format = 'GTiff', 
                    dstSRS= self.dstSRS,
                    xRes= self.pp.process.parameters.tr_xres,
                    yRes= self.pp.process.parameters.tr_yres,
                    outputBounds= (coordsD['minx'], coordsD['miny'], coordsD['maxx'], coordsD['maxy']),
                    resampleAlg= self.resampleD[self.pp.process.parameters.resample.lower()[0:4]], 
                    errorThreshold = 0.25)
            '''
            
            kwargs = {'format': 'GTiff', 
                    'dstSRS': self.dstSRS,
                    'xRes': self.pp.process.parameters.tr_xres,
                    'yRes': self.pp.process.parameters.tr_yres,
                    'outputBounds': (coordsD['minx'], coordsD['miny'], coordsD['maxx'], coordsD['maxy']),
                    'resampleAlg': self.resampleD[self.pp.process.parameters.resample.lower()[0:4]]}
            
            if hasattr(self.pp.process.parameters, 'cellnull'):
            
                kwargs['dstNodata'] = self.pp.process.parameters.cellnull

            if not self.pp.process.parameters.celltype == 'auto':
                
                outputType = self.cellTypeD[self.pp.process.parameters.celltype]
                
                kwargs['outputType'] = outputType
                
            if self.pp.process.parameters.errorthreshold > 0:
                                
                kwargs['errorTHreshold'] = self.pp.process.parameters.errorthreshold
                  
            ds = gdal.Warp(dstLayer.FPN, 
                           srcLayerFPN,
                           **kwargs)
            '''
            ds = gdal.Warp(dstLayer.FPN, 
                           srcLayerFPN,
                           outputBounds = (coordsD['minx'], coordsD['miny'], coordsD['maxx'], coordsD['maxy']),
                           dstSRS = self.dstSRS,
                           resampleAlg = self.resampleD[self.pp.process.parameters.resample.lower()[0:4]],
                           xRes = self.pp.process.parameters.tr_xres,
                           yRes = self.pp.process.parameters.tr_yres )
            ''' 
            ds = None
      
    def _MosaicTiles(self):
        '''
        '''
        
        for comp in self.pp.srcCompL:
            
            for datum in self.pp.srcPeriod.datumL:
                
                # Check if destination file is in place

                if self.pp.dstLayerD[self.pp.defregion][datum][comp]._Exists() and not self.pp.process.overwrite:
                    
                    continue
                
                # Check for missing tiles if the processes does not accept missing 
                if not self.pp.process.acceptmissing:
                        
                    missingL = []
                    
                    # Check that all input tiles are included
                    for tile in self.pp.srcLayerDateNonExistD:
                        
                        if datum in self.pp.srcLayerDateNonExistD[tile][comp]:
                            
                            missingL.append(tile)
                            
                    if len(missingL) > 0:
                        
                        printstr = 'Input tiles missing for mosaicking\n'
                        printstr += '    region: %s, datum: %s\n' %(self.pp.defregion, datum)
                        for item in missingL:
                            printstr += '        tile: %s' %(item)
                            
                        exit(printstr)
                 
                tileL = []
                    
                for tile in self.pp.srcLayerDateNonExistD:
                    
                    if datum in self.pp.srcLayerDateExistD[tile][comp]:
                        
                        tileL.append(self.pp.srcLayerD[tile][datum][comp].FPN)
                                                
                if len(tileL) == 0:  
                    
                    exitstr = 'EXITING - no tiles for creating mosaic'
                    
                    exit (exitstr)  
                     
                xRes = self.pp.process.parameters.tr_xres
                yRes =  self.pp.process.parameters.tr_yres
                resampleAlg = self.resampleD[self.pp.process.parameters.resample.lower()[0:4]]
                    
                vrtFPN = MakeMosaic(tileL, 
                                     self.pp.dstLayerD[self.pp.defregion][datum][comp].FPN,
                                     xRes,yRes,
                                     resampleAlg) 
                
                paramL = ['minx', 'maxy', 'maxx', 'miny']
                
                queryD = {'regionid': self.pp.defregion}
                
                projwin = self.session._SingleSearch(queryD, paramL, 'system', 'regions')
                
                self._Translate(vrtFPN, self.pp.dstLayerD[self.pp.defregion][datum][comp].FPN, 'GTiff', projwin)

                if self.pp.process.parameters.fillnodata:
                         
                    self._FillNodata(self.pp.dstLayerD[self.pp.defregion][datum][comp].FPN,
                                    self.pp.process.parameters.fillmaxdist,
                                    self.pp.process.parameters.fillsmooth
                                     )

    def _Translate(self,srcFPN, dstFPN, gdalformat='GTiff', projwin=False, scaleparams=False, exponents=False, xyres=False, resample=False, size=False):
        '''
        '''
        # see https://gdal.org/python/osgeo.gdal-module.html#TranslateOptions
        
        kwargs = {'format': gdalformat}
        
        if projwin:
            
            #projwin should be in the form [ulx, uly, lrx, lry]
            
            kwargs['projWin'] = projwin
            
        if scaleparams:
            
            #scaleParams should be in the form [src_min,src_max,dst_min,dst_max]
            kwargs['scalePArams'] = scaleparams
            
            if exponents:
                
                # exponents only possible combined with scaleparams
                kwargs['exponents'] = exponents
                            
        if xyres:
            
            #xyres should be in the form (xres,yrs)
            kwargs['xres'] = xyres[0]; kwargs['xres'] = xyres[1]
            
        if resample:
            
            kwargs['resampleAlg'] = self.resampleD[resample[0:4]]
            
        if size:
            
            kwargs['width']= size[0]
            
            kwarg['height'] = size[1]
            

        ds = gdal.Open(srcFPN)
            
        ds = gdal.Translate(dstFPN, ds, 
                            **kwargs)
        
        ds = None
        
    def _FillNodata(self,srcFPN, maxSearchDist=10, smoothingIterations=0):
        '''
        '''
        # see https://gdal.org/python/osgeo.gdal-module.html#TranslateOptions
        
        if self.verbose > 1:
            
            infostr = '            Filling No data'
            
            print (infostr)
            
        ds = gdal.Open(srcFPN, GA_Update)
        
        gdal.FillNodata(targetBand = ds.GetRasterBand(1), maxSearchDist = maxSearchDist, maskBand = None, smoothingIterations = smoothingIterations)
 
        '''
        kwargs = {'maxSearchDist':maxSearchDist,
                  'smoothingIterations':smoothingIterations}
        
        
            

        ds = gdal.Open(srcFPN)
            
        ds = gdal.FillNodata(dstFPN, ds, 
                            **kwargs)
        '''
        ds = None
                         
    def _SetConstantsAndDicts(self):
        '''
        '''
        
        self.dstSRS = 'EPSG:%(epsg)d' %{'epsg':self.pp.procsys.dstepsg}
            
        self.resampleD = {'near':gdal.GRA_NearestNeighbour, 'bili':gdal.GRA_Bilinear,
                     'cubi':gdal.GRA_Cubic,'aver':gdal.GRA_Average}
        
        self.cellTypeD = {'byte':'gdalconst.GDT_Byte','etc':'gdalconst.GDT_Byte'}
        
        self.gdalFormatD = {'.tif':'GTiff', '.vrt':'VRT'}
        
        
        
        prefix = {'slope':'slope','aspect':'aspect','hillshade':'hillshade','color-relief':'color-relief',
                  'TRI':'tri','TPI':'tpi','roughness':'roughness'}
        
        layerid = {'slope':'slope','aspect':'aspect','hillshade':'hillshade','color-relief':'color-relief',
                  'TRI':'tri','TPI':'tpi','roughness':'roughness'}
        
        suffix = {'slope':'slope','aspect':'aspect','hillshade':'hillshade','color-relief':'color-relief',
                  'TRI':'tri','TPI':'tpi','roughness':'roughness'}
        

        
        