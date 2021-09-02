
import numpy as np
from osgeo import gdal, osr
from pathlib import Path
import re


def get_tile_coords(f):
    '''returns tile coords part of Hansen GFC file name'''
    return re.search('[0-9]*[NS]_[0-9]*[EW]', f)[0]


def get_lossname(cover_file):
    '''returns name of the loss year file associated with a particular cover file'''
    tile_coords = re.search('[0-9]*[NS]_[0-9]*[EW]', cover_file)[0]
    return('Hansen_GFC-2020-v1.8_lossyear_' + get_tile_coords(cover_file) + '.tif')


def recode_cover_file(
    arguments: dict
):
    cover_file = arguments['cover_file']
    out_folder = arguments['out_folder']
    xres = arguments['xres']
    yres = arguments['yres']
    
    cover_ds = gdal.Open(str(cover_file))
    cover_band = cover_ds.GetRasterBand(1)

    cover = cover_band.ReadAsArray()
    cover = cover > 30

    loss_ds = gdal.Open(str(cover_file.parent / get_lossname(str(cover_file.name))))
    loss_band = loss_ds.GetRasterBand(1)
    loss = loss_band.ReadAsArray()

    driver = gdal.GetDriverByName("MEM")
    mem_ds = driver.Create(
        '',
        cover_band.XSize,
        cover_band.YSize,
        21,
        gdal.GDT_Byte
    )
    src_gt = cover_ds.GetGeoTransform()
    mem_ds.SetGeoTransform(src_gt)
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromWkt(cover_ds.GetProjectionRef())
    mem_ds.SetProjection(dst_srs.ExportToWkt())
    # Threshold forest to be greater than 30% cover
    for n in range(21):
        mem_ds.GetRasterBand(n + 1).WriteArray(cover * np.logical_or(loss == 0, loss > n) * 100)

    out_file = str(
        out_folder / Path(
            'Hansen_GFC-2020-v1.8_250m_coverbyyear_' + get_tile_coords(str(cover_file.name)) + '.tif'
        )
    )
    gdal.Warp(
        out_file,
        mem_ds,
        xRes=xres,
        yRes=yres,
        resampleAlg='average',
        options=['COMPRESS=LZW'],
        multithread=True,
        dstSRS="epsg:4326",
        outputType=gdal.GDT_Byte
    )
